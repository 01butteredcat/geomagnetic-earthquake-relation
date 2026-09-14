"""Parallel ETL: reduce every <station><date>dsec.sec file to (a) one row of
daily summary features and (b) a 1-minute-resampled series, without ever
holding more than one file's worth of 1Hz data in memory at a time.

Usage: build_daily_features.py --group G10

Outputs (under the group's own data/interim/<group>/):
  daily_features.csv
  minute_series_<station>.parquet   (one file per station)
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import list_day_refs, load_group_config  # noqa: E402
from parser import DayFileRef, _parse_day_file_from_text, parse_day_file  # noqa: E402

N_WORKERS = 16

# Second-to-second jump threshold (nT) beyond which a sample is treated as an
# instrument/telemetry glitch rather than real geomagnetic variation. Real
# secular + Sq + storm-time field changes are at most a few nT/s even during
# severe storms; this threshold is set far above that (confirmed empirically:
# twu shows post-quake glitches with jumps of 1e5 nT, ttn shows chronic
# smaller spikes up to ~7700 nT — both are instrumental, not geophysical).
SPIKE_THRESHOLD_NT = 300.0


def _despike(series: pd.Series) -> tuple[pd.Series, pd.Series, float]:
    """Return (cleaned series with newly-detected spikes set to NaN,
    boolean spike_mask that is True ONLY for glitch samples (never for
    pre-existing sentinel-NaN gaps, since diff() against NaN is NaN and
    NaN > threshold is False), max |jump| between two valid samples)."""
    d = series.diff().abs()
    max_jump = float(d.max()) if d.notna().any() else 0.0
    spike_mask = d > SPIKE_THRESHOLD_NT  # NaN comparisons are False, so pre-existing gaps are excluded
    cleaned = series.copy()
    cleaned[spike_mask] = np.nan
    return cleaned, spike_mask, max_jump


def _list_files(gdms_dir, stations: dict):
    """Split the group's day files into loose (plain .sec/.sec.gz) and
    tgz-sourced refs. Kept separate rather than one merged list because the
    two need completely different treatment in main(): loose files are
    cheap for a worker to open itself, but .tgz members are not (gzip has
    no random access -- see module docstring / build plan), so those are
    extracted once in the main process instead of per-worker."""
    refs = [r for r in list_day_refs(gdms_dir) if r.station in stations]
    loose = [r for r in refs if r.kind == "loose"]
    tgz = [r for r in refs if r.kind == "tgz"]
    return loose, tgz


def _features_from_df(df, station, date_str):
    n_total = len(df)
    if "F" in df.columns:  # scalar-only file, per its own header (see parser.is_scalar_only)
        f_raw = df["F"]
        f, spike_mask, max_jump = _despike(f_raw)
        n_valid = int(f.notna().sum())
        row = {
            "station": station,
            "date": date_str,
            "n_total": n_total,
            "n_valid": n_valid,
            "pct_missing": round(1 - n_valid / n_total, 4),
            "n_spikes": int(spike_mask.sum()),
            "max_abs_jump": round(max_jump, 2),
            "F_mean": float(f.mean()) if n_valid else np.nan,
            "F_std": float(f.std()) if n_valid else np.nan,
            "F_range": float(f.max() - f.min()) if n_valid else np.nan,
        }
        minute = df[["F"]].resample("1min").mean()
    else:
        x_raw, y_raw, z_raw = df["X"], df["Y"], df["Z"]
        x, spike_x, max_jump_x = _despike(x_raw)
        y, spike_y, max_jump_y = _despike(y_raw)
        z, spike_z, max_jump_z = _despike(z_raw)
        # a spike in any component invalidates that second's vector for H/D;
        # this is on top of (not instead of) pre-existing sentinel-NaN gaps
        spike_any = spike_x | spike_y | spike_z
        x, y, z = x.mask(spike_any), y.mask(spike_any), z.mask(spike_any)
        h = np.sqrt(x**2 + y**2)
        d = np.degrees(np.arctan2(y, x))
        n_valid = int(x.notna().sum())
        row = {
            "station": station,
            "date": date_str,
            "n_total": n_total,
            "n_valid": n_valid,
            "pct_missing": round(1 - n_valid / n_total, 4),
            "n_spikes": int(spike_any.sum()),
            "max_abs_jump": round(max(max_jump_x, max_jump_y, max_jump_z), 2),
            "H_mean": float(h.mean()) if n_valid else np.nan,
            "H_std": float(h.std()) if n_valid else np.nan,
            "H_range": float(h.max() - h.min()) if n_valid else np.nan,
            "Z_mean": float(z.mean()) if n_valid else np.nan,
            "Z_std": float(z.std()) if n_valid else np.nan,
            "Z_range": float(z.max() - z.min()) if n_valid else np.nan,
            "D_mean": float(d.mean()) if n_valid else np.nan,
        }
        minute = pd.DataFrame({"X": x, "Y": y, "Z": z, "H": h}).resample("1min").mean()

    return row, (station, minute)


def _process_one_ref(ref: DayFileRef):
    """Loose-file path: worker does its own I/O + parse, same as before
    this module gained .tgz support (no gzip random-access penalty for
    loose files, so there's no reason to change this)."""
    try:
        df = parse_day_file(ref, ref.station)
    except Exception as exc:  # noqa: BLE001
        return {"station": ref.station, "date": ref.date_str, "error": str(exc)}, None
    return _features_from_df(df, ref.station, ref.date_str)


def _process_one_text(text: str, ref: DayFileRef):
    """.tgz-sourced path: the member's text was already extracted once in
    the main process (see main()'s tgz dispatch loop), so this worker does
    pure CPU-bound parse/despike/resample and never touches tarfile itself."""
    try:
        df = _parse_day_file_from_text(text, ref, ref.station)
    except Exception as exc:  # noqa: BLE001
        return {"station": ref.station, "date": ref.date_str, "error": str(exc)}, None
    return _features_from_df(df, ref.station, ref.date_str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    loose_refs, tgz_refs = _list_files(cfg.gdms_dir, cfg.stations)
    n_files = len(loose_refs) + len(tgz_refs)
    print(
        f"[{args.group}] {n_files} files to process "
        f"({len(loose_refs)} loose, {len(tgz_refs)} from .tgz)",
        file=sys.stderr,
    )

    # Group the .tgz-sourced refs by archive so each archive is opened and
    # walked sequentially exactly once (see module docstring: gzip has no
    # random access, so per-member cold-opens -- e.g. one per worker -- are
    # ~600x slower than a single sequential pass, measured on real GDMS
    # batch downloads).
    by_tgz: dict[Path, dict[str, DayFileRef]] = {}
    for r in tgz_refs:
        by_tgz.setdefault(r.source_path, {})[r.member] = r

    rows = []
    errors = []
    minute_chunks: dict[str, list[pd.DataFrame]] = {s: [] for s in cfg.stations}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = [pool.submit(_process_one_ref, r) for r in loose_refs]
        for tgz_path, wanted in by_tgz.items():
            with tarfile.open(tgz_path, "r:gz") as tf:
                for member in tf:  # single sequential scan over this archive
                    if member.name not in wanted:
                        continue
                    ref = wanted[member.name]
                    extracted = tf.extractfile(member)
                    text = extracted.read().decode("ascii", errors="strict")
                    futures.append(pool.submit(_process_one_text, text, ref))
        done = 0
        for fut in as_completed(futures):
            row, minute_pair = fut.result()
            done += 1
            if "error" in row:
                errors.append(row)
                print(f"ERROR {row['station']} {row['date']}: {row['error']}", file=sys.stderr)
                continue
            rows.append(row)
            station, minute = minute_pair
            minute_chunks[station].append(minute)
            if done % 200 == 0:
                print(f"  {done}/{n_files}", file=sys.stderr)

    daily = pd.DataFrame(rows).sort_values(["station", "date"])
    daily.to_csv(cfg.interim_dir / "daily_features.csv", index=False)
    print(f"wrote daily_features.csv ({len(daily)} rows)", file=sys.stderr)

    if errors:
        pd.DataFrame(errors).to_csv(cfg.interim_dir / "parse_errors.csv", index=False)
        print(f"wrote parse_errors.csv ({len(errors)} rows)", file=sys.stderr)

    for station, chunks in minute_chunks.items():
        if not chunks:
            continue
        merged = pd.concat(chunks).sort_index()
        out_path = cfg.interim_dir / f"minute_series_{station}.parquet"
        merged.to_parquet(out_path)
        print(f"wrote {out_path.name} ({len(merged)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
