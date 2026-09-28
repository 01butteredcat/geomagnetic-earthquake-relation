"""ULF (Pc3/Pc4 band) polarization analysis for the near-vs-far station
comparison, per plan Section 5. Needs real X/Y/Z, so only runs when the
group's XYZ station pool is sufficient (see common.py's GroupConfig) --
scalar-only groups (e.g. G1) have no vector data to filter/Hilbert-transform
at all, so this step is skipped for them, not degraded.

For each station-day:
  1. Detrend X, Y, Z (subtract a 1-hour centered rolling mean).
  2. Zero-phase 4th-order Butterworth bandpass (scipy.signal.filtfilt) into
     Pc3 (10-45s period) and Pc4 (45-150s period) bands, per component.
  3. Hilbert envelope magnitude per component/band.
  4. Horizontal envelope = sqrt(envX^2 + envY^2); reduce to the confirmed
     local-night window (UTC hour in {17,18,19} == local 01:00-03:59) RMS.
  5. Z/H polarization ratio = RMS(envZ_night) / RMS(envH_night).

Near/far index (`ulf_near_far_index.csv`), station-normalized since 2026-09-28:
each station's nightly Z/H is first divided by that station's own median over
non-storm days, then the near and far pools take the median of those
normalized ratios, and diff = near - far. Before that the pools took the
median of the raw ratios over whichever stations had data that day, so a
station joining or leaving the pool shifted the index: G11's near station zbn
sat at Z/H 15-45 (the others 0.1-1.5) and stopped on 2024-12-23, which dropped
the near index ~20x on 12-24, right at the start of G11's pre-event window.
The raw-median version is kept as `<band>_*_zh_raw` for comparison, and
`<band>_n_near` / `<band>_n_far` record how many stations had data each day.
The reference median uses the whole series (one number per station), so it
sees future days too; that shifts a station's level, not its day-to-day shape.

Outputs (under the group's own data/interim/<group>/):
  ulf_daily.csv               -- full-window daily series
  ulf_near_far_index.csv      -- near-minus-far index (see above)
  ulf_spectrogram_<window>.json -- zoomed 1Hz spectrograms around the anchor event
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, hilbert, spectrogram

sys.path.insert(0, str(Path(__file__).parent))
from common import list_day_refs, load_group_config, resolve_day_ref  # noqa: E402
from parser import parse_day_file  # noqa: E402

FS = 1.0  # Hz
BANDS = {
    "pc3": (1 / 45, 1 / 10),   # 10-45s period
    "pc4": (1 / 150, 1 / 45),  # 45-150s period
}
NIGHT_HOURS_UTC = {17, 18, 19}
N_WORKERS = 10
ZOOM_WINDOW_DAYS = 12  # +/- this many days around the anchor event


def _bandpass(x: np.ndarray, low: float, high: float) -> np.ndarray:
    nyq = FS / 2
    b, a = butter(4, [low / nyq, high / nyq], btype="band")
    return filtfilt(b, a, x)


def _detrend(x: pd.Series) -> np.ndarray:
    trend = x.rolling(window=3601, center=True, min_periods=1).mean()
    return (x - trend).to_numpy()


def _process_day(ref):
    station, date_str = ref.station, ref.date_str
    try:
        df = parse_day_file(ref, station)
    except Exception as exc:  # noqa: BLE001
        return {"station": station, "date": date_str, "error": str(exc)}

    if "X" not in df.columns or df[["X", "Y", "Z"]].isna().to_numpy().any():
        # any gap breaks filtfilt's continuity assumption for this simple
        # pipeline; skip the day rather than risk filter artifacts near gaps
        # (also covers the defensive case of a scalar-only file slipping in)
        return {"station": station, "date": date_str, "error": "has_gap_or_scalar_skipped"}

    x = _detrend(df["X"])
    y = _detrend(df["Y"])
    z = _detrend(df["Z"])
    night_mask = df.index.hour.isin(NIGHT_HOURS_UTC)

    row = {"station": station, "date": date_str}
    for band_name, (low, high) in BANDS.items():
        xf, yf, zf = _bandpass(x, low, high), _bandpass(y, low, high), _bandpass(z, low, high)
        env_x, env_y, env_z = np.abs(hilbert(xf)), np.abs(hilbert(yf)), np.abs(hilbert(zf))
        env_h = np.sqrt(env_x**2 + env_y**2)

        rms_h_night = float(np.sqrt(np.mean(env_h[night_mask] ** 2)))
        rms_z_night = float(np.sqrt(np.mean(env_z[night_mask] ** 2)))
        row[f"{band_name}_H_night_rms"] = rms_h_night
        row[f"{band_name}_Z_night_rms"] = rms_z_night
        row[f"{band_name}_zh_ratio"] = rms_z_night / rms_h_night if rms_h_night > 1e-9 else np.nan

    return row


def build_daily_series(cfg, stations):
    refs = []
    for station in stations:
        refs.extend(list_day_refs(cfg.gdms_dir, station))
    print(f"[{cfg.group_id}] {len(refs)} station-days to process for ULF analysis", file=sys.stderr)

    rows, errors = [], []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = [pool.submit(_process_day, r) for r in refs]
        done = 0
        for fut in as_completed(futures):
            r = fut.result()
            done += 1
            if "error" in r:
                errors.append(r)
            else:
                rows.append(r)
            if done % 100 == 0:
                print(f"  {done}/{len(refs)}", file=sys.stderr)

    daily = pd.DataFrame(rows).sort_values(["station", "date"])
    daily.to_csv(cfg.interim_dir / "ulf_daily.csv", index=False)
    print(f"wrote ulf_daily.csv ({len(daily)} rows, {len(errors)} skipped)", file=sys.stderr)
    if errors:
        pd.DataFrame(errors).to_csv(cfg.interim_dir / "ulf_errors.csv", index=False)
    return daily


def build_zoom_spectrogram(cfg, station: str, start_date: str, end_date: str, label: str):
    dates = pd.date_range(start_date, end_date, freq="D").strftime("%Y%m%d")
    x_chunks = []
    for d in dates:
        ref = resolve_day_ref(cfg.gdms_dir, station, d)
        if ref is None:
            continue
        df = parse_day_file(ref, station)
        if "X" not in df.columns:
            continue
        if df["X"].isna().any():
            x_chunks.append(pd.Series(np.nan, index=df.index))
        else:
            x_chunks.append(df["X"])
    if not x_chunks:
        return
    x = pd.concat(x_chunks)
    x_detrended = _detrend(x.interpolate(limit=60))  # interpolate tiny gaps only

    f, t, sxx = spectrogram(x_detrended, fs=FS, nperseg=3600, noverlap=1800)
    band_mask = (f >= BANDS["pc4"][0]) & (f <= BANDS["pc3"][1])

    out = {
        "station": station,
        "start_date": start_date,
        "end_date": end_date,
        "freq_hz": f[band_mask].tolist(),
        "time_hours_from_start": (t / 3600).tolist(),
        "power_db": np.log10(sxx[band_mask, :] + 1e-6).tolist(),
    }
    out_path = cfg.interim_dir / f"ulf_spectrogram_{station}_{label}.json"
    out_path.write_text(json.dumps(out))
    print(f"wrote {out_path.name}", file=sys.stderr)


def near_far_table(daily: pd.DataFrame, near_stations, far_stations, storm_dates: set[str]) -> pd.DataFrame:
    """Near-minus-far differential of the nightly Z/H polarization ratio, per
    band: station-normalized (see module docstring) plus the old raw-median
    version as *_raw."""
    daily = daily.copy()
    quiet = ~daily["date"].astype(str).isin(storm_dates)
    for band in BANDS:
        col = f"{band}_zh_ratio"
        ref = daily[quiet].groupby("station")[col].median()
        ref = ref.where(ref > 1e-9)
        daily[f"{band}_zh_norm"] = daily[col] / daily["station"].map(ref)
    rows = []
    for date, g in daily.groupby("date"):
        row = {"date": date}
        near, far = g[g.station.isin(near_stations)], g[g.station.isin(far_stations)]
        for band in BANDS:
            for kind, col in (("", f"{band}_zh_norm"), ("_raw", f"{band}_zh_ratio")):
                nv, fv = near[col].dropna(), far[col].dropna()
                row[f"{band}_near_zh{kind}"] = nv.median() if len(nv) else np.nan
                row[f"{band}_far_zh{kind}"] = fv.median() if len(fv) else np.nan
                row[f"{band}_diff_zh{kind}"] = row[f"{band}_near_zh{kind}"] - row[f"{band}_far_zh{kind}"]
            row[f"{band}_n_near"] = int(near[f"{band}_zh_ratio"].notna().sum())
            row[f"{band}_n_far"] = int(far[f"{band}_zh_ratio"].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("date")


def build_near_far_differential(cfg, daily: pd.DataFrame, near_stations, far_stations):
    storm_path = cfg.interim_dir / "storm_days.csv"
    storm = set(pd.read_csv(storm_path, dtype={"date": str})["date"]) if storm_path.exists() else set()
    out = near_far_table(daily, near_stations, far_stations, storm)
    out.to_csv(cfg.interim_dir / "ulf_near_far_index.csv", index=False)
    print(f"wrote ulf_near_far_index.csv ({len(out)} rows)", file=sys.stderr)
    return out


def self_test() -> bool:
    """A high-level near station leaves the pool halfway: the raw-median near
    index jumps, the station-normalized one doesn't; raw columns keep the old
    formula."""
    rng = np.random.default_rng(0)
    dates = [f"202401{d:02d}" for d in range(1, 31)]
    levels = {"n1": 0.5, "n2": 20.0, "f1": 0.3, "f2": 0.4}
    rows = []
    for i, d in enumerate(dates):
        for st, lvl in levels.items():
            if st == "n2" and i >= 15:
                continue
            r = lvl * np.exp(rng.normal(0, 0.05))
            rows.append({"station": st, "date": d, "pc3_zh_ratio": r, "pc4_zh_ratio": r})
    daily = pd.DataFrame(rows)
    t = near_far_table(daily, ["n1", "n2"], ["f1", "f2"], set()).set_index("date")
    raw_jump = t.pc3_near_zh_raw.iloc[:15].median() / t.pc3_near_zh_raw.iloc[15:].median()
    norm_jump = t.pc3_near_zh.iloc[:15].median() / t.pc3_near_zh.iloc[15:].median()
    old = daily.groupby("date").apply(lambda g: g[g.station.isin(["n1", "n2"])].pc3_zh_ratio.median()
                                      - g[g.station.isin(["f1", "f2"])].pc3_zh_ratio.median())
    raw_same = np.allclose(t.pc3_diff_zh_raw.to_numpy(), old.reindex(t.index).to_numpy())
    ok = raw_jump > 5 and abs(norm_jump - 1) < 0.2 and raw_same and set(t.pc3_n_near) == {1, 2}
    print(f"[self-test] near index before/after dropout: raw {raw_jump:.1f}x, normalized {norm_jump:.2f}x; "
          f"raw columns = old formula: {raw_same}  {'PASS' if ok else 'FAIL'}", file=sys.stderr)
    return bool(ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--from-daily", action="store_true",
                    help="rebuild only ulf_near_far_index.csv from the existing ulf_daily.csv")
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not args.group:
        ap.error("--group is required")
    cfg = load_group_config(args.group)

    if not cfg.xyz_pool.sufficient:
        print(f"[{args.group}] XYZ pool insufficient ({len(cfg.xyz_pool.all_stations)} stations) "
              f"-- ULF polarization analysis needs real X/Y/Z and cannot run for this group", file=sys.stderr)
        return

    near, far = cfg.xyz_pool.near, cfg.xyz_pool.far
    stations = list(near) + list(far)

    if args.from_daily:
        daily = pd.read_csv(cfg.interim_dir / "ulf_daily.csv", dtype={"date": str})
        build_near_far_differential(cfg, daily, near, far)
        return
    daily = build_daily_series(cfg, stations)
    build_near_far_differential(cfg, daily, near, far)

    anchor_date = pd.to_datetime(cfg.anchor_event.time_utc.split(" ")[0])
    start = (anchor_date - pd.Timedelta(days=ZOOM_WINDOW_DAYS)).strftime("%Y-%m-%d")
    end = (anchor_date + pd.Timedelta(days=ZOOM_WINDOW_DAYS)).strftime("%Y-%m-%d")
    build_zoom_spectrogram(cfg, near[0], start, end, "eq_window")
    build_zoom_spectrogram(cfg, far[-1], start, end, "eq_window")


if __name__ == "__main__":
    main()
