"""IAGA-2002 parser for GDMSdata 1-second geomagnetic files.

Handles the two station conventions found across the full G1-G13 dataset:
  - vector stations: real data in X/Y/Z, F is the structural placeholder
    88888.00 (not reported).
  - scalar-only stations: real data in F, X/Y/Z always 88888.00.

Which convention applies is NOT a fixed property of a station code -- most
stations start scalar-only and are upgraded to vector partway through the
dataset's history (the transition always lands on a group-era boundary, per
the G1-G13 station registry survey; `ttn` is the one long-lived station that
never upgrades). Each file's own header states which applies via the
`Reported` line (`F` or `XYZF`), so `parse_day_file` reads that per file
rather than consulting a static registry -- see `parse_header`/`is_scalar_only`.

Both conventions also use 99999.00 as a distinct sentinel meaning genuine
data outage (as opposed to 88888.00's "channel not reported" meaning), which
is mapped to NaN in whichever column(s) are the "real" ones for that station.
"""
from __future__ import annotations

import gzip
import io
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

OUTAGE_SENTINEL = 99999.0
NOT_REPORTED_PLACEHOLDER = 88888.0

_HEADER_LINE_RE = re.compile(r"^DATE\s+TIME\s+DOY\s+")
_HEADER_FIELD_RE = re.compile(r"^\s*(.+?)\s{2,}(.*?)\s*\|\s*$")


@dataclass(frozen=True)
class DayFileRef:
    """Points at one station-day's worth of .sec data, regardless of whether
    it lives as a loose file or as a member inside a GDMS batch-download
    .tgz. Only str/Path fields -- safe to pickle across
    ProcessPoolExecutor process boundaries (unlike an open tarfile handle).
    `member` is None for loose files, the tar member name for .tgz sources.
    """

    station: str
    date_str: str  # YYYYMMDD
    source_path: Path  # loose: the .sec/.sec.gz itself; tgz: the .tgz file
    member: str | None = None

    @property
    def kind(self) -> str:
        return "tgz" if self.member is not None else "loose"

    @property
    def label(self) -> str:
        """For log/error messages, in place of just printing a bare path."""
        return f"{self.source_path}!{self.member}" if self.member else str(self.source_path)


def _as_ref(x: "str | Path | DayFileRef", station_code: str | None = None) -> DayFileRef:
    """Normalize a loose path (str/Path, as every caller historically passed)
    into a DayFileRef so the rest of this module only has one code path."""
    if isinstance(x, DayFileRef):
        return x
    path = Path(x)
    stem = path.name[:-3] if path.name.endswith(".sec.gz") else path.name
    return DayFileRef(station_code or stem[:3], stem[3:11], path, None)


def _read_text(ref: DayFileRef) -> str:
    """Read this ref's full file content into memory exactly once. A .sec
    file is ~6MB decompressed, small enough that reading it whole up front
    (instead of re-opening/re-decompressing per pass) is cheap for loose
    files and *necessary* for .tgz members: gzip doesn't support random
    access, so re-opening the same .tgz per read would be far more
    expensive than reading it once (see build_daily_features.py for the
    batch-extraction path that avoids doing this per file per worker)."""
    if ref.member is None:  # loose .sec or .sec.gz
        path = ref.source_path
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="ascii", errors="strict") as f:
                return f.read()
        with open(path, "r", encoding="ascii", errors="strict") as f:
            return f.read()
    # .tgz member
    with tarfile.open(ref.source_path, "r:gz") as tf:
        extracted = tf.extractfile(ref.member)
        if extracted is None:
            raise ValueError(f"{ref.label}: tar member is not a regular file")
        raw = extracted.read()
    return raw.decode("ascii", errors="strict")


def open_raw(ref: "str | Path | DayFileRef"):
    """Open a .sec/.sec.gz file or .tgz member for text reading, returning a
    context-manager-able, line-iterable object -- same interface as before
    (backward compatible with callers passing a bare str/Path, e.g.
    verify_pipeline.py's spot-check). Internally now always reads the whole
    file/member up front (see `_read_text`) rather than returning a live
    gzip/plain file handle."""
    return io.StringIO(_read_text(_as_ref(ref)))


def _find_data_start_from_text(text: str, label: str) -> int:
    """Return the 0-indexed line number of the first data row."""
    for i, line in enumerate(text.splitlines(keepends=True)):
        if _HEADER_LINE_RE.match(line):
            return i + 1
    raise ValueError(f"Could not find IAGA-2002 column header line in {label}")


def _parse_header_from_text(text: str, label: str) -> dict:
    fields: dict[str, str] = {}
    for line in text.splitlines(keepends=True):
        if _HEADER_LINE_RE.match(line):
            break
        m = _HEADER_FIELD_RE.match(line)
        if m:
            fields[m.group(1).strip()] = m.group(2).strip()

    if "Geodetic Latitude" not in fields or "Geodetic Longitude" not in fields:
        raise ValueError(f"{label}: could not find lat/lon in IAGA-2002 header")

    elevation = fields.get("Elevation", "")
    return {
        "station_name": fields.get("Station Name", ""),
        "lat": float(fields["Geodetic Latitude"]),
        "lon": float(fields["Geodetic Longitude"]),
        "elevation_m": float(elevation) if elevation else None,
        "reported": fields.get("Reported", ""),
        "source": fields.get("Source of Data", ""),
    }


def parse_header(ref: "str | Path | DayFileRef") -> dict:
    """Parse the IAGA-2002 header of one .sec file/tgz member into station
    metadata.

    Every file carries its own station name/lat/lon/elevation and a
    `Reported` field (`F` or `XYZF`) -- confirmed consistent across dates for
    a given station/era by spot-checking multiple files per station code
    across all 13 groups. This means station metadata and scalar-vs-vector
    classification can be derived directly from the data rather than
    hand-maintained in a separate lookup table.

    Returns dict with: station_name, lat, lon, elevation_m (None if the
    header field is blank -- seen for several older/retired stations),
    reported (raw string, e.g. "F" or "XYZF"), source.
    """
    ref = _as_ref(ref)
    return _parse_header_from_text(_read_text(ref), ref.label)


def is_scalar_only(ref: "str | Path | DayFileRef") -> bool:
    """True if this file's own header reports F only (no real X/Y/Z)."""
    return parse_header(ref)["reported"].strip().upper() == "F"


def _parse_day_file_from_text(text: str, ref: DayFileRef, station_code: str) -> pd.DataFrame:
    """Core parse logic, operating on already-read-into-memory text rather
    than re-reading from disk/tar. Exposed (not just an inline helper)
    because build_daily_features.py's batch .tgz path extracts member text
    once per file in the main process and hands it straight to worker
    processes, skipping the per-file I/O this function would otherwise do."""
    skiprows = _find_data_start_from_text(text, ref.label)
    scalar_only = _parse_header_from_text(text, ref.label)["reported"].strip().upper() == "F"

    df = pd.read_csv(
        io.StringIO(text),
        skiprows=skiprows,
        sep=r"\s+",
        header=None,
        names=["DATE", "TIME", "DOY", "X", "Y", "Z", "F"],
        dtype={"X": "float32", "Y": "float32", "Z": "float32", "F": "float32"},
        engine="c",
    )

    if len(df) != 86400:
        raise ValueError(
            f"{ref.label}: expected 86400 data rows, got {len(df)}"
        )

    ts = pd.to_datetime(df["DATE"] + " " + df["TIME"], format="%Y-%m-%d %H:%M:%S.%f")

    if scalar_only:
        real = df[["F"]].copy()
        real.loc[real["F"] == OUTAGE_SENTINEL, "F"] = pd.NA
    else:
        real = df[["X", "Y", "Z"]].copy()
        for col in ("X", "Y", "Z"):
            real.loc[real[col] == OUTAGE_SENTINEL, col] = pd.NA
            # Defensive: a stray not-reported placeholder in a vector column
            # would also be invalid data, though not expected per ground truth.
            real.loc[real[col] == NOT_REPORTED_PLACEHOLDER, col] = pd.NA

    real = real.astype("float32")
    real.index = ts
    real.index.name = "time"
    return real


def parse_day_file(ref: "str | Path | DayFileRef", station_code: str) -> pd.DataFrame:
    """Parse one <station><YYYYMMDD>dsec.sec file (loose or a .tgz member)
    into a DataFrame.

    Returns a DataFrame indexed by UTC-naive datetime (raw clock value as
    printed in the file; true UTC offset is resolved separately by
    timezone_check.py) with float32 columns depending on what this specific
    file's header reports (see `is_scalar_only`):
      - vector file: X, Y, Z
      - scalar-only file: F
    Missing/outage samples (sentinel 99999.00) are NaN. `station_code` is not
    used to determine this (kept for caller bookkeeping/error messages only)
    -- scalar-vs-vector status is read fresh from this file's own header.
    """
    ref = _as_ref(ref, station_code=station_code)
    return _parse_day_file_from_text(_read_text(ref), ref, station_code)


if __name__ == "__main__":
    import sys

    p = Path(sys.argv[1])
    station = p.name[:3]
    out = parse_day_file(p, station)
    print(out.describe())
    print(out.head())
    print("NaN count:\n", out.isna().sum())
