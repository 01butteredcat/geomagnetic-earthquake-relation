"""Shared loader for the expanded earthquake catalog (`fetch_earthquake_
catalog.py`'s output) + events.py's 20 hand-curated events, merged into one
event list without double-counting -- used by both
`superposed_epoch_analysis.py` and `backtest_rule.py` so their event
population definitions can't silently drift apart."""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from events import get_group


def load_extended_events(catalog_path: Path, group_ids: tuple[str, ...]) -> list[dict]:
    """Union of events.py's events (for group_ids) and the extended catalog's
    declustered/not-already-known rows. Each item: group, date (pd.Timestamp
    at day resolution), mag, source ("events.py", "USGS", or "CWA_GDMS" --
    the latter two per-row from fetch_earthquake_catalog.py's own "source"
    column; "usgs_catalog" as a fallback for catalogs fetched before that
    column existed)."""
    events: list[dict] = []
    for group_id in group_ids:
        for ev in get_group(group_id).events:
            events.append({
                "group": group_id,
                "date": pd.Timestamp(ev.time_utc.split(" ")[0]),
                "mag": ev.magnitude,
                "source": "events.py",
            })
    if catalog_path.exists():
        with catalog_path.open(newline="") as f:
            for row in csv.DictReader(f):
                if row["group"] not in group_ids:
                    continue
                if row["declustered"] != "True" or row["is_known_event"] == "True":
                    continue
                events.append({
                    "group": row["group"],
                    "date": pd.Timestamp(row["time_utc"].split(" ")[0]),
                    "mag": float(row["mag"]),
                    "source": row.get("source") or "usgs_catalog",
                })
    return events
