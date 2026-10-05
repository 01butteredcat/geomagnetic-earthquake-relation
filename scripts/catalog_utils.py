"""擴充地震目錄（`fetch_earthquake_
catalog.py` 的輸出）＋ events.py 中 20 起人工整理事件的共用載入器，合併成一份
不重複計算的事件清單——
`superposed_epoch_analysis.py` 和 `backtest_rule.py` 都用它，這樣兩者的
事件母體定義就不會默默分歧。"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from events import get_group


def load_extended_events(catalog_path: Path, group_ids: tuple[str, ...], min_mag: float) -> list[dict]:
    """events.py 的事件（用於 group_ids）和延伸目錄中
    去叢集後、尚未登錄的列的聯集，兩者都限制在規模 >=
    min_mag——2026-09-23 補登後 events.py 也含 M5 事件，
    所以沒有這個篩選的話，每個級距都會默默把它們算進去（一個「M>=6.0」的
    疊加有 99 起事件，其中 67 起 M<6）。每一項：group、date（pd.Timestamp，
    日解析度）、mag、source（"events.py"、"USGS" 或 "CWA_GDMS"——
    後兩者逐列取自 fetch_earthquake_catalog.py 自己的 "source"
    欄位；在那個欄位存在之前抓的目錄，則以 "usgs_catalog"
    作為備用值）。"""
    events: list[dict] = []
    for group_id in group_ids:
        for ev in get_group(group_id).events:
            if ev.magnitude < min_mag:
                continue
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
                if float(row["mag"]) < min_mag:
                    continue
                events.append({
                    "group": row["group"],
                    "date": pd.Timestamp(row["time_utc"].split(" ")[0]),
                    "mag": float(row["mag"]),
                    "source": row.get("source") or "usgs_catalog",
                })
    return events
