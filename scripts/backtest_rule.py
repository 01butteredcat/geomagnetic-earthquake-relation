"""回測特定的 z <= -4.1 規則（`stat_utils.FIXED_RULE_THRESHOLD`，
也就是 `report_template.html` 對 G10 2024-03-30 低谷引用的數字），把它當成
真正的地震前兆警報規則，依教授的第三個
建議：在全部 8 個向量站足夠的組別中，標出每一個
pc3_diff_zh / pc4_diff_zh 整條序列中位數／MAD z-score 跨過 -4.1 的日子，
再拿去和擴充地震目錄（`fetch_earthquake_
catalog.py` + events.py）比對——有多少個被標記的日子真的在 N 天內
有符合條件的地震（命中），有多少個沒有（誤
報）？使用和
`cross_group_analysis.py::WINDOWS_DAYS` 相同的 7/14/30 天窗口慣例（三者並列報告，
而不是事後挑看起來最好的那個）。

也重用 `cross_group_analysis.py::sliding_baseline_rate`，讓
回測觀察到的命中率可以和那支腳本已經建立的**同一個**經驗
基準比較，而不是另外發明第二個
不一致的基準定義。

用法：
  backtest_rule.py --catalog data/external/extended_catalog_m5.5.csv --label m5.5 --min-mag 5.5
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from catalog_utils import load_extended_events  # noqa: E402
from events import folder_events  # noqa: E402
from common import PROJECT_DIR, load_group_config  # noqa: E402
from cross_group_analysis import sliding_baseline_rate  # noqa: E402
from stat_utils import FIXED_RULE_THRESHOLD, mad_zscore  # noqa: E402

ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23", "G24")
BANDS = ("pc3", "pc4")
WINDOWS_DAYS = [7, 14, 30]


def load_events_for_group(group_id: str, catalog_path: Path, min_mag: float) -> list[pd.Timestamp]:
    events = load_extended_events(catalog_path, (group_id,), min_mag)
    return sorted({e["date"] for e in events})


def flagged_dates(group_id: str, band: str) -> tuple[list[str], list[str]]:
    """回傳 (all_dates_in_group, flagged_dates)，皆為 YYYYMMDD 字串。"""
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    df = pd.read_csv(path, dtype={"date": str}).sort_values("date").reset_index(drop=True)
    col = f"{band}_diff_zh"
    mask = df[col].notna()
    dates = df.loc[mask, "date"].tolist()
    vals = df.loc[mask, col].to_numpy(dtype=float)
    if len(vals) < 20:
        return df["date"].tolist(), []
    z, _, mad = mad_zscore(vals)
    if mad <= 1e-9:
        return df["date"].tolist(), []
    flagged = [d for d, zi in zip(dates, z) if zi <= FIXED_RULE_THRESHOLD]
    return df["date"].tolist(), flagged


def precursor_window_dates(event_date: pd.Timestamp, window_days: int) -> set[str]:
    return {(event_date - pd.Timedelta(days=d)).strftime("%Y%m%d") for d in range(1, window_days + 1)}


def run_band(band: str, catalog_path: Path, min_mag: float) -> dict:
    per_group = {}
    for group_id in ULF_GROUPS:
        all_dates, flagged = flagged_dates(group_id, band)
        events = load_events_for_group(group_id, catalog_path, min_mag)
        per_group[group_id] = {"all_dates": all_dates, "flagged": flagged, "events": events}

    windows_out = {}
    for w in WINDOWS_DAYS:
        total_flagged = 0
        total_hit_days = 0
        total_false_alarm_days = 0
        total_precursor_days = 0
        total_days = 0
        n_events = 0
        n_events_with_hit = 0
        per_group_detail = {}
        for group_id, g in per_group.items():
            precursor_days: set[str] = set()
            for ev in g["events"]:
                precursor_days |= precursor_window_dates(ev, w)
            # 共用原始資料夾的組別看到的是同樣的日子：另一組的真實事件有它
            # 自己的前兆窗口，所以那些日子對**這一組**的指標沒有任何意義——
            # 把它們從命中／誤報／分母的計算中拿掉，而不是讓兄弟組
            # 地震之前的旗標在這裡被算成誤報。
            sibling_dates = {pd.Timestamp(e.time_utc.split(" ")[0]) for e in folder_events(group_id)} - set(g["events"])
            sibling_only_days: set[str] = set()
            for sd in sibling_dates:
                sibling_only_days |= precursor_window_dates(sd, w)
            sibling_only_days -= precursor_days
            eval_dates = set(g["all_dates"]) - sibling_only_days
            flagged_all = set(g["flagged"])
            flagged_set = flagged_all - sibling_only_days
            hit_days = flagged_set & precursor_days
            false_alarm_days = flagged_set - precursor_days

            events_with_hit = sum(
                1 for ev in g["events"] if precursor_window_dates(ev, w) & flagged_set
            )

            baseline_rate, n_windows = sliding_baseline_rate(
                sorted(g["all_dates"]), flagged_all, w, exclude_window=set()
            )

            per_group_detail[group_id] = {
                "n_days": len(eval_dates),
                "n_flagged": len(flagged_set),
                "n_events": len(g["events"]),
                "n_events_with_hit": events_with_hit,
                "n_hit_days": len(hit_days),
                "n_false_alarm_days": len(false_alarm_days),
                "sliding_baseline_rate": None if baseline_rate != baseline_rate else round(baseline_rate, 4),
            }

            total_flagged += len(flagged_set)
            total_hit_days += len(hit_days)
            total_false_alarm_days += len(false_alarm_days)
            total_precursor_days += len(precursor_days & eval_dates)
            total_days += len(eval_dates)
            n_events += len(g["events"])
            n_events_with_hit += events_with_hit

        total_non_precursor_days = total_days - total_precursor_days
        windows_out[str(w)] = {
            "window_days": w,
            "total_flagged_days": total_flagged,
            "total_hit_days": total_hit_days,
            "total_false_alarm_days": total_false_alarm_days,
            "precision_hit_over_flagged": round(total_hit_days / total_flagged, 4) if total_flagged else None,
            "n_events": n_events,
            "n_events_with_hit": n_events_with_hit,
            "recall_events_with_hit": round(n_events_with_hit / n_events, 4) if n_events else None,
            "false_alarm_rate_per_non_precursor_day": (
                round(total_false_alarm_days / total_non_precursor_days, 5) if total_non_precursor_days else None
            ),
            "per_group": per_group_detail,
        }

    return {"band": band, "threshold": FIXED_RULE_THRESHOLD, "windows": windows_out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--min-mag", type=float, required=True,
                    help="magnitude threshold of this tier; must match the --catalog file's")
    args = ap.parse_args()

    out_dir = PROJECT_DIR / "data" / "interim" / "backtest"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for band in BANDS:
        r = run_band(band, args.catalog, args.min_mag)
        r["min_mag"] = args.min_mag
        results[band] = r
        for w in WINDOWS_DAYS:
            wr = r["windows"][str(w)]
            print(f"[{band}] window={w}d flagged={wr['total_flagged_days']} "
                  f"hits={wr['total_hit_days']} false_alarms={wr['total_false_alarm_days']} "
                  f"precision={wr['precision_hit_over_flagged']} "
                  f"recall={wr['recall_events_with_hit']} "
                  f"({wr['n_events_with_hit']}/{wr['n_events']} events)", file=sys.stderr)

    out_path = out_dir / f"z_rule_backtest_{args.label}.json"
    out_path.write_text(json.dumps(results, indent=2))
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
