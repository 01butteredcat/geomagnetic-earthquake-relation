"""抓取一組實際日期範圍內的歷史 Kp（GFZ Potsdam）和 Dst（Kyoto WDC）太空天氣
指數，並推出一份磁暴日旗標清單，
用來避免把磁暴造成的變動誤讀成
構造前兆訊號。

用法：fetch_space_weather.py --group G10

資料來源（純文字、不需驗證，每組抓一次並快取）：
  Kp:  https://kp.gfz.de/kpdata?startdate=...&enddate=...&format=kp2
       （kp.gfz-potsdam.de 會 301 轉址到這裡；-L 會跟著轉）。Kp 從
       1932 年就有紀錄，所以任何一組的日期範圍（2017-2026）都有涵蓋。
  Dst: https://wdc.kugi.kyoto-u.ac.jp/dst_<final|provisional|realtime>/<YYYYMM>/dst<YYMM>.for.request
       （需要像瀏覽器的 User-Agent，否則 Kyoto 會回 403）。final/provisional/realtime
       哪一種可用，取決於那個月份距今
       多久——每個月依這個順序嘗試，因為較舊的組別比較可能有定案資料，
       最新的那組（G13，2026）比較可能
       只有 realtime。

如果兩者都抓不到，就退回內部代理指標：該組遠端參考站
同時出現大偏差的日子
標記為「全球擾動」，標示為低信心。如果只缺一個來源
（或只缺部分 Dst 月份），仍會使用已抓到部分的磁暴日，
摘要的 "confidence" 會寫明缺了哪一部分
（"medium (...)"），而不是宣稱 Kp+Dst 完整涵蓋。
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import list_day_refs, load_group_config  # noqa: E402

KP_URL_TMPL = "https://kp.gfz.de/kpdata?startdate={start}&enddate={end}&format=kp2"
DST_URL_TMPL = "https://wdc.kugi.kyoto-u.ac.jp/dst_{kind}/{yyyymm}/dst{yymm}.for.request"
DST_KINDS = ("final", "provisional", "realtime")
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"

KP_STORM_THRESHOLD = 5.0  # G1+ 地磁暴
DST_STORM_THRESHOLD = -30.0  # nT 谷值


def _curl(url: str, extra_args: list[str] | None = None) -> str | None:
    args = ["curl", "-sL", "--fail", "-A", UA]
    if extra_args:
        args += extra_args
    args.append(url)
    try:
        out = subprocess.run(args, capture_output=True, timeout=30, check=True)
        return out.stdout.decode("utf-8", errors="replace")
    except Exception as exc:  # noqa: BLE001
        print(f"fetch failed for {url}: {exc}", file=sys.stderr)
        return None


def _group_date_range(cfg) -> tuple[str, str]:
    dates = sorted(r.date_str for r in list_day_refs(cfg.gdms_dir))
    d0, d1 = dates[0], dates[-1]
    return f"{d0[:4]}-{d0[4:6]}-{d0[6:8]}", f"{d1[:4]}-{d1[4:6]}-{d1[6:8]}"


def _month_range(start: str, end: str) -> list[str]:
    months = pd.period_range(start=start, end=end, freq="M")
    return [str(m).replace("-", "") for m in months]


def fetch_kp(start: str, end: str) -> pd.DataFrame | None:
    text = _curl(KP_URL_TMPL.format(start=start, end=end))
    if not text:
        return None
    rows = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) < 8:
            continue
        year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
        start_hr = float(parts[3])
        kp = float(parts[7])
        rows.append({"date": f"{year:04d}{month:02d}{day:02d}", "start_hour_utc": start_hr, "kp": kp})
    return pd.DataFrame(rows)


def fetch_dst(start: str, end: str) -> tuple[pd.DataFrame | None, list[str]]:
    """回傳（抓到的月份組成的 dataframe，一個都沒抓到則為 None，
    抓不到的 YYYYMM 月份清單）。每個月依序嘗試
    final/provisional/realtime，整個流程會再重試一次：以前
    一次暫時性失敗就會丟掉整個範圍（而某組
    快取的 "FETCH FAILED" 會因此卡好幾週——G1/G2/G3/G11）。"""
    rows = []
    missing: list[str] = []
    for yyyymm in _month_range(start, end):
        yymm = yyyymm[2:]
        text = None
        for _attempt in range(2):
            for kind in DST_KINDS:
                url = DST_URL_TMPL.format(kind=kind, yyyymm=yyyymm, yymm=yymm)
                text = _curl(url)
                if text:
                    break
            if text:
                break
        if not text:
            missing.append(yyyymm)
            continue
        for line in text.strip().splitlines():
            if not line.startswith("DST"):
                continue
            # 例如 "DST2404*01PPX120   0  -4  -5  -9 ... -7"
            head, rest = line[:20], line[20:]
            day = int(head[8:10])
            vals = rest.split()
            if len(vals) < 25:
                continue
            hourly = [float(v) for v in vals[:24]]
            daily_mean = float(vals[24])
            date_str = f"{yyyymm}{day:02d}"
            for hr, v in enumerate(hourly, start=1):
                rows.append({"date": date_str, "hour_utc": hr % 24, "dst": v})
            rows.append({"date": date_str, "hour_utc": -1, "dst": daily_mean})  # -1 = 日平均標記
    return (pd.DataFrame(rows) if rows else None), missing


def internal_proxy_storm_days(cfg) -> pd.DataFrame:
    """備援：把遠端參考站同時出現
    大 H/F 偏差的日子標記出來，作為內部推斷（信心較低）的磁暴
    代理指標。有 XYZ 遠站池就用該組的 XYZ 遠站池，否則用 F 遠站池。"""
    far = cfg.xyz_pool.far or cfg.f_pool.far
    range_col = "H_range" if cfg.xyz_pool.far else "F_range"
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[daily.station.isin(far)].copy()
    sub["z"] = sub.groupby("station")[range_col].transform(
        lambda s: (s - s.median()) / (1.4826 * (s - s.median()).abs().median() + 1e-9)
    )
    agg = sub.groupby("date")["z"].median().reset_index()
    agg["storm_flag"] = agg["z"] > 3.0
    agg["source"] = f"internal_proxy (far-station {range_col} MAD z>3, LOW CONFIDENCE)"
    return agg[agg.storm_flag][["date", "source"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    ap.add_argument("--check-cache", action="store_true",
                    help="exit 0 if storm_days.csv is high-confidence AND covers the folder's whole date range, else 1")
    args = ap.parse_args()
    cfg = load_group_config(args.group)
    start, end = _group_date_range(cfg)
    if args.check_cache:
        # 資料夾資料延長之前建的快取（G6/G7/G8、G17，直到 2026-09-25）仍
        # 寫 "high"，但後來的事件就沒有磁暴旗標，所以也要檢查涵蓋範圍。
        summary_path = cfg.interim_dir / "storm_days_summary.json"
        ok = (cfg.interim_dir / "storm_days.csv").exists() and summary_path.exists()
        if ok:
            summary = json.loads(summary_path.read_text())
            ok = summary.get("confidence", "").startswith("high") and summary["date_range"] == [start, end]
        sys.exit(0 if ok else 1)
    print(f"[{args.group}] fetching space weather for {start}..{end}", file=sys.stderr)

    kp_df = fetch_kp(start, end)
    dst_df, dst_missing_months = fetch_dst(start, end)

    source_note = []
    storm_days = set()

    if kp_df is not None and len(kp_df):
        kp_df.to_csv(cfg.external_dir / "kp.csv", index=False)
        kp_daily_max = kp_df.groupby("date")["kp"].max()
        kp_storm_days = set(kp_daily_max[kp_daily_max >= KP_STORM_THRESHOLD].index)
        storm_days |= kp_storm_days
        source_note.append(f"Kp (GFZ): {len(kp_df)} 3h records, {len(kp_storm_days)} storm days (max Kp>={KP_STORM_THRESHOLD})")
    else:
        source_note.append("Kp (GFZ): FETCH FAILED")

    if dst_df is not None and len(dst_df):
        dst_df.to_csv(cfg.external_dir / "dst.csv", index=False)
        dst_hourly = dst_df[dst_df.hour_utc >= 0]
        dst_daily_min = dst_hourly.groupby("date")["dst"].min()
        dst_storm_days = set(dst_daily_min[dst_daily_min <= DST_STORM_THRESHOLD].index)
        storm_days |= dst_storm_days
        source_note.append(f"Dst (Kyoto WDC): {len(dst_df)} records, {len(dst_storm_days)} storm days (min Dst<={DST_STORM_THRESHOLD}nT)")
        if dst_missing_months:
            source_note.append(f"Dst (Kyoto WDC): PARTIAL, missing months {dst_missing_months}")
    else:
        source_note.append("Dst (Kyoto WDC): FETCH FAILED")

    kp_ok = kp_df is not None and len(kp_df) > 0
    dst_ok = dst_df is not None and len(dst_df) > 0
    if kp_ok and dst_ok and not dst_missing_months:
        confidence = "high (official Kp/Dst indices)"
    elif kp_ok or dst_ok:
        gaps = []
        if not kp_ok:
            gaps.append("Kp unavailable")
        if not dst_ok:
            gaps.append("Dst unavailable")
        elif dst_missing_months:
            gaps.append(f"Dst missing {len(dst_missing_months)} month(s)")
        confidence = f"medium (official indices incomplete: {', '.join(gaps)})"
    else:
        proxy = internal_proxy_storm_days(cfg)
        storm_days |= set(proxy["date"])
        source_note.append(f"FALLBACK: internal far-station proxy, {len(proxy)} storm days")
        confidence = "low (internal proxy only, official sources unavailable)"

    # 每個磁暴日之後加 1-2 天恢復期（環電流衰減）
    all_dates = sorted(storm_days)
    extended = set(storm_days)
    for d in all_dates:
        dt = pd.to_datetime(d, format="%Y%m%d")
        for delta in (1, 2):
            extended.add((dt + pd.Timedelta(days=delta)).strftime("%Y%m%d"))

    storm_df = pd.DataFrame(sorted(extended), columns=["date"])
    storm_df["is_storm_or_recovery"] = True
    storm_df["is_storm_onset"] = storm_df["date"].isin(storm_days)
    storm_df.to_csv(cfg.interim_dir / "storm_days.csv", index=False)

    summary = {
        "group": args.group,
        "date_range": [start, end],
        "sources": source_note,
        "confidence": confidence,
        "dst_missing_months": dst_missing_months,
        "n_storm_onset_days": len(storm_days),
        "n_storm_or_recovery_days": len(extended),
        "storm_onset_dates": sorted(storm_days),
    }
    (cfg.interim_dir / "storm_days_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
