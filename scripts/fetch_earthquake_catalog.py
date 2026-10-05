"""抓取地磁測站網附近的**擴充**地震目錄，
涵蓋每一起 M>=6.0 事件**以及**其間較小的事件——這是教授建議方法中
「把一次地震變成很多次」的那一塊：
`events.py` 只有用來定義 13 組的 20 起人工挑選 M>=6.0 事件，
不足以做真正的疊加時間
疊加，也不足以構成有意義的回測母體。

來源：使用者提供的 CWA GDMS 目錄匯出檔（ML，UTC）——本專案
以 CWA 作為所有事件和規模的權威來源。兩份匯出檔接續涵蓋
2009-01-01 ~ 2026-07-31（CWA_CATALOGS）；每組的窗口會同時從兩份讀取，所以
橫跨兩個檔案的窗口（G19）仍然是純 CWA。USGS FDSN 查詢
只保留在 --allow-usgs 之後，給這個範圍以外的窗口用：在
2026-09-30 之前，它是 2024-09~2026-07 以外每一組都會默默退回的來源，
因而把 Mw/mb 規模混進了 M>=5.0/5.5 級距。

範圍：
  - 只有 8 個有 `ulf_near_far_index.csv` 的組別（向量站足夠的
    XYZ 測站池；G1/G2/G3 只有純量，根本跑不了極化方法，
    所以不屬於疊加時間／回測母體）。
  - 每組的日期窗口 = 該組 ULF 指標
    實際涵蓋的窗口（直接從 CSV 讀取，不從
    events.py／文件重新推導，所以永遠不會和實際可用的資料脫節）。
  - 邊界框涵蓋整個測站網並留有餘裕：緯度
    20.5-27（測站網從 hcn 21.94N 到 mtu 26.17N），經度 117.5-123.5
    （kma 118.35E 到 122E 以外的外海事件）。
  - 去叢集：依時間排序；一個事件如果落在某個已保留、規模
    >= 它自己的事件的 `decluster_days` 天**且** `decluster_km` km 之內，
    就丟掉（也就是每個緊密的
    時空叢集只保留最大的事件）——避免把一個餘震序列當成許多
    獨立樣本，和 events.py 自己的
    docstring 針對 13 組設計提出的偽重複疑慮相同。
  - 和 events.py 已登錄的事件交叉比對（發震時間
    在 KNOWN_EVENT_SEC 內、震央在 KNOWN_EVENT_KM 內；**不**比較
    規模，因為同一起地震的 USGS 和 CWA ML 常常相差超過
    0.3），標記為 `is_known_event`，而不是
    重複計入。

用法：
  fetch_earthquake_catalog.py --min-mag 5.5
  fetch_earthquake_catalog.py --min-mag 5.0 --output data/external/extended_catalog_m5.0.csv
  fetch_earthquake_catalog.py --min-mag 5.5 --reflag   # 只離線重算 is_known_event
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR, haversine_km, load_group_config  # noqa: E402
from events import ALL_GROUP_IDS, assign_group_for_time, folder_events, get_group, sibling_group_ids  # noqa: E402

USGS_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BBOX = {"minlatitude": 20.5, "maxlatitude": 27.0, "minlongitude": 117.5, "maxlongitude": 123.5}
ULF_GROUPS = ("G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23", "G24")

DECLUSTER_DAYS = 3
DECLUSTER_KM = 100
# 目錄和 events.py 中的同一起地震：發震時間相差一分鐘以內，震央
# 相距 50 km 以內。舊規則（6 小時、規模相差 0.3 以內）會漏掉 USGS 規模
# 和 CWA ML 相差超過 0.3 的重複事件（G4 2020-12-10：6.1 vs 6.64，相差 0 秒），也會配到
# 距離已登錄事件數小時的不同餘震。
KNOWN_EVENT_SEC = 60
KNOWN_EVENT_KM = 50

# 使用者提供的 CWA GDMS 區域規模報告匯出檔（以空白分隔，標題列
# "date time lat lon depth ML nstn dmin gap trms ERH ERZ fixed nph quality"），M>=5.0，
# 格式為（路徑, 第一天, 最後一天），對應每份匯出檔請求的範圍。兩份都是 UTC：
# GDMScatalog.txt 由 G11 的 2025-01-21 錨點確認（它的 2025-01-20 16:17 UTC 那一列只有在
# 台灣當地時間才落在 01-21）；GDMScatalog_2009-2024.txt 由 events.py 在 2024-09 之前的
# 全部 62 起事件都在 2 秒內對到一列確認。匯出檔涵蓋的範圍比 BBOX 大（經度延伸到
# 125.6），所以 fetch_cwa() 會自己套用 BBOX。
CWA_CATALOGS = (
    (PROJECT_DIR.parent / "GDMScatalog_2009-2024.txt", "2009-01-01", "2024-08-31"),
    (PROJECT_DIR.parent / "GDMScatalog.txt", "2024-09-01", "2026-07-31"),
)
CWA_CATALOG_START = CWA_CATALOGS[0][1]
CWA_CATALOG_END = CWA_CATALOGS[-1][2]


def group_date_window(group_id: str) -> tuple[str, str] | None:
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_near_far_index.csv"
    if not path.exists():
        return None
    dates = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            dates.append(row["date"])
    if not dates:
        return None
    dates.sort()
    d0, d1 = dates[0], dates[-1]
    return (f"{d0[:4]}-{d0[4:6]}-{d0[6:8]}", f"{d1[:4]}-{d1[4:6]}-{d1[6:8]}")


def fetch_usgs(start: str, end: str, min_mag: float) -> list[dict]:
    params = {**BBOX, "format": "geojson", "starttime": start, "endtime": end, "minmagnitude": min_mag}
    url = USGS_URL + "?" + "&".join(f"{k}={v}" for k, v in params.items())
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 2:
                raise
            print(f"  發生錯誤後重試：{exc}", file=sys.stderr)
            time.sleep(2)
    events = []
    for feat in data["features"]:
        p = feat["properties"]
        lon, lat, depth = feat["geometry"]["coordinates"]
        t = datetime.fromtimestamp(p["time"] / 1000, tz=timezone.utc)
        events.append({
            "time_utc": t.strftime("%Y-%m-%d %H:%M:%S"),
            "_dt": t,
            "lat": lat, "lon": lon, "depth_km": depth,
            "mag": p["mag"], "place": p.get("place") or "",
            "usgs_id": feat.get("id", ""),
        })
    return events


def in_bbox(lat: float, lon: float) -> bool:
    return (BBOX["minlatitude"] <= lat <= BBOX["maxlatitude"]
            and BBOX["minlongitude"] <= lon <= BBOX["maxlongitude"])


def fetch_cwa_all(paths: list[Path], start: str, end: str, min_mag: float) -> list[dict]:
    """對每一份匯出檔執行 fetch_cwa()，讓橫跨兩份的窗口被完整讀取。"""
    return [e for path in paths for e in fetch_cwa(path, start, end, min_mag)]


def fetch_cwa(path: Path, start: str, end: str, min_mag: float) -> list[dict]:
    """解析一份使用者提供的 CWA GDMS 目錄匯出檔，篩選到 [start, end]
    （含端點，YYYY-MM-DD）、mag >= min_mag 和 BBOX（和 fetch_usgs() 查詢的框相同）。
    回傳和 fetch_usgs() 相同結構的 dict，讓 decluster()/flag_known_events() 不用改
    就能用。"""
    start_dt = datetime.strptime(start, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    end_dt = datetime.strptime(end, "%Y-%m-%d").replace(tzinfo=timezone.utc) + timedelta(days=1)
    events = []
    with path.open() as f:
        header = f.readline()
        assert header.split()[:6] == ["date", "time", "lat", "lon", "depth", "ML"], \
            f"非預期的 GDMS 目錄標題列：{header!r}"
        for line in f:
            parts = line.split()
            if not parts:
                continue
            date_s, time_s, lat_s, lon_s, depth_s, mag_s = parts[:6]
            mag = float(mag_s)
            if mag < min_mag:
                continue
            t = datetime.strptime(f"{date_s} {time_s}", "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
            if not (start_dt <= t < end_dt):
                continue
            if not in_bbox(float(lat_s), float(lon_s)):
                continue
            events.append({
                "time_utc": t.strftime("%Y-%m-%d %H:%M:%S"),
                "_dt": t,
                "lat": float(lat_s), "lon": float(lon_s), "depth_km": float(depth_s),
                "mag": mag, "place": "", "usgs_id": "",
            })
    return events


def decluster(events: list[dict], days: float = DECLUSTER_DAYS, km: float = DECLUSTER_KM) -> list[dict]:
    kept: list[dict] = []
    for e in sorted(events, key=lambda x: x["_dt"]):
        suppressed = False
        for k in kept:
            dt_days = abs((e["_dt"] - k["_dt"]).total_seconds()) / 86400
            if dt_days > days:
                continue
            if haversine_km(e["lat"], e["lon"], k["lat"], k["lon"]) > km:
                continue
            if k["mag"] >= e["mag"]:
                suppressed = True
                break
        if not suppressed:
            kept.append(e)
    kept_set = {id(k) for k in kept}
    for e in events:
        e["declustered"] = id(e) in kept_set
    return events


def flag_known_events(events: list[dict], group_id: str) -> None:
    # folder_events：兄弟組的已登錄事件和這一組自己的一樣算「已知」。
    known = folder_events(group_id)
    for e in events:
        e["is_known_event"] = False
        for ev in known:
            known_dt = datetime.strptime(ev.time_utc, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            if (abs((e["_dt"] - known_dt).total_seconds()) <= KNOWN_EVENT_SEC
                    and haversine_km(e["lat"], e["lon"], ev.lat, ev.lon) <= KNOWN_EVENT_KM):
                e["is_known_event"] = True
                break


def reflag_existing(path: Path) -> None:
    """只在已經抓好的目錄 CSV 上離線重算 is_known_event：資料列、
    去叢集和組別歸屬都維持抓取時的樣子。"""
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    n_changed = 0
    for group_id in dict.fromkeys(r["group"] for r in rows):
        grp = [r for r in rows if r["group"] == group_id]
        evs = [{"_dt": datetime.strptime(r["time_utc"][:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc),
                "lat": float(r["lat"]), "lon": float(r["lon"])} for r in grp]
        flag_known_events(evs, group_id)
        for r, e in zip(grp, evs):
            new = str(e["is_known_event"])
            if r["is_known_event"] != new:
                n_changed += 1
                print(f"  {group_id} {r['time_utc']} M{r['mag']} {r['source']}: "
                      f"is_known_event {r['is_known_event']} -> {new}", file=sys.stderr)
                r["is_known_event"] = new
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"已重新標記 {path}：{n_changed} 列有變動", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-mag", type=float, default=5.5)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--groups", nargs="*", default=list(ULF_GROUPS))
    ap.add_argument("--allow-usgs", action="store_true",
                    help="窗口落在 CWA 匯出檔範圍"
                         f"（{CWA_CATALOG_START} ~ {CWA_CATALOG_END}）以外的組改查 USGS；不加的話，"
                         "這種組會報錯，絕不默默退回")
    ap.add_argument("--reflag", action="store_true",
                    help="只在既有的輸出 CSV 上重算 is_known_event（不抓取）")
    args = ap.parse_args()

    out_path = args.output or (PROJECT_DIR / "data" / "external" / f"extended_catalog_m{args.min_mag}.csv")
    if args.reflag:
        reflag_existing(out_path)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for group_id in args.groups:
        window = group_date_window(group_id)
        if window is None:
            print(f"[{group_id}] 沒有 ulf_near_far_index.csv——跳過（只有純量或還沒跑）", file=sys.stderr)
            continue
        start, end = window
        # 實務上 USGS 的 endtime 在日界線差不多是不含端點的；多補 1 天
        end_padded = (datetime.strptime(end, "%Y-%m-%d")).strftime("%Y-%m-%d")
        if CWA_CATALOG_START <= start and end <= CWA_CATALOG_END:
            paths = [p for p, _, _ in CWA_CATALOGS]
            missing = [str(p) for p in paths if not p.exists()]
            if missing:
                sys.exit(f"缺少 CWA 目錄匯出檔：{missing}")
            print(f"[{group_id}] 讀取 CWA GDMS 目錄 {start} ~ {end_padded}，M>={args.min_mag}", file=sys.stderr)
            events = fetch_cwa_all(paths, start, end_padded, args.min_mag)
            source = "CWA_GDMS"
        elif not args.allow_usgs:
            sys.exit(f"[{group_id}] 窗口 {start} ~ {end} 在 CWA 匯出檔範圍"
                     f"（{CWA_CATALOG_START} ~ {CWA_CATALOG_END}）之外；請加入 CWA 匯出檔，或傳入 --allow-usgs")
        else:
            print(f"[{group_id}] 查詢 USGS {start} ~ {end_padded}，M>={args.min_mag}", file=sys.stderr)
            events = fetch_usgs(start, end_padded, args.min_mag)
            source = "USGS"
        for e in events:
            e["source"] = source
        events = decluster(events)
        flag_known_events(events, group_id)
        if len(sibling_group_ids(group_id)) > 1:
            # 兄弟組（共用原始資料夾）都查詢同樣的日期範圍：每個
            # 目錄事件只保留在其中一組（上面的去叢集／已知標記是先對
            # 整個窗口跑的，所以橫跨兩組的前震／餘震對仍會
            # 一起去叢集），否則它會在每個兄弟組各被算一次。
            n_all = len(events)
            events = [e for e in events if assign_group_for_time(group_id, e["_dt"]) == group_id]
            print(f"  和 {sibling_group_ids(group_id)} 共用資料夾：保留 {len(events)}/{n_all} 起"
                  f"最靠近 {group_id} 錨點的事件", file=sys.stderr)
        n_kept = sum(1 for e in events if e["declustered"])
        print(f"  {len(events)} raw -> {n_kept} declustered", file=sys.stderr)
        for e in events:
            e["group"] = group_id
            all_rows.append(e)

    fields = ["group", "time_utc", "lat", "lon", "depth_km", "mag", "place", "usgs_id",
              "declustered", "is_known_event", "source"]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r[k] for k in fields})

    n_total = len(all_rows)
    n_decl = sum(1 for r in all_rows if r["declustered"])
    print(f"已寫入 {out_path}（{n_total} 列原始資料，{n_decl} 起去叢集後的獨立事件）", file=sys.stderr)


if __name__ == "__main__":
    main()
