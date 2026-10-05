"""多事件地磁前兆分析流程的逐組設定
載入器。

這個模組原本放的是固定常數（GDMS_DIR/EQ_*/STATIONS/
SCALAR_ONLY_STATIONS/NEAR_STATIONS/FAR_STATIONS/KNOWN_OUTAGE_WINDOWS），
只適用於單一 G10（2024-04-03 M7.2）事件。`load_group_config`
用逐組載入器取代了這一切，由 `events.py`（事件
中繼資料）和每個檔案自己的 IAGA-2002 檔頭（測站中繼資料 + 純量
vs 向量狀態，透過 `parser.parse_header`）驅動，所以在
全部 13 組上都能一致運作，不需要人工維護的分年代測站表。
"""
from __future__ import annotations

import json
import math
import re
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from events import Event, get_group  # noqa: E402
from parser import DayFileRef, parse_header  # noqa: E402

_TGZ_MEMBER_RE = re.compile(r"([a-z]{3})(\d{8})dsec\.sec$")

PROJECT_DIR = Path(__file__).resolve().parent.parent
GX_DATA_ROOT = PROJECT_DIR  # 2026-09-14 搬家：G1..G23 現在放在 geomag_precursor/ 裡，不在上一層
OUTPUT_DIR = PROJECT_DIR / "output"

# 已由 timezone_check.py 實際確認（Sq 日變化曲線 + 工業
# 雜訊日變化曲線，在 G10/2024 資料上互相印證）原始
# TIME 欄位是 UTC。所有組別都是同一個 CWA IAGA-2002 來源／格式，所以
# 當成全資料集共用的常數，不逐組重驗。
DATA_TIMEZONE = "UTC"
LOCAL_UTC_OFFSET_HOURS = 8

# 一個通道類型測站池（只有 F 或 XYZ）至少要有這麼多測站，才會嘗試
# 近站／遠站共模迴歸篩檢方法；低於這個數字，
# 測站池標記為不足，該組就跳過這個方法。
MIN_STATIONS_FOR_METHOD = 5
N_NEAR_STATIONS = 3
N_FAR_STATIONS = 2

# 一個測站–日缺漏樣本比例超過這個值，就從
# 「乾淨」基準日中排除——這個從資料自動偵測的機制
# 取代舊的人工整理 KNOWN_OUTAGE_WINDOWS 清單，作為通用的
# 逐組中斷排除邏輯（見下面的 auto_outage_dates），因為
# 替 13 組 x 最多 19 站人工整理窗口清單，無法像
# 只有一組時那樣擴展。
OUTAGE_PCT_MISSING_THRESHOLD = 0.05

# G10 原本人工整理的中斷清單（附逐次說明，例如
# "twu ~86% missing on 2024-04-23"）。保留在這裡，作為
# 單一事件 2024 流程當時人工確認內容的歷史紀錄——
# 通用的多組邏輯已不再查詢它（下面的 auto_outage_dates
# 處理 G10 的方式和其他組完全相同），但仍保留，
# 因為 G10/NOTES.md 有引用它，而且這些逐次說明
# 提供了光看 pct_missing 得不到的脈絡。
G10_KNOWN_OUTAGE_WINDOWS = [
    {"station": "ALL", "start": "2024-01-04 23:20:00", "end": "2024-01-04 23:59:59",
     "note": "全網短暫中斷"},
    {"station": "ALL", "start": "2024-03-02 23:20:00", "end": "2024-03-02 23:59:59",
     "note": "全網短暫中斷"},
    {"station": "ALL", "start": "2024-03-05 23:20:00", "end": "2024-03-05 23:59:59",
     "note": "全網短暫中斷"},
    {"station": "ALL", "start": "2024-03-06 23:20:00", "end": "2024-03-06 23:59:59",
     "note": "全網短暫中斷"},
    {"station": "twu", "start": "2024-04-19 00:00:00", "end": "2024-04-19 23:59:59",
     "note": "缺 404 秒"},
    {"station": "twu", "start": "2024-04-23 00:00:00", "end": "2024-04-23 23:59:59",
     "note": "約 86% 缺漏"},
    {"station": "twu", "start": "2024-04-24 00:00:00", "end": "2024-04-25 23:59:59",
     "note": "完全缺漏"},
    {"station": "zbn", "start": "2024-04-13 16:10:00", "end": "2024-04-13 23:59:59",
     "note": "震後缺口"},
    {"station": "kma", "start": "2024-01-10 00:00:00", "end": "2024-01-10 23:59:59",
     "note": "約 6 小時缺口"},
    {"station": "yhg", "start": "2024-01-15 00:00:00", "end": "2024-01-15 23:59:59",
     "note": "約 11 小時缺口"},
]


def list_tgz_files(gdms_dir: Path) -> list[Path]:
    """組資料夾內直接放著的所有 GDMS 批次下載壓縮檔，
    依檔名排序。同一測站＋日期出現在多個 .tgz 時，
    「第一個 .tgz 優先」指的就是這個排序（它們的
    檔名不能可靠地反映日期範圍內容，所以這只是
    一個決定性的同分處理，不代表哪個壓縮檔
    「比較新」）。"""
    return sorted(gdms_dir.glob("*.tgz"))


def _tgz_member_index(tgz_path: Path) -> list[dict]:
    """列出 .tgz 內每一個 <station><YYYYMMDD>dsec.sec 成員，不
    取出任何檔案內容——只讀 tarfile 自己的成員標頭。
    gzip 不支援隨機存取，所以完整掃過一個大壓縮檔
    （實測一個 276MB／481 檔的樣本約 30 秒）並不便宜；結果會
    快取到壓縮檔旁邊的 `<name>.tgz.idx.json` 附屬檔，以
    (size, mtime) 為鍵，所以用同名但內容不同的檔案替換 .tgz
    會自動讓快取失效。"""
    idx_path = tgz_path.with_name(tgz_path.name + ".idx.json")
    st = tgz_path.stat()
    if idx_path.exists():
        try:
            cached = json.loads(idx_path.read_text())
        except (json.JSONDecodeError, OSError):
            cached = {}
        if cached.get("size") == st.st_size and cached.get("mtime") == st.st_mtime:
            return cached["members"]

    members = []
    with tarfile.open(tgz_path, "r:gz") as tf:
        for m in tf:  # 只走標頭；不解壓成員內容
            if not m.isfile():
                continue
            match = _TGZ_MEMBER_RE.search(m.name)
            if match:
                members.append({"station": match.group(1), "date_str": match.group(2), "member": m.name})
    idx_path.write_text(json.dumps({"size": st.st_size, "mtime": st.st_mtime, "members": members}))
    return members


def list_day_refs(gdms_dir: Path, station: str = "*") -> list[DayFileRef]:
    """列出一個測站（或用預設的 '*' 列出所有測站）的 DayFileRef，
    涵蓋零散的 .sec/.sec.gz 檔和 gdms_dir 內直接放著的任何 *.tgz 批次
    壓縮檔——這是唯一定義優先順序的地方，所以
    每支列出日檔的腳本都會以同樣方式支援 .tgz。

    同一測站＋日期出現在多個來源時的優先順序：
    零散的 .sec/.sec.gz 永遠優先於 .tgz 內容（不出聲——這是
    預期的常態，不值得警告）；多個
    .tgz 壓縮檔之間，依 list_tgz_files() 順序先遇到的優先，
    並印出警告，讓重疊的批次下載被注意到。
    """
    by_key: dict[tuple[str, str], DayFileRef] = {}

    pattern = "*dsec.sec*" if station == "*" else f"{station}*dsec.sec*"
    for p in sorted(gdms_dir.glob(pattern)):
        if p.suffix not in (".sec", ".gz"):
            continue
        stem = p.name[:-3] if p.name.endswith(".sec.gz") else p.name
        st, date_str = stem[:3], stem[3:11]
        key = (st, date_str)
        cur = by_key.get(key)
        # 如果 .sec 和 .sec.gz 都存在，優先用一般 .sec。
        if cur is None or (cur.source_path.suffix == ".gz" and p.suffix == ".sec"):
            by_key[key] = DayFileRef(st, date_str, p, None)

    n_loose = len(by_key)
    n_tgz = 0
    n_dup = 0
    for tgz in list_tgz_files(gdms_dir):
        for rec in _tgz_member_index(tgz):
            st, date_str, member = rec["station"], rec["date_str"], rec["member"]
            if station != "*" and st != station:
                continue
            key = (st, date_str)
            existing = by_key.get(key)
            if existing is None:
                by_key[key] = DayFileRef(st, date_str, tgz, member)
                n_tgz += 1
            elif existing.member is not None:
                n_dup += 1
                print(
                    f"[list_day_refs] WARNING: {st}{date_str} duplicate -- "
                    f"{tgz.name}:{member} 已略過，已經在使用 "
                    f"{existing.source_path.name}:{existing.member}",
                    file=sys.stderr,
                )
            # 否則：既有的是零散檔，永遠不出聲地優先。

    if n_tgz or n_dup:
        print(
            f"[list_day_refs] {gdms_dir.name}: {n_loose} loose, {n_tgz} from .tgz, "
            f"{n_dup} 個重複的 tgz 條目已略過",
            file=sys.stderr,
        )
    return sorted(by_key.values(), key=lambda r: (r.station, r.date_str))


def resolve_day_ref(gdms_dir: Path, station: str, date_str: str) -> DayFileRef | None:
    """回傳一個測站–日的 DayFileRef，有一般 .sec 就優先用，
    否則用 .sec.gz，再不然就搜尋 gdms_dir 內每一個 .tgz；哪裡都
    找不到則回傳 None。"""
    plain = gdms_dir / f"{station}{date_str}dsec.sec"
    if plain.exists():
        return DayFileRef(station, date_str, plain, None)
    gz = gdms_dir / f"{station}{date_str}dsec.sec.gz"
    if gz.exists():
        return DayFileRef(station, date_str, gz, None)
    for tgz in list_tgz_files(gdms_dir):
        for rec in _tgz_member_index(tgz):
            if rec["station"] == station and rec["date_str"] == date_str:
                return DayFileRef(station, date_str, tgz, rec["member"])
    return None


def to_local_hour(utc_index):
    """給定 UTC DatetimeIndex，回傳台灣當地（UTC+8）的時（0-23）。"""
    import pandas as pd

    return (utc_index + pd.Timedelta(hours=LOCAL_UTC_OFFSET_HOURS)).hour


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(min(1.0, a ** 0.5))


def auto_outage_dates(daily_features_df, threshold: float = OUTAGE_PCT_MISSING_THRESHOLD) -> dict[str, set[str]]:
    """每個測站從「乾淨」基準日中排除的日期集合，由
    daily_features.csv 自己的 pct_missing 欄位推出，而不是
    人工整理的清單——這是舊的
    只限 G10 的 KNOWN_OUTAGE_WINDOWS 機制的通用逐組替代品。daily_features_df 必須有
    'station'、'date'、'pct_missing' 欄位（和
    build_daily_features.py 寫出的一樣）。"""
    dirty = daily_features_df[daily_features_df["pct_missing"] > threshold]
    out: dict[str, set[str]] = {}
    for station, g in dirty.groupby("station"):
        out[station] = set(g["date"])
    return out


@dataclass(frozen=True)
class StationPool:
    """一組內單一通道類型測站池（只有 F 或 XYZ）的近站／遠站劃分。"""

    channel: str  # "F" 或 "XYZ"
    all_stations: tuple[str, ...]  # 這個測站池在該組中出現的每個測站代碼，由近到遠
    near: tuple[str, ...]
    far: tuple[str, ...]
    sufficient: bool  # 如果這個測站池的測站數足以嘗試近站／遠站方法則為 True


@dataclass(frozen=True)
class GroupConfig:
    group_id: str
    gdms_dir: Path
    interim_dir: Path
    external_dir: Path
    anchor_event: Event
    all_events: tuple[Event, ...]
    stations: dict  # 測站代碼 -> {"name","lat","lon","elevation_m","reported","distance_km"}
    f_pool: StationPool
    xyz_pool: StationPool


def _discover_stations(gdms_dir: Path, anchor_event: Event) -> dict:
    """用 glob 找出這組資料夾中出現的不同測站代碼前綴，
    每個代碼讀一個代表性檔案的檔頭，取得名稱／緯度／經度／
    高程／Reported 狀態。純量 vs 向量狀態在
    單一組內是固定的——已在全部 13 組中實際確認，F->XYZF
    升級總是落在組別年代的邊界上，從不在組內中途發生——所以每個
    代碼讀一個代表性檔案就夠了。"""
    codes = sorted({r.station for r in list_day_refs(gdms_dir)})
    stations = {}
    for code in codes:
        sample = list_day_refs(gdms_dir, code)[0]
        h = parse_header(sample)
        stations[code] = {
            "name": h["station_name"],
            "lat": h["lat"],
            "lon": h["lon"],
            "elevation_m": h["elevation_m"],
            "reported": h["reported"].strip().upper(),
            "distance_km": round(haversine_km(h["lat"], h["lon"], anchor_event.lat, anchor_event.lon), 1),
        }
    return stations


def _build_pool(stations: dict, channel: str) -> StationPool:
    wanted_reported = "F" if channel == "F" else "XYZF"
    ranked = sorted(
        (c for c, m in stations.items() if m["reported"] == wanted_reported),
        key=lambda c: stations[c]["distance_km"],
    )
    near = tuple(ranked[:N_NEAR_STATIONS])
    far = tuple(c for c in ranked[-N_FAR_STATIONS:] if c not in near)
    return StationPool(
        channel=channel,
        all_stations=tuple(ranked),
        near=near,
        far=far,
        sufficient=len(ranked) >= MIN_STATIONS_FOR_METHOD and len(near) >= 1 and len(far) >= 1,
    )


def load_group_config(group_id: str) -> GroupConfig:
    group = get_group(group_id)
    gdms_dir = GX_DATA_ROOT / group.folder
    interim_dir = PROJECT_DIR / "data" / "interim" / group_id
    external_dir = PROJECT_DIR / "data" / "external" / group_id
    interim_dir.mkdir(parents=True, exist_ok=True)
    external_dir.mkdir(parents=True, exist_ok=True)

    anchor = group.anchor_event
    stations = _discover_stations(gdms_dir, anchor)
    f_pool = _build_pool(stations, "F")
    xyz_pool = _build_pool(stations, "XYZ")

    return GroupConfig(
        group_id=group_id,
        gdms_dir=gdms_dir,
        interim_dir=interim_dir,
        external_dir=external_dir,
        anchor_event=anchor,
        all_events=group.events,
        stations=stations,
        f_pool=f_pool,
        xyz_pool=xyz_pool,
    )


if __name__ == "__main__":
    group_ids = sys.argv[1:] or ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13"]
    for gid in group_ids:
        cfg = load_group_config(gid)
        a = cfg.anchor_event
        print(f"== {gid} == anchor {a.date} {a.magnitude_type}{a.magnitude} ({a.lat},{a.lon}) "
              f"confidence={a.coord_confidence}")
        print(f"  F   測站池：{len(cfg.f_pool.all_stations):2d} 站  near={cfg.f_pool.near}  "
              f"far={cfg.f_pool.far}  sufficient={cfg.f_pool.sufficient}")
        print(f"  XYZ 測站池：{len(cfg.xyz_pool.all_stations):2d} 站  near={cfg.xyz_pool.near}  "
              f"far={cfg.xyz_pool.far}  sufficient={cfg.xyz_pool.sufficient}")
