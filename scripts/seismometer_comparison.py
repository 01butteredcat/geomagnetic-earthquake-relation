"""同址／附近地震儀 vs 地磁的比較，針對
`coseismic_step_analysis.py` / `coseismic_stacking_analysis.py`
在每起地震發震秒附近找到的同震異常。

## 為什麼有這支腳本

那兩支腳本都分不出「磁力儀外殼被
實際震動」和「磁場本身改變了」——光看 1Hz 磁場資料，
兩者產生相同的特徵（發震秒附近的階躍／突波）。
這支腳本加上唯一**能**區分兩者的東西：同址或附近
真實地震儀／加速度儀的獨立地動資料。
邏輯（依使用者自己的框架，這支
腳本直接實作）：

  - 如果地磁異常的起始和持續時間和
    地震儀自己的強震窗口對齊，那和
    震動造成的儀器雜訊一致（雖然不能證明）。
  - 如果地磁異常在地面開始動**之前**就明顯開始，
    或在地動平息之後仍**持續**，那就傾向
    真實的地球物理機制（壓磁效應等）。

## 使用的資料（由使用者另外下載，不經由這個流程）

`events.py` 117 起事件中有 31 起有 SAC PoleZero 儀器響應
檔（`<GROUP_MMDD>/SAC_PZs_TW_<STA>_<CHAN>_...`）和對應的 miniSEED
波形檔（`<GROUP_MMDD>/<GROUP_MMDD>_w.mseed`，約 event_utc-60s 到
event_utc+600s，100Hz），都放在 `seismometer/<GROUP_MMDD>/` 底下——見
下面的 `SEISMIC_DATA_DIRS`，以 (group_id, event.date) 為鍵。（2026-08-17
重新整理：之前 mseed 檔散放在這個
專案的上層目錄，比對應的 PZ 資料夾高一層；
現在兩者一起放在 `seismometer/` 底下，每個事件一個資料夾，不和
Gx 地磁資料夾混在一起。）未涵蓋的 19 起事件中，6 起是
永久性的結構缺口：G14 的 5 起（最早的事件 2009-07-14 早於
地震資料來源 2012 年的起始）加上 G21 的 1 起（2010-11-21，同樣原因
——2026-09-19 確認資料來源可下載範圍從
2012-01-01 才開始）。其餘 13 起（分散在 G5、G8、G9、G10、G11、G12、
G15、G17）只是還沒下載。（G22 的 1 起、G23 的 1 起 + G24 的
1 起——G23/G24 在 2026-09-22 從原本有 2 起事件的單一合併 "G23" 拆出——
mseed 在 2026-09-19 下載、PZ 檔在 2026-09-20 加入，現在
已經接上；
PZ 檔組缺了幾個 mseed 裡有的測站——G22：CHK/ELD loc 11
和 HEN，G23：HEN 和 SSH——只影響那些測站的波形，它們都
不是被選中的最近測站。）
（截至 2026-08-13/14，只下載了 16 個錨點事件；
其餘 11 起非錨點事件加上先前漏掉的 G9 2022-09-17
前震在 2026-08-16 下載並確認——見 `coverage_summary.json`，
它列出 `events.py` 目前登錄表中的所有事件，附上真實的
`seismic_data_status`，讓報告永遠不會暗示比實際更多的涵蓋。）
（2026-09-23：登錄表從 49 起成長到 117 起，68 起新事件全是非錨點，
也都沒有下載地震儀資料，所以上面「19 起未涵蓋」的拆解已經過時
——它只描述最初的 19 起，不是現在約 86 起未涵蓋的；新一批的來源
見 `docs/candidate_events_gdms_2024_2026.md`。
新事件中有一起，G11 的 2025-01-21b，結果已經有真實
資料——它落在 2025-01-21 錨點已下載的 mseed 窗口內
——已接進下面的 `SEISMIC_DATA_DIRS`；見那個鍵的註解。）

（2026-09-25：其餘 80 起可下載的事件一批下載完成——
14 起 M>=6、66 起 M5——依同一個命名規則配置，
由 `_register_convention_dirs()` 推出而不是逐一列出。涵蓋率現在是 117 起中的
111 起；另外 6 起早於資料來源（G14 的 5 起、G21 的 1 起），回報
`no_data_pre_2012`。這一批沒有附 PoleZero 檔；每個新資料夾放的是
既有資料夾中時段涵蓋該事件的 PZ 檔副本，因為
同名的 PZ 檔只差在 CREATED 那一行。）

這個模組處理的、已實際確認的怪癖（在哪裡處理見下面的
函式）：(1) PZ 資料夾名稱 <-> mseed 檔名的對應**不是**
可推導的規則——G9 的 PZ 資料夾是 `G9_0918`，但它的 mseed 檔是
`G09_0918_w.mseed`——所以下面的 SEISMIC_DATA_DIRS 是一張寫死的小表，
和這個程式庫 events.py 寫死而不推斷的先例
一致；(2) miniSEED 波形**不一定**都是 -60s/+600s——有些測站
（例如 G9 的 `ECS`，否則它是 csg 持續偏移案例極佳的 1.9km 同址
候選）是較短的觸發式強震儀紀錄，
在 +600s 之前很早就結束，用 `triggered_short_trace` 旗標處理，
而不是假設不存在；(3) 有些 (station, channel) 配對在同一個 mseed 檔裡有重複
波形（因缺口而分段的紀錄，或兩次相近的
觸發），用 `Stream.merge()` 處理，而不是天真的 `select()[0]`。

## 方法

計算兩個獨立的強震窗口偵測器，兩者都保留
（不默默偏好其中一個），因為已確認的短觸發波形
情況沒有足夠的震前基準讓 STA/LTA 運作：
  - `detect_window_sta_lta`：經典的 STA/LTA 觸發（需要震前
    基準；波形太短時回傳 `no_pre_event_baseline`，而不是捏造的
    結果）。
  - `detect_window_envelope_threshold`：帶通 + 包絡，以相對於
    自己峰值的門檻判斷——即使在短觸發波形上也能用，代價是
    對門檻敏感。
地磁那一側直接重用 `coseismic_step_analysis.py` 自己的去趨勢 +
step30 統計量機制（不用 `coseismic_stacking_analysis.py` 的
快取輸出——兩支新腳本沒有執行順序的相依），並和疊加
腳本完全一樣，以該事件自己的事件外雜訊底做 z-score
標準化。

用法：
  seismometer_comparison.py --self-test                         # PZ 解析 + 去除儀器響應的健全性檢查
  seismometer_comparison.py --group G9 --geomag-station csg      # 抽查代表性的 G9 案例
  seismometer_comparison.py --all                                 # 每個已下載地震資料的組
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from events import GROUPS, folder_events, get_group  # noqa: E402
from stat_utils import mad_zscore  # noqa: E402
from coseismic_step_analysis import (  # noqa: E402
    DETREND_WINDOW_SEC,
    EXCLUSION_BUFFER_SEC,
    _build_channels,
    _detrend,
    _load_station_days,
    _pos,
    _rank_stations_for_event,
    _step_statistic,
)

# 2026-08-17 重新整理：每個 GXX_MMDD 資料夾（mseed 檔已移進
# 資料夾裡，不再是 DATA_ROOT 底下零散的 <GROUP_MMDD>_w.mseed）現在放在
# 這個子目錄底下，而不是直接和 Gx 地磁
# 資料夾並列。
# 2026-09-14 搬家：seismometer/ 現在放在 geomag_precursor/ 本身裡面
# （和 Gx 一起，Gx 當天稍早也這樣搬過），不再在上一層的
# common.DATA_ROOT——所以下面用的是 GX_DATA_ROOT（== PROJECT_DIR），不是
# DATA_ROOT。
SEISMIC_ROOT = common.GX_DATA_ROOT / "seismometer"

# 寫死是因為資料夾名稱 <-> group-id 的對應無法可靠地
# 推導——G9 的不一致（pz_dir "G9_0918" vs mseed "G09_0918_w.mseed"）
# 就是具體的反例。G14 刻意不列：它的 2009-12-19
# 錨點早於地震資料來源 2012 年的起始。
#
# 以 (group_id, event.date) 為鍵——截至 2026-09-20，這涵蓋了當時登錄表 49 起事件中的
# 30 起（除了 G14、G21 和 13 起還沒下載的都有；不只是最初每組單一事件下載的 16 個
# 錨點），已確認存在於磁碟上，mseed 窗口也正確地
# 包住每個事件的發震秒。G10 的 2024-04-23a/2024-04-23b 共用
# 一個 PZ 資料夾（同一個 UTC 日曆天，測站中繼資料不會
# 每分鐘變），但各有自己的 mseed 檔——下面任何其他同組、同日曆天的事件
# 也用同樣的 a/b/c 尾碼 `date` 模式（原因見 events.py 的
# 2026-09-23 說明：`event.date` 在這裡以及
# coseismic_step_analysis.py/coseismic_stacking_analysis.py 的逐事件輸出命名中都被當成同組唯一鍵，
# 2026-09-23 那批 68 起新的非錨點事件——登錄表現在 117 起——否則
# 會在 G11/G13/G20 的 7 對同日事件上默默撞名；全部在
# events.py 源頭就加上尾碼，而不是在這裡繞過）。截至 2026-09-23，117 起事件中有 31 起
# 有真實地震資料（原有 30 起 + 2025-01-21b，從下面錨點已下載的窗口中
# 救回——見那個鍵自己的註解）。
SEISMIC_DATA_DIRS: dict[str, dict[str, dict]] = {
    "G1": {
        "2018-02-04": {"pz_dir": "G01_0204", "mseed": "G01_0204_w.mseed"},
        "2018-02-06": {"pz_dir": "G01_0206", "mseed": "G01_0206_w.mseed"},
    },
    # G2/G3 和 G6/G7/G8 在 2026-09-20 拆分前各是一個合併組；地震儀
    # 資料夾／檔名（用舊的合併名稱下載）維持不變。
    "G2": {
        "2019-04-18": {"pz_dir": "G02_G03_0418", "mseed": "G02_G03_0418_w.mseed"},
    },
    "G3": {
        "2019-08-08": {"pz_dir": "G02_G03_0808", "mseed": "G02_G03_0808_w.mseed"},
    },
    "G4": {
        "2020-12-10": {"pz_dir": "G04_1210", "mseed": "G04_1210_w.mseed"},
    },
    "G5": {
        "2021-04-18": {"pz_dir": "G05_0418", "mseed": "G05_0418_w.mseed"},
    },
    "G6": {
        "2021-10-24": {"pz_dir": "G06_G07_G08_1024", "mseed": "G06_G07_G08_1024_w.mseed"},
    },
    "G7": {
        "2022-01-03": {"pz_dir": "G06_G07_G08_0103", "mseed": "G06_G07_G08_0103_w.mseed"},
    },
    "G8": {
        "2022-03-23": {"pz_dir": "G06_G07_G08_0323", "mseed": "G06_G07_G08_0323_w.mseed"},
    },
    "G9": {
        "2022-09-17": {"pz_dir": "G09_0917", "mseed": "G09_0917_w.mseed"},
        "2022-09-18": {"pz_dir": "G9_0918", "mseed": "G09_0918_w.mseed"},
    },
    "G10": {
        "2024-04-03": {"pz_dir": "G10_0403", "mseed": "G10_0403_w.mseed"},
        "2024-04-23a": {"pz_dir": "G10_0422", "mseed": "G10_0422_w.mseed"},
        "2024-04-23b": {"pz_dir": "G10_0422", "mseed": "G10_0422_02_w.mseed"},
        "2024-05-10": {"pz_dir": "G10_0510", "mseed": "G10_0510_w.mseed"},
    },
    "G11": {
        "2025-01-21": {"pz_dir": "G11_0121", "mseed": "G11_0121_w.mseed"},
        # 2025-01-21b（當地 00:26:25 = 2025-01-20 16:26:25 UTC）落在錨點
        # 已下載的 event_utc-60s~+600s 窗口內（16:16:26~16:27:26 UTC）——同一個 mseed 檔，
        # 確實也涵蓋這個事件。2025-01-21c（當地 01:42:31）**沒有**（超出
        # 窗口 1 小時以上），這裡正確地不對應（回報 not_fetched）。
        "2025-01-21b": {"pz_dir": "G11_0121", "mseed": "G11_0121_w.mseed"},
    },
    "G12": {
        "2025-08-27": {"pz_dir": "G12_0827", "mseed": "G12_0827_w.mseed"},
    },
    "G13": {
        "2026-05-01": {"pz_dir": "G13_0501", "mseed": "G13_0501_w.mseed"},
    },
    "G15": {
        "2013-06-02": {"pz_dir": "G15_0602", "mseed": "G15_0602_w.mseed"},
    },
    "G16": {
        "2013-10-31": {"pz_dir": "G16_1031", "mseed": "G16_1031_w.mseed"},
    },
    "G17": {
        "2014-12-11": {"pz_dir": "G17_1211", "mseed": "G17_1211_w.mseed"},
        "2015-02-14": {"pz_dir": "G17_0214", "mseed": "G17_0214_w.mseed"},
    },
    "G18": {
        "2016-02-06": {"pz_dir": "G18_0206", "mseed": "G18_0206_w.mseed"},
        "2016-05-31": {"pz_dir": "G18_0531", "mseed": "G18_0531_w.mseed"},
    },
    "G19": {
        "2024-08-16": {"pz_dir": "G19_0816", "mseed": "G19_0816_w.mseed"},
    },
    "G20": {
        "2025-12-24": {"pz_dir": "G20_1224", "mseed": "G20_1224_w.mseed"},
        "2025-12-27": {"pz_dir": "G20_1227", "mseed": "G20_1227_w.mseed"},
    },
    "G22": {
        "2012-06-10": {"pz_dir": "G22_0610", "mseed": "G22_0610_w.mseed"},
    },
    "G23": {
        "2020-07-26": {"pz_dir": "G23_0726", "mseed": "G23_0726_w.mseed"},
    },
    "G24": {
        "2020-06-14": {"pz_dir": "G23_0614", "mseed": "G23_0614_w.mseed"},
    },
}

# 2026-09-25 起下載的東西都遵循同一個命名規則，所以用推導的
# 而不是逐一列出：G<2 位數組號>_<event.date 的 MMDD，也就是台灣當地
# 日期><如果有則加 a/b/c 尾碼>/，裡面放 <那個名稱>_w.mseed 加上
# 時段涵蓋該事件的 SAC PoleZero 檔。（上面人工列出的條目早於
# 這個規則：補零不一致、有些資料夾用 UTC 日期、還有拆分前的 G06_G07_G08 名稱。）
# 2026-09-25 那一批以這種方式加入 80 起事件（14 起 M>=6、66 起 M5）。
SEISMIC_DATA_SOURCE_START_UTC = pd.Timestamp("2012-01-01")  # 更早的都無法下載


def convention_dir_name(group_id: str, event_date: str) -> str:
    m = re.fullmatch(r"\d{4}-(\d{2})-(\d{2})([a-z]?)", event_date)
    return f"G{int(group_id[1:]):02d}_{m.group(1)}{m.group(2)}{m.group(3)}"


def _register_convention_dirs() -> None:
    for group_id, group in GROUPS.items():
        for event in group.events:
            if event.date in SEISMIC_DATA_DIRS.get(group_id, {}):
                continue
            name = convention_dir_name(group_id, event.date)
            if (SEISMIC_ROOT / name / f"{name}_w.mseed").exists():
                SEISMIC_DATA_DIRS.setdefault(group_id, {})[event.date] = {"pz_dir": name, "mseed": f"{name}_w.mseed"}


_register_convention_dirs()

GEOMAG_HALF_SEC = 240        # 目標：地磁那一側 z-score 化 step30 剖面的窗口
SEARCH_HALF_SEC = 180        # 目標：峰值搜尋子窗口，和 coseismic_step_analysis.py 的
                              # SCAN_HALF_SEC 相同慣例（2026-08-16 從 120 放寬到 180——之前 G12/G13/G20 的
                              # obs_lag 都卡在舊的 ±120 秒邊界。2026-08-19：研究過再放寬
                              # 到 300/360 秒，檢查 G20 的 2025-12-24 obs_lag=-178s（距
                              # 這個 180 秒邊界 2 秒）是不是截斷假象——見 plan
                              # artifact-wobbly-kettle.md。結果**不是**（300 秒時 obs_lag 仍是 -178s），
                              # 但較寬的窗口改變了其他好幾個**沒有**卡在邊界的事件
                              # 回報的 obs_lag（搜尋範圍越寬，即使只有雜訊，
                              # N 個中的最大值也越可能偶然變大；完整說明見 coseismic_step_analysis.py 的 SCAN_HALF_SEC 註解）
                              # ——退回 180/240 秒。兩者都是逐事件的**目標**，
                              # 不是固定值——見 `_effective_half_sec`：G10 的 2024-04-23a/b
                              # 只相隔 357 秒，所以它們實際的窗口被限制在這些目標以下，
                              # 避免一個事件的剖面／搜尋伸進另一個事件的真實異常。）
MIN_OFF_EVENT_SAMPLES = 200
MIN_POST_EVENT_SEC = 60      # 發震後資料少於這麼多的波形會被標記為 triggered_short_trace
PERSISTENCE_Z_THRESHOLD = 2.0
PERSISTENCE_TAIL_FRACTION = 0.3
ALIGNMENT_TOLERANCE_SEC = 5
# alignment_verdict_gated 的門檻：來自
# coseismic_step_analysis.py 的逐事件 step30 p 值（D7，2000 個隨機參考時間，同樣的 ±180 秒
# 搜尋）。沒有它的話，判定只是在對窗口中剛好最大的那個值
# 計時——2026-09-25 時 93 起事件中只有 7 起達到 p < 0.05。
ANOMALY_P_THRESHOLD = 0.05
STEP_SUMMARY_CSV = common.PROJECT_DIR / "data" / "interim" / "coseismic_step_analysis" / "summary.csv"
NOISE_PRE_LAGS = (-660, -60)   # 地磁 1Hz 差分雜訊比值的震前參考
# 磁力儀本身在強烈震動時可能會停：2025-01-21 大埔主震時的 twu
# （附近 2065 gal）在發震時變平，接下來 240 秒中有 146 秒
# 缺值——那裡的「沒有異常」是感測器中斷，不是平靜的磁場。
DROPOUT_WINDOW_SEC = (0, 300)
DROPOUT_MISSING_FRACTION = 0.2

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison"

_MSEED_CACHE: dict[str, object] = {}


# ---------------------------------------------------------------------------
# PZ 目錄：解析每個 PZ 檔自己的虛線註解檔頭（沒有另外的
# 檔名解析路徑——已確認檔頭本身就帶有所需的每一個
# 欄位）。
# ---------------------------------------------------------------------------

_PZ_FIELD_RE = {
    "station": re.compile(r"STATION\s*\(KSTNM\):\s*(\S+)"),
    "location": re.compile(r"LOCATION\s*\(KHOLE\):\s*(\S+)"),
    "channel": re.compile(r"CHANNEL\s*\(KCMPNM\):\s*(\S+)"),
    "start": re.compile(r"\bSTART\s*:\s*(\S+)"),
    "end": re.compile(r"\bEND\s*:\s*(\S+)"),
    "lat": re.compile(r"LATITUDE\s*:\s*([\-0-9.]+)"),
    "lon": re.compile(r"LONGITUDE\s*:\s*([\-0-9.]+)"),
    "elevation": re.compile(r"ELEVATION\s*:\s*([\-0-9.]+)"),
}
# 選用（PZ 檔解析不需要它）：每 m/s**2 的 counts——這裡每個 PZ
# 檔都是這個單位的 HL 頻段加速度儀（2026-09-25 檢查過）。
_PZ_SENSITIVITY_RE = re.compile(r"SENSITIVITY\s*:\s*([0-9.eE+\-]+)")


def _parse_pz_header(text: str) -> dict | None:
    vals = {}
    for key, pat in _PZ_FIELD_RE.items():
        m = pat.search(text)
        if m is None:
            return None
        vals[key] = m.group(1)
    return {
        "station": vals["station"], "location": vals["location"], "channel": vals["channel"],
        "start_utc": pd.Timestamp(vals["start"]), "end_utc": pd.Timestamp(vals["end"]),
        "lat": float(vals["lat"]), "lon": float(vals["lon"]), "elevation_m": float(vals["elevation"]),
        "sensitivity": (float(m.group(1)) if (m := _PZ_SENSITIVITY_RE.search(text)) else None),
    }


def _load_pz_catalog(pz_dir: Path) -> dict[tuple[str, str], list[dict]]:
    """key = (station, channel)。value = 時段 dict 的清單（location、
    start_utc、end_utc、lat、lon、elevation_m、path），每個 PZ 檔一個——
    一個測站／通道可以有不只一個時段（由
    `location` 及／或 `start_utc` 區分）。"""
    catalog: dict[tuple[str, str], list[dict]] = {}
    for path in sorted(pz_dir.glob("SAC_PZs_*")):
        text = path.read_text(errors="replace")
        h = _parse_pz_header(text)
        if h is None:
            continue
        h["path"] = path
        catalog.setdefault((h["station"], h["channel"]), []).append(h)
    return catalog


def _select_pz_epoch(catalog: dict, station: str, channel: str, location: str,
                      event_utc: pd.Timestamp) -> dict | None:
    """先比對 `location`（直接讀自 miniSEED 波形自己的
    stats.location——逐波形明確的中繼資料，不是猜的），再
    挑出在 event_utc 時有效（start_utc <= event_utc <= end_utc）
    且 start_utc 最晚的條目。只有在沒有任何有效時段符合 location 時，
    才退而忽略 location 比對（正常情況不應發生）。"""
    entries = catalog.get((station, channel), [])
    # 時段必須在事件發生時仍然有效，不只是已經開始——
    # 否則一個在事件前更換過儀器的測站會拿到它舊的響應
    in_force = [e for e in entries if e["start_utc"] <= event_utc <= e["end_utc"]]
    candidates = [e for e in in_force if e["location"] == location]
    if not candidates:
        candidates = in_force
    if not candidates:
        return None
    return max(candidates, key=lambda e: e["start_utc"])


def find_colocated_seismic_stations(pz_dir: Path, ref_lat: float, ref_lon: float, top_n: int = 3) -> list[dict]:
    """依距離 `ref_lat`/`ref_lon` 排序地震站——那是
    **地磁**測站自己的座標，不是震央。目標是
    和被檢定的磁力儀同址，而不是靠近
    地震震源。"""
    catalog = _load_pz_catalog(pz_dir)
    seen: dict[str, float] = {}
    for (station, _channel), entries in catalog.items():
        if station in seen:
            continue
        e = entries[0]
        seen[station] = common.haversine_km(ref_lat, ref_lon, e["lat"], e["lon"])
    ranked = sorted(seen.items(), key=lambda kv: kv[1])
    return [{"station": s, "distance_km": round(d, 2)} for s, d in ranked[:top_n]]


# ---------------------------------------------------------------------------
# miniSEED 存取
# ---------------------------------------------------------------------------

def _read_mseed_cached(mseed_path: Path):
    from obspy import read
    key = str(mseed_path)
    if key not in _MSEED_CACHE:
        _MSEED_CACHE[key] = read(str(mseed_path))
    return _MSEED_CACHE[key]


def load_trace(mseed_path: Path, pz_catalog: dict, station: str, channel: str,
                event_utc: pd.Timestamp, remove_response: bool = True) -> dict:
    """讀取（並快取）整個 mseed 檔，選出 (station, channel)
    波形，合併任何重複（已確認存在——有些 (station,
    channel) 配對在同一個檔案裡有 2 條波形，因缺口而分段的紀錄或
    兩次相近的觸發），需要時用對應的 PZ 時段去除儀器
    響應，並標記較短的觸發式強震儀波形，
    而不是假設窗口都一樣。"""
    from obspy import UTCDateTime

    st = _read_mseed_cached(mseed_path)
    sel = st.select(station=station, channel=channel).copy()
    if len(sel) == 0:
        return {"status": "no_trace"}
    try:
        sel.merge(method=1, fill_value=None)
    except Exception:
        pass

    ev = UTCDateTime(event_utc.isoformat())
    overlapping = [tr for tr in sel if tr.stats.starttime <= ev <= tr.stats.endtime]
    tr = (min(overlapping, key=lambda t: t.stats.endtime - t.stats.starttime)
          if overlapping else min(sel, key=lambda t: abs(t.stats.starttime - ev)))
    tr = tr.copy()

    triggered_short = (tr.stats.endtime - ev) < MIN_POST_EVENT_SEC
    response_removed = False
    if remove_response:
        from obspy.io.sac.sacpz import attach_paz
        epoch = _select_pz_epoch(pz_catalog, station, channel, tr.stats.location, event_utc)
        if epoch is not None:
            try:
                attach_paz(tr, str(epoch["path"]))
                tr.simulate(paz_remove=tr.stats.paz)
                response_removed = True
            except Exception as exc:
                return {"status": "response_removal_failed", "error": str(exc)}
        else:
            return {"status": "no_pz_epoch"}

    return {"status": "ok", "trace": tr, "response_removed": response_removed,
            "triggered_short_trace": bool(triggered_short)}


PGA_BAND_HZ = (0.1, 20.0)
PGA_WINDOW_SEC = (-10, 120)   # 只在事件自己的發震時間周圍這個窗口內取峰值


def peak_ground_acceleration(mseed_path: Path, pz_catalog: dict, station: str,
                              event_utc: pd.Timestamp, window_end_sec: float = PGA_WINDOW_SEC[1]) -> dict:
    """`station` 處的 PGA，單位 gal（cm/s**2）：兩個水平分量上的
    最大 |加速度|（只有兩個水平分量都不能用時才用 HLZ）。

    加速度 = 去平均、經 PGA_BAND_HZ 帶通的 counts / PZ 檔頭的
    SENSITIVITY（每 m/s**2 的 counts）——這些加速度儀標準的強震
    換算。峰值只在這個事件發震時間周圍 PGA_WINDOW_SEC[0] 到
    window_end_sec 之間取：有些事件共用另一個
    事件的 mseed 檔（G11 2025-01-21b 落在錨點窗口開始後 9 分鐘），
    不加窗的最大值會回傳另一個事件的震動。刻意不用 load_trace() 去除響應後的
    波形：PZ 檔的 INPUT UNIT 是 M，所以去除完整響應得到的是
    位移，不是加速度。"""
    from obspy import Stream, UTCDateTime

    ev = UTCDateTime(event_utc.isoformat())
    t0, t1 = ev + PGA_WINDOW_SEC[0], ev + window_end_sec
    per_comp: dict[str, float] = {}
    for comp in ("HLN", "HLE", "HLZ"):
        info = load_trace(mseed_path, pz_catalog, station, comp, event_utc, remove_response=False)
        if info["status"] != "ok":
            continue
        tr = info["trace"]
        epoch = _select_pz_epoch(pz_catalog, station, comp, tr.stats.location, event_utc)
        if epoch is None or not epoch.get("sensitivity"):
            continue
        peak = 0.0
        for seg in Stream([tr]).split():  # merge() 在有缺口的地方會留下遮罩陣列
            if seg.stats.npts < 2 * seg.stats.sampling_rate:
                continue
            seg.data = seg.data.astype(float)
            seg.detrend("demean")
            seg.filter("bandpass", freqmin=PGA_BAND_HZ[0],
                       freqmax=min(PGA_BAND_HZ[1], seg.stats.sampling_rate / 2 - 0.5), zerophase=True)
            win = seg.slice(t0, t1)  # 先濾整段再加窗，避免邊緣效應
            if win.stats.npts:
                peak = max(peak, float(np.max(np.abs(win.data))))
        if peak > 0:
            per_comp[comp] = peak / epoch["sensitivity"] * 100.0  # m/s**2 -> gal 換算
    horizontal = {c: v for c, v in per_comp.items() if c != "HLZ"}
    use = horizontal or per_comp
    if not use:
        return {"pga_gal": None, "pga_component": None}
    comp = max(use, key=use.get)
    return {"pga_gal": round(use[comp], 4), "pga_component": comp,
            "pga_by_component_gal": {c: round(v, 4) for c, v in per_comp.items()}}


# ---------------------------------------------------------------------------
# 強震窗口偵測器（兩者都計算，都不默默偏好
# ——原因見模組 docstring：已確認的短觸發波形情況
# 沒有足夠的震前基準讓 STA/LTA 運作）。
# ---------------------------------------------------------------------------

def detect_window_sta_lta(tr, event_utc: pd.Timestamp, sta_sec: float = 1.0, lta_sec: float = 10.0,
                           trigger_on: float = 3.5, trigger_off: float = 1.0) -> dict:
    from obspy import UTCDateTime
    from obspy.signal.trigger import classic_sta_lta, trigger_onset

    sr = tr.stats.sampling_rate
    ev = UTCDateTime(event_utc.isoformat())
    pre_event_sec = ev - tr.stats.starttime
    if pre_event_sec < lta_sec + 2:
        return {"status": "no_pre_event_baseline", "pre_event_sec": round(float(pre_event_sec), 2)}

    cft = classic_sta_lta(tr.data.astype(float), max(1, int(sta_sec * sr)), max(2, int(lta_sec * sr)))
    onsets = trigger_onset(cft, trigger_on, trigger_off)
    if len(onsets) == 0:
        return {"status": "no_trigger"}

    ev_idx = int(round((ev - tr.stats.starttime) * sr))
    on_i, off_i = min(onsets, key=lambda w: abs(w[0] - ev_idx))
    return {
        "status": "ok",
        "onset_lag_sec": round(float(on_i / sr - (ev - tr.stats.starttime)), 2),
        "offset_lag_sec": round(float(off_i / sr - (ev - tr.stats.starttime)), 2),
    }


def detect_window_envelope_threshold(tr, event_utc: pd.Timestamp, threshold_fraction: float = 0.1,
                                      band: tuple[float, float] = (1.0, 20.0)) -> dict:
    from obspy import UTCDateTime
    from obspy.signal.filter import envelope

    ev = UTCDateTime(event_utc.isoformat())
    tr2 = tr.copy()
    nyquist = tr.stats.sampling_rate / 2 - 0.5
    try:
        tr2.filter("bandpass", freqmin=band[0], freqmax=min(band[1], nyquist))
    except Exception:
        pass
    env = envelope(tr2.data.astype(float))
    peak = float(np.max(env)) if len(env) else 0.0
    if peak <= 0:
        return {"status": "flat_trace"}
    above = np.where(env > threshold_fraction * peak)[0]
    if len(above) == 0:
        return {"status": "no_signal_above_threshold"}

    sr = tr.stats.sampling_rate
    start_offset = tr.stats.starttime - ev  # 秒，波形起點相對於發震時間
    return {
        "status": "ok",
        "onset_lag_sec": round(float(above[0] / sr + start_offset), 2),
        "offset_lag_sec": round(float(above[-1] / sr + start_offset), 2),
        "peak_amplitude": peak,
    }


def _effective_half_sec(target_half_sec: int, event_utc: pd.Timestamp,
                         exclude_centers: list[pd.Timestamp]) -> int:
    """coseismic_step_analysis.py 同名輔助函式的本地副本
    （完整理由見那個模組）——把逐事件窗口限制在
    到同組最近*另一個*真實事件間隔的一半，讓
    G10 的 2024-04-23a/b（相隔 357 秒）永遠不會讓一個事件取出的
    剖面／搜尋伸進另一個事件的真實異常。在其他地方
    都不起作用（其他每組的事件都相隔數小時到數年）。"""
    others = [c for c in exclude_centers if c != event_utc]
    if not others:
        return target_half_sec
    nearest_gap_sec = min(abs((event_utc - c).total_seconds()) for c in others)
    return int(min(target_half_sec, nearest_gap_sec // 2))


# ---------------------------------------------------------------------------
# 地磁那一側：直接重用 coseismic_step_analysis.py 的去趨勢／統計量
# 機制（獨立計算，不依賴
# coseismic_stacking_analysis.py 的快取輸出）。
# ---------------------------------------------------------------------------

def geomag_profile(cfg: "common.GroupConfig", group, event, station: str, channel_type: str,
                    half_sec: int = GEOMAG_HALF_SEC) -> dict | None:
    event_utc = pd.Timestamp(event.time_utc)
    sibling_utcs = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]
    half_sec = _effective_half_sec(half_sec, event_utc, sibling_utcs)
    effective_search_half_sec = min(SEARCH_HALF_SEC, half_sec)
    buffer_sec = half_sec + EXCLUSION_BUFFER_SEC + DETREND_WINDOW_SEC // 2 + 60
    df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, buffer_sec)
    if df is None:
        return None
    channels = _build_channels(df)
    ch_label = "H" if channel_type == "XYZ" else "F"
    if ch_label not in channels:
        return None
    raw = channels[ch_label]
    idx = raw.index
    if idx.min() + pd.Timedelta(seconds=half_sec) > event_utc or \
       idx.max() - pd.Timedelta(seconds=half_sec) < event_utc:
        return None

    d = _detrend(raw)
    arr = _step_statistic(d, 30)  # step30：和 coseismic_stacking_analysis.py 相同的「主要」統計量，用於**起始**計時

    exclude_centers = sibling_utcs  # 和上面計算 half_sec 上限用的是同一個清單；這裡另取名稱，因為
                                      # 從這裡開始它的角色是虛無／基準排除，而不是決定窗口大小

    def _off_event_baseline(series: np.ndarray) -> tuple[float, float] | None:
        mask = np.ones(len(series), dtype=bool)
        for c in exclude_centers:
            p = _pos(idx, c)
            lo, hi = max(0, p - EXCLUSION_BUFFER_SEC), min(len(series), p + EXCLUSION_BUFFER_SEC + 1)
            mask[lo:hi] = False
        vals = series[mask]
        vals = vals[~np.isnan(vals)]
        if len(vals) < MIN_OFF_EVENT_SAMPLES:
            return None
        _, med, mad = mad_zscore(vals)
        if not np.isfinite(mad) or mad <= 1e-9:
            return None
        return med, mad

    def _extract_profile(series: np.ndarray) -> np.ndarray:
        p_center = _pos(idx, event_utc)
        lo_i, hi_i = max(0, p_center - half_sec), min(len(series), p_center + half_sec + 1)
        profile = np.full(2 * half_sec + 1, np.nan)
        dst_lo = lo_i - (p_center - half_sec)
        profile[dst_lo:dst_lo + (hi_i - lo_i)] = series[lo_i:hi_i]
        return profile

    step_base = _off_event_baseline(arr)
    if step_base is None:
        return None
    step_med, step_mad = step_base
    z_profile = (_extract_profile(arr) - step_med) / step_mad
    lags = np.arange(-half_sec, half_sec + 1)

    # `arr`（移動窗口的階躍**差值**）很適合找出轉變
    # **何時**發生，但在結構上不適合拿來
    # 問磁場之後是否**維持**偏移：一旦
    # 滑動窗口兩側都落在同一個新平台上，階躍差值
    # 統計量依定義就會回到約 0，不管那個
    # 平台是永久的（真正持續的偏移），還是本身就快要
    # 回復。持續性必須在去趨勢後的**磁場水準** `d`
    # 本身上判斷，用同樣方式做 z-score——這是和
    # 起始計時用的剖面不同的另一條剖面。
    level_base = _off_event_baseline(d)
    level_z_profile = (_extract_profile(d) - level_base[0]) / level_base[1] if level_base is not None else None

    # 寬窗口上的原始 1Hz 一階差分，給震動雜訊比值用，
    # compare_event 知道震動窗口後才計算（輸出前會移除）
    p0 = _pos(idx, event_utc)
    lo_d, hi_d = max(1, p0 + NOISE_PRE_LAGS[0]), min(len(raw), p0 + 601)
    raw_vals = raw.to_numpy(dtype=float)
    diff_lags = np.arange(lo_d, hi_d) - p0
    raw_diff = raw_vals[lo_d:hi_d] - raw_vals[lo_d - 1:hi_d - 1]

    search_mask = np.abs(lags) <= effective_search_half_sec
    sub, sub_lags = z_profile[search_mask], lags[search_mask]
    if np.all(np.isnan(sub)):
        obs_lag_sec, obs_peak_z = None, None
    else:
        j = int(np.nanargmax(np.abs(sub)))
        obs_lag_sec, obs_peak_z = int(sub_lags[j]), round(float(sub[j]), 4)

    return {
        "station": station, "statistic": "step30", "half_sec": half_sec,
        "search_half_sec": effective_search_half_sec,
        "lags_sec": lags.tolist(),
        "z_profile": [None if np.isnan(v) else round(float(v), 4) for v in z_profile],
        "level_z_profile": (None if level_z_profile is None else
                             [None if np.isnan(v) else round(float(v), 4) for v in level_z_profile]),
        "obs_lag_sec": obs_lag_sec, "obs_peak_z": obs_peak_z,
        "_diff_lags": diff_lags, "_diff": raw_diff,
    }


def geomag_dropout(diff_lags: np.ndarray, diff: np.ndarray) -> bool | None:
    """如果發震後 DROPOUT_WINDOW_SEC 內磁力儀的樣本
    缺漏超過 DROPOUT_MISSING_FRACTION 則為 True。"""
    sel = (diff_lags >= DROPOUT_WINDOW_SEC[0]) & (diff_lags <= DROPOUT_WINDOW_SEC[1])
    if not sel.any():
        return None
    return bool(np.mean(np.isnan(diff[sel])) > DROPOUT_MISSING_FRACTION)


def geomag_noise_ratio(diff_lags: np.ndarray, diff: np.ndarray,
                        onset: float | None, offset: float | None) -> float | None:
    """震動期間地磁 1Hz 一階差分的 RMS 除以
    事件前 10 分鐘（NOISE_PRE_LAGS）的 RMS——G10 試行中看到的
    震動雜訊特徵（約 11.6 倍）。震動窗口 = [onset, offset]，
    如果沒偵測到 offset，則為從 onset 起 60 秒。"""
    if onset is None:
        return None
    end = offset if offset is not None else onset + 60
    pre = diff[(diff_lags >= NOISE_PRE_LAGS[0]) & (diff_lags < NOISE_PRE_LAGS[1])]
    shake = diff[(diff_lags >= max(onset, 0)) & (diff_lags <= end)]
    pre, shake = pre[~np.isnan(pre)], shake[~np.isnan(shake)]
    if len(pre) < 300 or len(shake) < 5:
        return None
    pre_rms = float(np.sqrt(np.mean(pre ** 2)))
    return round(float(np.sqrt(np.mean(shake ** 2))) / pre_rms, 4) if pre_rms > 0 else None


# ---------------------------------------------------------------------------
# 比較／判定
# ---------------------------------------------------------------------------

def assess_persistence(lags_sec: list[int], z_profile: list[float | None],
                        shaking_offset_lag_sec: float | None,
                        z_threshold: float = PERSISTENCE_Z_THRESHOLD,
                        tail_frac: float = PERSISTENCE_TAIL_FRACTION) -> bool | None:
    """如果地磁 z 剖面在震動停止後（不知道震動結束時間時則為
    lag=0 之後），至少有 `tail_frac` 比例的樣本
    維持升高（|z| > z_threshold）則為 True——這是
    「G9 csg 不回復、G10 xcg 會回復」的量化形式。"""
    lags = np.array(lags_sec, dtype=float)
    z = np.array([np.nan if v is None else v for v in z_profile], dtype=float)
    ref = shaking_offset_lag_sec if shaking_offset_lag_sec is not None else 0.0
    tail = z[lags > ref]
    tail = tail[~np.isnan(tail)]
    if len(tail) == 0:
        return None
    return bool(np.mean(np.abs(tail) > z_threshold) >= tail_frac)


def alignment_verdict(geomag_lag_sec: float | None, shaking_onset_lag_sec: float | None,
                       shaking_offset_lag_sec: float | None, persists: bool | None,
                       tolerance_sec: float = ALIGNMENT_TOLERANCE_SEC) -> str:
    """以下之一：aligned_with_shaking / leads_shaking / persists_after_shaking_ends
    / insufficient_data——直接實作使用者自己的框架：
    對齊支持儀器雜訊的解釋；領先或
    比震動持續更久支持真實的地球物理機制。"""
    if geomag_lag_sec is None or shaking_onset_lag_sec is None:
        return "insufficient_data"
    if persists:
        return "persists_after_shaking_ends"
    if geomag_lag_sec < shaking_onset_lag_sec - tolerance_sec:
        return "leads_shaking"
    return "aligned_with_shaking"


def compare_event(cfg: "common.GroupConfig", group, event, geomag_station: str, channel_type: str,
                   seismic_station: str, mseed_path: Path, pz_catalog: dict) -> dict:
    event_utc = pd.Timestamp(event.time_utc)
    geomag = geomag_profile(cfg, group, event, geomag_station, channel_type)
    if geomag is None:
        return {"status": "no_geomag_data"}
    diff_lags, diff = geomag.pop("_diff_lags"), geomag.pop("_diff")

    seismic_channels: dict[str, dict] = {}
    best_comp = None
    for comp in ("HLZ", "HLN", "HLE"):
        info = load_trace(mseed_path, pz_catalog, seismic_station, comp, event_utc)
        if info["status"] != "ok":
            seismic_channels[comp] = {"status": info["status"]}
            continue
        sta_lta = detect_window_sta_lta(info["trace"], event_utc)
        env = detect_window_envelope_threshold(info["trace"], event_utc)
        seismic_channels[comp] = {
            "status": "ok",
            "response_removed": info["response_removed"],
            "triggered_short_trace": info["triggered_short_trace"],
            "trace_start_utc": str(info["trace"].stats.starttime),
            "trace_end_utc": str(info["trace"].stats.endtime),
            "sta_lta": sta_lta, "envelope": env,
        }
        if best_comp is None or comp == "HLZ":
            best_comp = comp

    if best_comp is None:
        return {"status": "no_seismic_trace", "geomag": geomag, "seismic_channels": seismic_channels}

    primary = seismic_channels[best_comp]
    shaking_onset = (primary["sta_lta"]["onset_lag_sec"] if primary["sta_lta"]["status"] == "ok"
                      else primary["envelope"].get("onset_lag_sec") if primary["envelope"]["status"] == "ok"
                      else None)
    shaking_offset = (primary["sta_lta"]["offset_lag_sec"] if primary["sta_lta"]["status"] == "ok"
                       else primary["envelope"].get("offset_lag_sec") if primary["envelope"]["status"] == "ok"
                       else None)

    persists = (assess_persistence(geomag["lags_sec"], geomag["level_z_profile"], shaking_offset)
                if geomag.get("level_z_profile") is not None else None)
    verdict = alignment_verdict(geomag["obs_lag_sec"], shaking_onset, shaking_offset, persists)
    others = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]
    pga = peak_ground_acceleration(mseed_path, pz_catalog, seismic_station, event_utc,
                                    _effective_half_sec(PGA_WINDOW_SEC[1], event_utc, others))

    return {
        "status": "ok",
        **pga,
        "geomag_noise_ratio": geomag_noise_ratio(diff_lags, diff, shaking_onset, shaking_offset),
        "geomag_dropout": geomag_dropout(diff_lags, diff),
        "primary_seismic_channel": best_comp,
        "geomag": geomag,
        "seismic_channels": seismic_channels,
        "shaking_onset_lag_sec": shaking_onset,
        "shaking_offset_lag_sec": shaking_offset,
        "geomag_persists_after_shaking": persists,
        "alignment_verdict": verdict,
    }


# ---------------------------------------------------------------------------
# 個案研究：compare_event 的結果加上降取樣的原始波形，給
# 報告畫圖用。重用 compare_event，而不是另外重算起始
# 偵測。
# ---------------------------------------------------------------------------

def _downsample(data: np.ndarray, factor: int) -> list[float]:
    if factor <= 1:
        return [round(float(v), 4) for v in data]
    n = (len(data) // factor) * factor
    if n == 0:
        return []
    return [round(float(v), 4) for v in data[:n].reshape(-1, factor).mean(axis=1)]


def build_case_study(group_id: str, event_date: str, geomag_station: str) -> dict:
    cfg = common.load_group_config(group_id)
    group = get_group(group_id)
    event = next((e for e in group.events if e.date == event_date), None)
    if event is None:
        return {"status": "unknown_event", "group": group_id, "event_date": event_date}

    if group_id not in SEISMIC_DATA_DIRS or event_date not in SEISMIC_DATA_DIRS[group_id]:
        return {"status": "no_seismic_data", "group": group_id, "event_date": event_date}

    channel_type = "XYZ" if cfg.xyz_pool.all_stations else "F"
    geomag_lat, geomag_lon = cfg.stations[geomag_station]["lat"], cfg.stations[geomag_station]["lon"]
    pz_dir = SEISMIC_ROOT / SEISMIC_DATA_DIRS[group_id][event_date]["pz_dir"]
    mseed_path = pz_dir / SEISMIC_DATA_DIRS[group_id][event_date]["mseed"]
    pz_catalog = _load_pz_catalog(pz_dir)

    colocated = find_colocated_seismic_stations(pz_dir, geomag_lat, geomag_lon, top_n=1)
    if not colocated:
        return {"status": "no_seismic_station", "group": group_id, "event_date": event_date}
    seismic_station = colocated[0]["station"]

    result = compare_event(cfg, group, event, geomag_station, channel_type, seismic_station, mseed_path, pz_catalog)
    result["group"] = group_id
    result["event_date"] = event_date
    result["geomag_station"] = geomag_station
    result["seismic_station"] = seismic_station
    result["seismic_colocation_km"] = colocated[0]["distance_km"]

    if result.get("status") != "ok":
        return result

    event_utc = pd.Timestamp(event.time_utc)
    primary_comp = result["primary_seismic_channel"]
    info = load_trace(mseed_path, pz_catalog, seismic_station, primary_comp, event_utc)
    waveform_payload = None
    if info["status"] == "ok":
        from obspy.signal.filter import envelope as _envelope
        tr = info["trace"]
        sr = tr.stats.sampling_rate
        t0_offset = (tr.stats.starttime.datetime - event_utc.to_pydatetime()).total_seconds()
        factor = max(1, int(round(sr / 20)))  # 降取樣到約 20Hz，讓 JSON 精簡
        wave_ds = _downsample(tr.data.astype(float), factor)
        env_ds = _downsample(_envelope(tr.data.astype(float)), factor)
        lags_ds = [round(t0_offset + i * factor / sr, 3) for i in range(len(wave_ds))]
        waveform_payload = {"sampling_rate_downsampled_hz": round(sr / factor, 3),
                             "lags_sec": lags_ds, "waveform": wave_ds, "envelope": env_ds}
    result["seismic_waveform_downsampled"] = waveform_payload
    return result


# ---------------------------------------------------------------------------
# 自我測試：在一個已知的真實檔案上做 PZ 解析 + 去除響應的來回檢查
# （不是合成資料——這裡值得健全性檢查的是 obspy/PZ
# 管線，而不是統計，統計已由 coseismic_stacking_analysis.py 的
# 自我測試涵蓋）。
# ---------------------------------------------------------------------------

def self_test() -> bool:
    ok = True
    pz_dir = SEISMIC_ROOT / "G10_0403"
    catalog = _load_pz_catalog(pz_dir)
    key = ("ALS", "HLZ")
    status1 = "PASS" if catalog.get(key) else "FAIL"
    print(f"[self-test] PZ catalog parse ALS/HLZ: {len(catalog.get(key, []))} epoch(s)  {status1}")
    ok = ok and status1 == "PASS"

    event_utc = pd.Timestamp("2024-04-02 23:58:11")
    mseed_path = pz_dir / "G10_0403_w.mseed"
    info = load_trace(mseed_path, catalog, "ALS", "HLZ", event_utc, remove_response=True)
    status2 = "PASS" if info.get("status") == "ok" and info.get("response_removed") else "FAIL"
    print(f"[self-test] load_trace ALS/HLZ status={info.get('status')} "
          f"response_removed={info.get('response_removed')}  {status2}")
    ok = ok and status2 == "PASS"

    if info.get("status") == "ok":
        data = info["trace"].data
        finite = bool(np.all(np.isfinite(data)))
        spread = float(np.std(data))
        status3 = "PASS" if finite and spread > 0 else "FAIL"
        print(f"[self-test] response-removed trace sanity: finite={finite} std={spread:.6g}  {status3}")
        ok = ok and status3 == "PASS"

        sta_lta = detect_window_sta_lta(info["trace"], event_utc)
        env = detect_window_envelope_threshold(info["trace"], event_utc)
        status4 = ("PASS" if sta_lta["status"] in ("ok", "no_trigger", "no_pre_event_baseline")
                   and env["status"] in ("ok", "no_signal_above_threshold", "flat_trace") else "FAIL")
        print(f"[self-test] detectors ran without error: sta_lta={sta_lta['status']} "
              f"envelope={env['status']}  {status4}")
        ok = ok and status4 == "PASS"

    colocated = find_colocated_seismic_stations(pz_dir, 24.038, 121.609, top_n=3)  # xcg 的座標
    status5 = "PASS" if colocated and colocated[0]["station"] == "HWA" else "FAIL"
    print(f"[self-test] find_colocated_seismic_stations(xcg): nearest={colocated[:1]}  {status5}")
    ok = ok and status5 == "PASS"

    return ok


# ---------------------------------------------------------------------------
# 涵蓋率摘要 + 真實資料統籌
# ---------------------------------------------------------------------------

def build_coverage_summary() -> dict:
    items = []
    for group_id, group in GROUPS.items():
        group_dirs = SEISMIC_DATA_DIRS.get(group_id, {})
        for event in group.events:
            if event.date in group_dirs:
                status = "available"
            elif pd.Timestamp(event.time_utc) < SEISMIC_DATA_SOURCE_START_UTC:
                status = "no_data_pre_2012"
            else:
                status = "not_fetched"
            items.append({"group": group_id, "date": event.date, "anchor": event.anchor,
                           "seismic_data_status": status})
    n_available = sum(1 for it in items if it["seismic_data_status"] == "available")
    return {"n_events_total": len(items), "n_available": n_available, "events": items}


def gate_verdicts(rows: list[dict]) -> list[dict]:
    """附上該事件自己的 step30 p 值（同一測站，H 或 F），來自
    coseismic_step_analysis.py，並推出 alignment_verdict_gated：異常顯著時
    沿用原本的判定，否則為 no_significant_anomaly
    （沒有對應的 p 值時為 anomaly_p_unavailable）。未經門檻篩選的
    alignment_verdict 維持原樣。"""
    step = pd.read_csv(STEP_SUMMARY_CSV) if STEP_SUMMARY_CSV.exists() else pd.DataFrame()
    if len(step):
        step = step[(step.statistic_type == "step") & (step.window_sec == 30)]
        lookup = {(r.group, r.event_date, r.station, r.channel): r.p_value for r in step.itertuples()}
    else:
        lookup = {}
    for row in rows:
        ch = "H" if row.get("channel_type") == "XYZ" else "F"
        p = lookup.get((row["group"], row["date"], row["geomag_station"], ch))
        p = None if p is None or p != p else float(p)
        row["geomag_step30_p"] = p
        row["anomaly_significant"] = None if p is None else bool(p < ANOMALY_P_THRESHOLD)
        verdict = row.get("alignment_verdict")
        if verdict in (None, "insufficient_data") or row.get("status") != "ok":
            row["alignment_verdict_gated"] = verdict
        elif p is None:
            row["alignment_verdict_gated"] = "anomaly_p_unavailable"
        else:
            row["alignment_verdict_gated"] = verdict if p < ANOMALY_P_THRESHOLD else "no_significant_anomaly"
    return rows


def run_available_events(group_ids: tuple[str, ...] | None = None) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)
    (OUT_DIR / "case_studies").mkdir(exist_ok=True)

    coverage = build_coverage_summary()
    (OUT_DIR / "coverage_summary.json").write_text(json.dumps(coverage, indent=2))

    targets = group_ids if group_ids else tuple(SEISMIC_DATA_DIRS.keys())
    rows: list[dict] = []
    run_summary = {"events": []}

    for group_id in targets:
        group_dirs = SEISMIC_DATA_DIRS.get(group_id, {})
        if not group_dirs:
            print(f"[{group_id}] no seismic data fetched, skipping", file=sys.stderr)
            continue
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)

        for event in group.events:
            if event.date not in group_dirs:
                continue  # 例如 G14 的事件（完全沒有地震資料來源涵蓋）
            pz_dir = SEISMIC_ROOT / group_dirs[event.date]["pz_dir"]
            mseed_path = pz_dir / group_dirs[event.date]["mseed"]
            pz_catalog = _load_pz_catalog(pz_dir)

            channel_type = "XYZ" if cfg.xyz_pool.all_stations else "F"
            near = _rank_stations_for_event(cfg, event, channel_type, 1)
            if not near:
                run_summary["events"].append({"group": group_id, "date": event.date, "status": "no_geomag_station"})
                continue
            geomag_station, geomag_distance_km = near[0]
            geomag_lat, geomag_lon = cfg.stations[geomag_station]["lat"], cfg.stations[geomag_station]["lon"]

            colocated = find_colocated_seismic_stations(pz_dir, geomag_lat, geomag_lon, top_n=3)
            if not colocated:
                run_summary["events"].append({"group": group_id, "date": event.date, "status": "no_seismic_station"})
                continue
            seismic_station = colocated[0]["station"]

            result = compare_event(cfg, group, event, geomag_station, channel_type,
                                    seismic_station, mseed_path, pz_catalog)
            result.update({"group": group_id, "geomag_station": geomag_station,
                            "geomag_distance_km": round(geomag_distance_km, 1),
                            "seismic_station": seismic_station,
                            "seismic_colocation_km": colocated[0]["distance_km"]})

            out_path = OUT_DIR / "events" / f"{group_id}__{event.date}__{geomag_station}.json"
            out_path.write_text(json.dumps(result, indent=2, default=str))

            rows.append({
                "group": group_id, "date": event.date, "anchor": event.anchor,
                "geomag_station": geomag_station,
                "geomag_distance_km": round(geomag_distance_km, 1),
                "seismic_station": seismic_station, "seismic_colocation_km": colocated[0]["distance_km"],
                "status": result.get("status"),
                "geomag_obs_lag_sec": (result.get("geomag") or {}).get("obs_lag_sec"),
                "shaking_onset_lag_sec": result.get("shaking_onset_lag_sec"),
                "shaking_offset_lag_sec": result.get("shaking_offset_lag_sec"),
                "geomag_persists_after_shaking": result.get("geomag_persists_after_shaking"),
                "alignment_verdict": result.get("alignment_verdict"),
                "geomag_obs_peak_z": (result.get("geomag") or {}).get("obs_peak_z"),
                "geomag_noise_ratio": result.get("geomag_noise_ratio"),
                "geomag_dropout": result.get("geomag_dropout"),
                "pga_gal": result.get("pga_gal"), "pga_component": result.get("pga_component"),
                "channel_type": channel_type,
            })
            run_summary["events"].append({"group": group_id, "date": event.date, "status": result.get("status"),
                                           "alignment_verdict": result.get("alignment_verdict")})
            print(f"[{group_id} {event.date}] geomag={geomag_station} seismic={seismic_station} "
                  f"status={result.get('status')} verdict={result.get('alignment_verdict')}", file=sys.stderr)

    pd.DataFrame(gate_verdicts(rows)).to_csv(OUT_DIR / "comparison_summary.csv", index=False)
    (OUT_DIR / "all_comparisons_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    print(f"[run] {len(rows)} events compared -> {OUT_DIR}", file=sys.stderr)

    if "G9" in targets:
        cs = build_case_study("G9", "2022-09-18", "csg")
        (OUT_DIR / "case_studies" / "G9_csg_vs_ECS.json").write_text(json.dumps(cs, indent=2, default=str))
        print(f"[case-study] G9 csg -> verdict={cs.get('alignment_verdict')}", file=sys.stderr)
    if "G10" in targets:
        cs = build_case_study("G10", "2024-04-03", "xcg")
        (OUT_DIR / "case_studies" / "G10_xcg_vs_HWA.json").write_text(json.dumps(cs, indent=2, default=str))
        print(f"[case-study] G10 xcg -> verdict={cs.get('alignment_verdict')}", file=sys.stderr)

    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="PZ-parse + response-removal sanity check only")
    ap.add_argument("--group", action="append", dest="groups",
                     help="restrict to this group ID (repeatable); default is every group with seismic data fetched")
    ap.add_argument("--geomag-station", dest="geomag_station", default=None,
                     help="(informational, single-group spot checks) which geomag station to compare against")
    ap.add_argument("--all", action="store_true", help="explicit alias for the default (no --group filter)")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    group_ids = tuple(args.groups) if args.groups else None
    run_available_events(group_ids)
