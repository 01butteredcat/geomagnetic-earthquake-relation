"""同震（瞬間、發震時刻）地磁階躍／突波偵測器。

這個程式庫中其他所有東西（`build_daily_features.py`、`ulf_analysis.py`、
`compute_indices.py`、`superposed_epoch_analysis.py`……）都在日或
重取樣成 1 分鐘的粒度上運作，問的是「地震前後幾天到幾週
有沒有前兆／餘波特徵」。這支腳本問的是一個
更窄的問題：**就在地震發震那一秒，原始 1Hz X/Y/Z/F 磁場
有沒有階躍型或突波型的跳動**，也就是壓磁效應或
地震波感應機制原則上可能產生的那種
訊號（例如 2011 年東北 Mw9.0 有近場測站出現幾 nT 同震
脈衝的報告）——或者，同樣可能的是，震動造成的
感測器外殼雜訊，而不是真正的磁場變化（見下面「G10
試行發現了什麼」）。

## 沿革／範圍

第一階段（單一代表事件試行）：只有 G10 的錨點事件（台灣時間 2024-04-03 07:58:11
／2024-04-02 23:58:11 UTC，M7.2），對照它最近的兩個 XYZ 測站。
第二階段（目前）：透過 `run_all()` / `--all` 擴展到 `events.py` 117 起事件
登錄表中的每一個事件——到底改了什麼、為什麼改，見下面
「擴展到全部 117 起事件」。

### G10 試行發現了什麼（這說明了為什麼在下任何
### 通則性結論之前，必須在不只一個事件上執行）

G10 最近的兩個測站（xcg 24.8km、cnu 66.3km）都在發震時間
約一分鐘內出現統計顯著的階躍／突波異常
（好幾個 p <= 0.0025，在那個規模下撐過 Bonferroni 校正）。但
直接檢視原始波形後發現，這*不是*乾淨的直流階躍：
逐樣本雜訊振幅跳到震前
基準的約 11.6-11.7 倍，近站比遠站
晚幾秒開始——和強烈地動抵達、實際
搖晃磁力儀外殼一致，而不是壓磁型的磁場變化。
這種「震動雜訊」型態是不是整個資料集的普遍情況，
還是 G10（最大、儀器最完整的事件）只是特例，
正是 117 起事件執行要回答的。

## 為什麼需要新的機制

`parser.parse_day_file` 和 `common.resolve_day_ref` 只處理完整的
測站–日——這個程式庫中任何地方都沒有子範圍讀取，也沒有跨午夜拼接
（每個既有的呼叫者一次只處理一個日曆
天）。這支腳本補上了這點：它載入涵蓋發震時間戳所需緩衝區的
完整一天（或幾天）原始資料（對於像 G10 這樣
接近 UTC 午夜的發震時間，就是兩個測站–日），把它們串接起來，
直接在得到的連續 1Hz 序列上運作——沒有
計算上的理由要先裁成子窗口，因為 `parse_day_file`
不管之後用到多少，都會解析整個 86400 列的
檔案。

## 統計設計（仿照 `surrogate_test.py` / `superposed_epoch_
analysis.py` 已建立的模式：觀測統計量 + 重抽建立的
虛無分布 + p 值 + 固定種子，而不是發明新的
方法）

1. 每個通道用和 `ulf_analysis.py::_detrend` 相同的 1 小時置中滑動平均
   相減去趨勢（在這裡保留一份本地副本，
   並註明如此，而不是從一支目的是 ULF 近站／遠站流程的
   腳本匯入——為什麼 `ulf_analysis.py` 完全沒被修改，見
   計畫檔）。
2. 在去趨勢序列 `d` 上算兩個互補的統計量：
   - **step**：`S_M(t) = mean(d[t+1 .. t+M]) - mean(d[t-M+1 .. t])`，
     `M in {10, 30, 90}` 秒，一個滑動的階躍／偏移偵測器。
   - **spike**：`|diff(d)[t]|`，單一樣本的跳動偵測器（和
     `build_daily_features.py::_despike` 的門檻同樣的想法，這裡把它
     本身當成統計量，而不是雜訊濾波器）。
3. **觀測值**：回報發震秒 ±SCAN_HALF_SEC（180 秒，逐事件
   設上限——見 `_effective_half_sec`）內的 `max(|statistic|)`（一個
   搜尋窗口，不是單一點——容許回報時間的捨入
   和任何真實的傳播延遲）。
4. **虛無分布**：在同一個已載入窗口中其他地方（距真實發震時間
   >=600 秒）抽出的 N_NULL=2000 個隨機參考時間上，重複*同樣的*
   「搜尋 ±SCAN_HALF_SEC 並取極值」程序，
   讓虛無分布和觀測值在對等的基礎上比較，
   而不是用單一點的虛無分布。
5. p 值 = (1 + #{|null| >= |obs|}) / (N_NULL + 1)，單尾，因為兩個
   統計量本來就取了 |.|。
6. 每個結果也回報 `approx_min_detectable_nt`（虛無
   分布的第 95 百分位數）——讓不顯著的結果可以
   誠實地解讀成「檢定力足夠、沒有訊號」（這個數字小）
   或「這個測站–日的雜訊底太高，什麼都說不了」（這個
   數字大），而不是把兩者都混成「不顯著」。

## 已知限制（事先說明，不是事後的藉口）

- 物理上，同震地磁瞬變（壓磁效應、
  地震波感應電流）在文獻中真正有充分紀錄的，只有
  非常大（M8+）事件的近場測站（個位數 km）。G10 的 M7.2 主震和它
  最近約 10-30km 的
  測站遠在那個範圍之外——**這裡得到虛無結果是
  預期中的結果**，不代表這個機制在別處不存在。
- 單一事件、2 個測站、總共大約 30 個檢定（見輸出中的
  `multiple_comparisons_context`）——沒有做多重比較
  校正（這個規模不需要），但會報告檢定數，
  以便和未來的多事件擴展比較。
- 刻意不排除地磁暴：磁暴主相／恢復相
  在數小時到數天的尺度上變化，大多會被 1 小時去趨勢去掉，而且
  虛無分布是從同一個已載入窗口抽的，所以被磁暴抬高的雜訊
  底會讓觀測和虛無一起升高。剩下的風險是搜尋窗口內有單一次
  外部瞬變（SSC、亞暴開始、Pi2）。
  改成標記每個事件（`is_storm_day` / `is_storm_onset`，取自該
  組的 `storm_days.csv`，見 `_storm_status`），`--all` 會在
  `all_events_run_summary.json::storm_sensitivity` 中報告磁暴 vs. 平靜的拆分。
  直接排除被標記的事件不切實際：依日尺度
  流程的 Kp>=5 / Dst<=-30nT（+恢復期）定義，大約一半的
  事件都落在被標記的日子。
- 階躍統計量「後減前」的實作（見
  `_step_statistic`）有一個已知約 1 個樣本的邊緣對齊近似，來自
  計算前瞻平均所用的反向滾動技巧；這在這裡使用的
  ±SCAN_HALF_SEC（180 秒）搜尋窗口尺度下無關緊要。

## 擴展到全部 117 起事件（已完成——`run_all()` / `--all`）

- 直接使用 `events.py` 的 `GROUPS` dict 和每個 `Event.time_utc`（**不是**
  `catalog_utils.load_extended_events()`，那會截斷到只有日期的
  解析度）——每一組的每一個事件，不只 20 個 `anchor=True`
  事件，因為每個事件的同震檢定都是它自己獨立的
  量測（不像日尺度的跨組統計，這裡沒有
  前震／主震／餘震偽重複的疑慮）。
- 測站距離只相對於每組單一錨點事件預先算好
  （`common.py::_discover_stations`）；`_rank_stations_for_event`
  為每個個別事件重新計算 `common.haversine_km(station_lat, station_lon, ev.lat, ev.lon)`，
  讓非錨點事件得到正確的
  近站排名。
- 當一組在兩種測站池都有測站時，每個事件都會對**兩者**——只有 F 和 XYZ 測站池——
  做檢定（例如 G2_G3 有 10 個 F 測站和 2 個
  XYZ 測站——兩者都檢定，不會犧牲其中一個）；大多數
  組別只有一個非空的測站池（見 `common.py` 逐組的
  F/XYZ 測站池列印），所以大多就變成「有哪個測站池就用哪個」。
- 虛無候選排除涵蓋同一組中的每一個事件（不只
  被檢定的那個）在 `EXCLUSION_BUFFER_SEC` 以內的範圍，因為有些組
  的事件在時間上夠近（G9 的兩起事件相隔 17 小時，G10 的
  2024-04-23a/b 相隔六分鐘），天真的單一事件排除
  可能讓虛無分布被落在同一個已載入窗口中的另一起真實
  地震汙染。
- *觀測值*的搜尋窗口和虛無那一側有同樣的汙染風險，
  所以逐事件設上限來防範（`_effective_half_sec`）：
  G10 的 2024-04-23a/b 只相隔 357 秒，所以一旦 SCAN_HALF_SEC 超過
  約 178 秒，在其中一個事件周圍不設上限的 ±SCAN_HALF_SEC 搜尋就會越過
  中點，伸進另一個事件自己的異常。`_effective_half_sec`
  把每個事件實際的搜尋窗口限制在 `min(SCAN_HALF_SEC, 同組最近
  事件間隔 // 2)`，所以只有這一對會被縮窄（其他
  每組的事件都相隔數小時到數年）；限制後的值
  逐事件記錄在 `scan_half_sec`，以求透明。
- 原本 `_extremum_in_window` 中 O(n) 的布林遮罩掃描（每次虛無抽樣
  都遮罩整個約 172,800 個樣本的陣列）是單一事件試行的
  執行時間瓶頸（2 個測站約 20 分鐘）；改成 O(window) 的
  位置索引切片（`_pos` + 現在的 `_extremum_in_window`），
  因為一天的資料保證是完整的 86400 列網格，
  串接的相鄰日之間沒有缺口，所以對陣列自己起點
  做時間戳運算可以得到精確的整數陣列位置。
- 正向對照注入測試**沒有**在全部 117 起事件上重複
  （見計畫檔）——方法已經在 G10 上驗證過一次
  （`injection_power_curve.json`，`run_all()` 不會動它）；`--all` 一律
  以 `run_injection=False` 執行。
- 跨事件疊加（用
  `stat_utils.mad_zscore` 以每個事件自己的事件外統計量值標準化它的窗口，再像
  `superposed_epoch_analysis.py` 疊加日尺度序列那樣跨事件平均，
  並以每個事件的隨機參考
  時間建立虛無帶）是自然的下一步，如果 117 起事件的執行顯示許多
  個別檢定力不足、但有提示性的事件；但這一輪刻意不
  實作——見計畫檔。

用法：
  coseismic_step_analysis.py                  # 只跑 G10 錨點、兩種測站池、完整執行（冒煙測試）
  coseismic_step_analysis.py --all             # events.py 登錄表中全部 117 起事件
  coseismic_step_analysis.py --self-test       # 只做合成白雜訊健全性檢查
  coseismic_step_analysis.py --no-injection    # 跳過正向對照注入測試（只限 G10 路徑）
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
import parser as sec_parser  # noqa: E402
from events import GROUPS, folder_events, get_group  # noqa: E402
from stat_utils import histogram_summary  # noqa: E402

SEED = 20260805


def keyed_rng(*keys) -> np.random.Generator:
    """每個鍵 tuple 一條獨立、穩定的隨機亂數流，例如 (group, event,
    station, channel)。所有事件共用一條亂數流時，每次虛無
    抽樣都取決於前面事件消耗了多少次抽樣，所以修改
    任何一個事件（或只跑單一事件而不是 --all）都會重洗
    之後每個事件的 p 值。"""
    return np.random.default_rng([SEED, zlib.crc32("|".join(map(str, keys)).encode())])

GROUP_ID = "G10"
N_STATIONS = 2

STEP_WINDOWS_SEC = (10, 30, 90)
SCAN_HALF_SEC = 180          # 目標：在回報發震秒的兩側各搜尋這麼遠
                              # （2026-08-16 從 120 放寬到 180：之前 G12/G13/G20 的 obs_lag 都卡在
                              # 舊的 ±120 秒邊界。2026-08-19：研究過再放寬到 300 秒，
                              # 檢查 G20 的 2025-12-24 obs_lag=-178s（距這個 180 秒
                              # 邊界 2 秒）本身是不是截斷假象——見 plan
                              # artifact-wobbly-kettle.md。結果**不是**：用 300 秒重跑，G20 的
                              # obs_lag 仍是 -178s，確認那是真正的局部極值，
                              # 不是邊界截斷。但 300 秒的嘗試也*改變*了其他好幾個
                              # **沒有**卡在邊界的事件回報的 obs_lag（例如 G1
                              # 2018-02-06：180 秒時 -133s -> 300 秒時 243s，z=2.84 -> z=3.43）——
                              # 搜尋窗口越寬，即使只有雜訊，N 個中的最大值也越可能偶然變大
                              # （虛無分布也會相應變寬，所以 p 值仍然有效，
                              # 但回報成「那個」峰值的延遲會變得取決於窗口大小，
                              # 離發震秒越遠，物理上越難
                              # 解讀）。既然 G20 終究不需要較寬的窗口，就退回
                              # 180 秒。這是逐事件的**目標**，不是固定值——見
                              # `_effective_half_sec`，它從 300 秒嘗試中保留下來，是真正（雖然範圍很窄）的
                              # 修正：G10 的 2024-04-23a/b 只相隔 357 秒，所以不管這個目標，它們的
                              # 搜尋窗口都被限制在約 178 秒，避免一個事件的搜尋
                              # 伸進另一個事件的真實異常。）
DETREND_WINDOW_SEC = 3601    # 和 ulf_analysis.py::_detrend 相同（1 小時置中滑動平均）
REQUIRED_MARGIN_SEC = SCAN_HALF_SEC + max(STEP_WINDOWS_SEC)  # 270 秒：兩端各需要的最少資料支撐
                              # （用的是 SCAN_HALF_SEC *目標值*，不是逐事件限制後的值，所以這個
                              # 餘裕檢查不管有沒有逐事件上限，都保持保守且一致）
N_NULL = 2000
EXCLUSION_BUFFER_SEC = 600   # 虛無參考時間要離真實發震時間這麼遠
INJECTION_EXCLUSION_BUFFER_SEC = 2000  # 較寬，給正向對照執行用（見 run_injection_power_curve）
INJECTION_AMPLITUDES_SIGMA = (1, 2, 3, 5, 10)  # 局部虛無 MAD 的倍數
MISSING_FRACTION_THRESHOLD = 0.10  # ±SCAN_HALF_SEC 窗口內 NaN 比例超過這個，該檢定作廢

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_step_analysis"


# ---------------------------------------------------------------------------
# 訊號處理（本地副本／改寫——為什麼這些不從
# ulf_analysis.py / build_daily_features.py 匯入，見模組 docstring）
# ---------------------------------------------------------------------------

def _detrend(x: pd.Series) -> np.ndarray:
    """和 ulf_analysis.py::_detrend 相同的公式：減掉 1 小時置中
    滑動平均。本地副本，不是匯入——見模組 docstring。"""
    trend = x.rolling(window=DETREND_WINDOW_SEC, center=True, min_periods=1).mean()
    return (x - trend).to_numpy(dtype=float)


def _step_statistic(d: np.ndarray, half_win: int) -> np.ndarray:
    """S_M(t) = mean(d[t+1..t+M]) - mean(d[t-M+1..t])。前瞻平均所用
    反向滾動技巧約 1 個樣本的邊緣對齊注意事項，見模組 docstring；
    在這支腳本使用的 ±SCAN_HALF_SEC
    （180 秒）搜尋窗口尺度下無關緊要。"""
    s = pd.Series(d)
    min_p = max(3, int(0.6 * half_win))
    pre = s.rolling(half_win, min_periods=min_p).mean()
    post = (
        s[::-1].reset_index(drop=True)
        .rolling(half_win, min_periods=min_p).mean()
        [::-1].reset_index(drop=True).shift(-1)
    )
    post.index = s.index
    return (post - pre).to_numpy(dtype=float)


def _spike_statistic(d: np.ndarray) -> np.ndarray:
    """|diff(d)[t]|，和 build_daily_features.py::
    _despike 同樣的跳動偵測想法，這裡本身當成統計量，而不是雜訊濾波器。"""
    return np.abs(np.diff(d, prepend=d[0]))


# ---------------------------------------------------------------------------
# 資料載入
# ---------------------------------------------------------------------------

def _load_station_days(gdms_dir: Path, station: str, event_utc: pd.Timestamp,
                        buffer_sec: int) -> tuple[pd.DataFrame | None, list[str], list[str]]:
    """載入並串接涵蓋
    [event_utc - buffer_sec, event_utc + buffer_sec] 的測站–日檔案。parser.py 沒有子範圍
    讀取，所以每個涵蓋到的日子都完整解析（86400 列），
    結果再由呼叫者切片。回傳
    (concatenated_df_or_None, missing_date_strs, wanted_date_strs)。"""
    lo = event_utc - pd.Timedelta(seconds=buffer_sec)
    hi = event_utc + pd.Timedelta(seconds=buffer_sec)
    wanted = sorted({d.strftime("%Y%m%d") for d in pd.date_range(lo.normalize(), hi.normalize(), freq="D")})
    frames = []
    missing = []
    for date_str in wanted:
        ref = common.resolve_day_ref(gdms_dir, station, date_str)
        if ref is None:
            missing.append(date_str)
            continue
        frames.append(sec_parser.parse_day_file(ref, station))
    if not frames:
        return None, missing, wanted
    df = pd.concat(frames).sort_index()
    return df, missing, wanted


# ---------------------------------------------------------------------------
# 窗口內極值的輔助函式（觀測值和虛無抽樣共用，所以
# 比較的兩側都經過完全相同的程式碼）。用位置
# （整數秒偏移）切片，而不是對整個陣列做布林遮罩：
# 一天的資料保證是完整的 86400 列網格，串接的
# 相鄰日之間沒有缺口，所以 `index[0]` + 整數秒可以得到精確的
# 陣列位置——每次呼叫 O(window)，而不是 O(len(index))。
# 原本的布林遮罩版本是這支腳本在 117 起事件規模下的
# 執行時間瓶頸（見模組 docstring 的「擴展到全部 117 起事件」）。
# ---------------------------------------------------------------------------

def _pos(index: pd.DatetimeIndex, t: pd.Timestamp) -> int:
    return int(round((t - index[0]).total_seconds()))


def _extremum_in_window(stat_values: np.ndarray, index: pd.DatetimeIndex,
                         center: pd.Timestamp, half_sec: int) -> tuple[float, int] | None:
    p = _pos(index, center)
    lo, hi = max(0, p - half_sec), min(len(stat_values), p + half_sec + 1)
    sub = stat_values[lo:hi]
    if len(sub) == 0 or np.all(np.isnan(sub)):
        return None
    j = int(np.nanargmax(np.abs(sub)))
    return float(sub[j]), int((lo + j) - p)


def _draw_null_centers(rng: np.random.Generator, valid_lo: pd.Timestamp, valid_hi: pd.Timestamp,
                        exclude_centers: list[pd.Timestamp], exclude_buffer_sec: int, n: int) -> list[pd.Timestamp]:
    """exclude_centers 是清單（不是單一時間戳），所以候選時間只要
    落在被檢定組中*任何*真實事件附近就會被拒絕——
    需要這樣是因為有些組的事件在時間上夠近（見模組
    docstring），單一事件排除可能讓虛無抽樣落在
    另一起真實地震的正上方。"""
    total_sec = (valid_hi - valid_lo).total_seconds()
    centers: list[pd.Timestamp] = []
    attempts = 0
    while len(centers) < n and attempts < n * 5:
        attempts += 1
        cand = valid_lo + pd.Timedelta(seconds=rng.uniform(0, total_sec))
        if any(abs((cand - ec).total_seconds()) < exclude_buffer_sec for ec in exclude_centers):
            continue
        centers.append(cand)
    return centers


def _effective_half_sec(target_half_sec: int, event_utc: pd.Timestamp,
                         exclude_centers: list[pd.Timestamp]) -> int:
    """限制搜尋／掃描半窗口，讓它永遠不會越過到同組最近
    *另一個*真實事件的中點。需要這樣是因為
    SCAN_HALF_SEC 是全部 117 起事件共用的目標，但 G10 的
    2024-04-23a/b 只相隔 357 秒——一旦目標超過約 178 秒，在其中一個周圍
    不設上限的 ±SCAN_HALF_SEC 搜尋就會伸進
    另一個自己的真實異常，汙染「觀測值」的搜尋
    （虛無那一側已經透過 `_draw_null_centers` 中的 exclude_centers
    處理；這是觀測那一側對應的防護）。
    其他每組的事件都相隔數小時到數年，所以除了那一對之外，
    這在任何地方都不起作用。"""
    others = [c for c in exclude_centers if c != event_utc]
    if not others:
        return target_half_sec
    nearest_gap_sec = min(abs((event_utc - c).total_seconds()) for c in others)
    return int(min(target_half_sec, nearest_gap_sec // 2))


# ---------------------------------------------------------------------------
# 逐（測站, 通道）檢定
# ---------------------------------------------------------------------------

def _test_channel(raw: pd.Series, event_utc: pd.Timestamp, exclude_centers: list[pd.Timestamp],
                   rng: np.random.Generator, channel_label: str, scan_half_sec: int) -> dict:
    idx = raw.index
    nan_frac_at_event = float(
        raw[(idx >= event_utc - pd.Timedelta(seconds=scan_half_sec)) &
            (idx <= event_utc + pd.Timedelta(seconds=scan_half_sec))].isna().mean()
    ) if len(raw) else 1.0

    if idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) > event_utc or \
       idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) < event_utc:
        return {"channel": channel_label, "data_status": "insufficient_margin"}
    if nan_frac_at_event > MISSING_FRACTION_THRESHOLD:
        return {"channel": channel_label, "data_status": "insufficient_data",
                "nan_fraction_at_event": round(nan_frac_at_event, 4)}

    d = _detrend(raw)
    valid_lo = idx.min() + pd.Timedelta(seconds=scan_half_sec)
    valid_hi = idx.max() - pd.Timedelta(seconds=scan_half_sec)
    null_centers = _draw_null_centers(rng, valid_lo, valid_hi, exclude_centers, EXCLUSION_BUFFER_SEC, N_NULL)

    out = {"channel": channel_label, "data_status": "ok",
           "nan_fraction_at_event": round(nan_frac_at_event, 4), "scan_half_sec": scan_half_sec, "step": {}}

    for m in STEP_WINDOWS_SEC:
        stat = _step_statistic(d, m)
        obs = _extremum_in_window(stat, idx, event_utc, scan_half_sec)
        null_vals = [r for c in null_centers if (r := _extremum_in_window(stat, idx, c, scan_half_sec)) is not None]
        out["step"][str(m)] = _summarize_test(obs, null_vals)

    spike_stat = _spike_statistic(d)
    obs = _extremum_in_window(spike_stat, idx, event_utc, scan_half_sec)
    null_vals = [r for c in null_centers if (r := _extremum_in_window(spike_stat, idx, c, scan_half_sec)) is not None]
    out["spike"] = _summarize_test(obs, null_vals)

    return out


def _summarize_test(obs: tuple[float, int] | None, null_vals: list[tuple[float, int]]) -> dict:
    if obs is None or not null_vals:
        return {"data_status": "no_observed_or_null_value"}
    obs_val, obs_lag = obs
    null_abs = np.array([abs(v) for v, _ in null_vals])
    n_used = len(null_abs)
    p_value = (1 + int(np.sum(null_abs >= abs(obs_val)))) / (n_used + 1)
    return {
        "obs_value_nt": round(obs_val, 4),
        "obs_lag_sec": obs_lag,
        "p_value": round(p_value, 5),
        "n_null_samples": n_used,
        "approx_min_detectable_nt": round(float(np.percentile(null_abs, 95)), 4),
        "null_summary": histogram_summary(null_abs),
    }


# ---------------------------------------------------------------------------
# 正向對照：把合成脈衝注入真實背景雜訊，
# 確認偵測器的 p 值隨注入振幅單調下降。
# ---------------------------------------------------------------------------

def run_injection_power_curve(raw: pd.Series, exclude_centers: list[pd.Timestamp],
                               rng: np.random.Generator, channel_label: str) -> dict:
    idx = raw.index
    valid_lo = idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC + INJECTION_EXCLUSION_BUFFER_SEC)
    valid_hi = idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC + INJECTION_EXCLUSION_BUFFER_SEC)
    # 挑一個遠離真實事件和陣列邊緣的注入點。
    candidates = _draw_null_centers(rng, valid_lo, valid_hi, exclude_centers, INJECTION_EXCLUSION_BUFFER_SEC, 1)
    if not candidates:
        return {"channel": channel_label, "status": "no_valid_injection_point"}
    inject_at = candidates[0]
    inject_exclude = exclude_centers + [inject_at]

    # 從這個位置上未注入的虛無執行建立局部雜訊底（sigma）。
    baseline_d = _detrend(raw)
    baseline_null_centers = _draw_null_centers(rng, valid_lo, valid_hi, inject_exclude,
                                                INJECTION_EXCLUSION_BUFFER_SEC, N_NULL)
    m0 = STEP_WINDOWS_SEC[0]
    baseline_stat = _step_statistic(baseline_d, m0)
    baseline_null = [abs(r[0]) for c in baseline_null_centers
                      if (r := _extremum_in_window(baseline_stat, idx, c, SCAN_HALF_SEC)) is not None]
    sigma = float(np.median(baseline_null)) if baseline_null else 1.0
    if sigma <= 0:
        sigma = 1.0

    duration_sec = 2 * m0  # 短平台：先升後降，兩者都像階躍
    curve = []
    for k in INJECTION_AMPLITUDES_SIGMA:
        amplitude = k * sigma
        injected = raw.copy()
        mask = (idx >= inject_at) & (idx < inject_at + pd.Timedelta(seconds=duration_sec))
        injected.loc[mask] = injected.loc[mask] + amplitude
        d = _detrend(injected)
        stat = _step_statistic(d, m0)
        obs = _extremum_in_window(stat, idx, inject_at, SCAN_HALF_SEC)
        null_centers = _draw_null_centers(rng, valid_lo, valid_hi, inject_exclude,
                                           INJECTION_EXCLUSION_BUFFER_SEC, N_NULL)
        null_vals = [r for c in null_centers
                     if (r := _extremum_in_window(stat, idx, c, SCAN_HALF_SEC)) is not None]
        result = _summarize_test(obs, null_vals)
        result["injected_amplitude_nt"] = round(amplitude, 4)
        result["injected_amplitude_sigma_multiple"] = k
        curve.append(result)

    return {
        "channel": channel_label,
        "status": "ok",
        "inject_at_utc": inject_at.strftime("%Y-%m-%d %H:%M:%S"),
        "local_sigma_nt": round(sigma, 4),
        "step_window_sec": m0,
        "duration_sec": duration_sec,
        "curve": curve,
    }


# ---------------------------------------------------------------------------
# 合成資料自我測試（不用真實資料）：確認統計量 + 虛無／p 值
# 機制本身正確，透過 --self-test 執行。
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n = 3600
    noise = rng.normal(0, 1.0, size=n)
    step_at = 1800
    amplitude = 20.0
    injected = noise.copy()
    injected[step_at:] += amplitude

    ok = True
    for m in STEP_WINDOWS_SEC:
        stat = _step_statistic(injected, m)
        peak_idx = int(np.nanargmax(np.abs(stat)))
        peak_val = stat[peak_idx]
        near_step = abs(peak_idx - step_at) <= m + 2
        magnitude_ok = abs(peak_val) > amplitude * 0.5  # boxcar 平滑會壓低峰值；只是健全性界限
        status = "PASS" if (near_step and magnitude_ok) else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"[self-test] M={m:3d}s  peak_idx={peak_idx} (expected ~{step_at})  "
              f"peak_val={peak_val:.2f} (injected {amplitude})  {status}")

    spike = _spike_statistic(injected)
    spike_peak = int(np.nanargmax(spike))
    spike_status = "PASS" if abs(spike_peak - step_at) <= 2 else "FAIL"
    if spike_status == "FAIL":
        ok = False
    print(f"[self-test] spike peak_idx={spike_peak} (expected ~{step_at})  {spike_status}")
    return ok


# ---------------------------------------------------------------------------
# 逐事件統籌：對每一個非空的
# （F/XYZ）測站池、它的所有測站、所有通道檢定一個 Event。由
# 單一事件的 `run()` 入口和 117 起事件的 `run_all()` 掃描共用。
# ---------------------------------------------------------------------------

def _rank_stations_for_event(cfg: "common.GroupConfig", event, channel_type: str, n: int) -> list[tuple[str, float]]:
    """每個事件重新做 haversine 排名——cfg.stations[*]['distance_km']
    只相對於該組的錨點事件有效（見 common.py::
    _discover_stations），對多事件組中的其他任何事件都是錯的。"""
    wanted = "F" if channel_type == "F" else "XYZF"
    ranked = sorted(
        ((code, common.haversine_km(m["lat"], m["lon"], event.lat, event.lon))
         for code, m in cfg.stations.items() if m["reported"] == wanted),
        key=lambda x: x[1],
    )
    return ranked[:n]


def _build_channels(df: pd.DataFrame) -> dict[str, pd.Series]:
    if "F" in df.columns:
        return {"F": df["F"]}
    return {"X": df["X"], "Y": df["Y"], "Z": df["Z"], "H": np.sqrt(df["X"] ** 2 + df["Y"] ** 2)}


def _iter_test_payloads(result: dict):
    for m, payload in result.get("step", {}).items():
        yield "step", int(m), payload
    if "spike" in result:
        yield "spike", None, result["spike"]


def _storm_status(cfg: "common.GroupConfig", event_utc: pd.Timestamp) -> dict:
    """在該組日尺度的 storm_days.csv 中查詢事件的 UTC 日期
    （fetch_space_weather.py 的輸出，以 UTC 檔案日期為鍵）。只是旗標，給
    敏感度分析用——這裡什麼都不排除。窗口內的局部虛無分布
    已經吸收了被磁暴均勻抬高的雜訊底；它無法
    吸收的是落在 ±SCAN_HALF_SEC 搜尋窗口內的單一瞬變（SSC、亞暴開始），
    而這正是旗標讓你可以檢查的。
    檔案不存在或事件落在已抓取的日期範圍之外時，旗標為 None（未知），
    而不是默默讀成「平靜」。"""
    unknown = {"is_storm_day": None, "is_storm_onset": None, "storm_flag_confidence": None}
    csv_path = cfg.interim_dir / "storm_days.csv"
    summary_path = cfg.interim_dir / "storm_days_summary.json"
    if not csv_path.exists() or not summary_path.exists():
        return unknown
    summary = json.loads(summary_path.read_text())
    lo, hi = summary["date_range"]
    if not (pd.Timestamp(lo) <= event_utc.normalize() <= pd.Timestamp(hi)):
        return unknown
    storm_days = pd.read_csv(csv_path, dtype={"date": str})
    row = storm_days[storm_days.date == event_utc.strftime("%Y%m%d")]
    return {
        "is_storm_day": bool(row.is_storm_or_recovery.iloc[0]) if len(row) else False,
        "is_storm_onset": bool(row.is_storm_onset.iloc[0]) if len(row) else False,
        "storm_flag_confidence": summary.get("confidence"),
    }


def process_event(cfg: "common.GroupConfig", group, event,
                   run_injection: bool = False) -> tuple[dict, list[dict], list[dict]]:
    """對一個 Event 跑每一個 (channel_type, station, channel) 檢定。
    exclude_centers 涵蓋 `group` 中的每一個事件（不只這一個）——
    原因見模組 docstring 的「擴展到全部 117 起事件」說明。"""
    event_utc = pd.Timestamp(event.time_utc)
    # 用 folder_events 而不是 group.events：G2/G3 和 G6/G7/G8 共用原始資料夾，而
    # 同一份資料中另一組的真實事件也必須排除在這組的虛無抽樣之外。
    exclude_centers = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]
    scan_half_sec = _effective_half_sec(SCAN_HALF_SEC, event_utc, exclude_centers)
    storm = _storm_status(cfg, event_utc)

    event_result = {
        "group": cfg.group_id, "event_date": event.date, "event_time_local": event.time_local,
        "event_time_utc": event.time_utc, "magnitude": f"{event.magnitude_type}{event.magnitude}",
        "coord_confidence": event.coord_confidence, "anchor": event.anchor,
        **storm,
        "seed": SEED, "n_null": N_NULL, "step_windows_sec": list(STEP_WINDOWS_SEC),
        "scan_half_sec": scan_half_sec, "scan_half_sec_target": SCAN_HALF_SEC,
        "exclusion_buffer_sec": EXCLUSION_BUFFER_SEC,
        "results": [],
    }
    summary_rows: list[dict] = []
    injection_runs: list[dict] = []
    buffer_sec = REQUIRED_MARGIN_SEC + DETREND_WINDOW_SEC // 2 + 60  # 也替去趨勢的邊緣留充足餘裕

    for channel_type in ("F", "XYZ"):
        stations = _rank_stations_for_event(cfg, event, channel_type, N_STATIONS)
        for station, distance_km in stations:
            df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, buffer_sec)
            if df is None:
                event_result["results"].append({
                    "station": station, "channel_type": channel_type,
                    "distance_km": round(distance_km, 1),
                    "data_status": "no_data", "missing_dates": missing,
                })
                continue

            for ch_label, series in _build_channels(df).items():
                rng = keyed_rng(cfg.group_id, event.date, channel_type, station, ch_label)
                result = _test_channel(series, event_utc, exclude_centers, rng, ch_label, scan_half_sec)
                result["station"] = station
                result["channel_type"] = channel_type
                result["distance_km"] = round(distance_km, 1)
                result["missing_source_dates"] = missing
                event_result["results"].append(result)

                for stat_type, window, payload in _iter_test_payloads(result):
                    summary_rows.append({
                        "group": cfg.group_id, "event_date": event.date,
                        "magnitude": f"{event.magnitude_type}{event.magnitude}",
                        "coord_confidence": event.coord_confidence,
                        "is_storm_day": storm["is_storm_day"], "is_storm_onset": storm["is_storm_onset"],
                        "channel_type": channel_type, "station": station,
                        "distance_km": round(distance_km, 1), "channel": ch_label,
                        "statistic_type": stat_type, "window_sec": window,
                        "data_status": result.get("data_status"),
                        "obs_value_nt": payload.get("obs_value_nt"),
                        "obs_lag_sec": payload.get("obs_lag_sec"),
                        "p_value": payload.get("p_value"),
                        "n_null_samples": payload.get("n_null_samples"),
                        "approx_min_detectable_nt": payload.get("approx_min_detectable_nt"),
                    })

                if run_injection and result.get("data_status") == "ok":
                    injection_runs.append(
                        run_injection_power_curve(series, exclude_centers,
                                                  keyed_rng(cfg.group_id, event.date, channel_type, station, ch_label, "injection"),
                                                  f"{station}:{ch_label}")
                    )

    n_tests = len(summary_rows)
    event_result["multiple_comparisons_context"] = {
        "n_total_tests": n_tests,
        "bonferroni_alpha_at_p05": round(0.05 / n_tests, 6) if n_tests else None,
        "note": "No correction applied per-event -- see all_events_run_summary.json for the "
                "dataset-wide test count when running --all, and the module docstring.",
    }
    return event_result, summary_rows, injection_runs


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

def run(group_id: str = GROUP_ID, run_injection: bool = True) -> dict:
    """單一事件的冒煙測試／回歸檢查入口：只檢定
    `group_id` 的錨點事件（預設 G10），對照每一個非空的測站
    池。這就是 `coseismic_step_analysis.py`（不加旗標）跑的東西。"""
    cfg = common.load_group_config(group_id)
    group = get_group(group_id)
    event = group.anchor_event

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)

    event_result, summary_rows, injection_runs = process_event(cfg, group, event, run_injection=run_injection)

    out_path = OUT_DIR / "events" / f"{group_id}__{event.date}.json"
    out_path.write_text(json.dumps(event_result, indent=2))
    print(f"wrote {out_path} ({len(summary_rows)} tests)", file=sys.stderr)

    if run_injection:
        inj_path = OUT_DIR / "injection_power_curve.json"
        inj_path.write_text(json.dumps({"group": group_id, "event_date": event.date, "runs": injection_runs},
                                        indent=2))
        print(f"wrote {inj_path} ({len(injection_runs)} channel runs)", file=sys.stderr)

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "summary.csv", index=False)
    print(f"wrote summary.csv ({len(summary_rows)} rows)", file=sys.stderr)
    return event_result


def _storm_sensitivity(items: list[dict]) -> dict:
    """把 run_all 的逐事件彙總依磁暴旗標拆開，比較
    各子集中 p<0.05 的檢定比例。如果是磁暴驅動了
    偵測結果，磁暴子集的比例會明顯高於平靜
    子集；比例相近代表結果對磁暴不敏感。"""
    def _rollup(subset: list[dict]) -> dict:
        n_tests = sum(i["n_tests"] for i in subset)
        n_sig = sum(i["n_significant_p_lt_05"] for i in subset)
        return {
            "n_events": len(subset),
            "n_tests": n_tests,
            "n_significant_p_lt_05": n_sig,
            "frac_tests_p_lt_05": round(n_sig / n_tests, 4) if n_tests else None,
            "n_events_with_any_p_lt_05": sum(1 for i in subset if i["n_significant_p_lt_05"]),
        }

    out = {"all": _rollup(items)}
    for flag in ("is_storm_day", "is_storm_onset"):
        out[f"{flag}=True"] = _rollup([i for i in items if i[flag] is True])
        out[f"{flag}=False"] = _rollup([i for i in items if i[flag] is False])
    out["storm_flag_unknown"] = _rollup([i for i in items if i["is_storm_day"] is None])
    return out


def run_all(run_injection: bool = False) -> dict:
    """掃過 events.py 117 起事件登錄表中的每一個事件（全部 20 組、
    每組全部事件——不只 20 個錨點）。在 events/ 底下每個事件寫一個 JSON，
    寫一份跨所有事件的 summary.csv，以及
    all_events_run_summary.json（run_all_groups.sh 風格的逐項
    狀態彙總，所以不會有任何東西被默默跳過）。"""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "events").mkdir(exist_ok=True)

    all_summary_rows: list[dict] = []
    all_injection_runs: list[dict] = []
    run_summary = {"seed": SEED, "n_null": N_NULL, "events": []}

    for group_id, group in GROUPS.items():
        cfg = common.load_group_config(group_id)
        for event in group.events:
            event_result, summary_rows, injection_runs = process_event(
                cfg, group, event, run_injection=run_injection)

            out_path = OUT_DIR / "events" / f"{group_id}__{event.date}.json"
            out_path.write_text(json.dumps(event_result, indent=2))
            all_summary_rows.extend(summary_rows)
            all_injection_runs.extend(injection_runs)

            p_values = [r["p_value"] for r in summary_rows if r["p_value"] is not None]
            data_statuses = sorted({r["data_status"] for r in event_result["results"]})
            n_stations = len({(r["station"], r.get("channel_type")) for r in event_result["results"]
                               if "station" in r})
            item = {
                "group": group_id, "event_date": event.date, "anchor": event.anchor,
                "magnitude": f"{event.magnitude_type}{event.magnitude}",
                "coord_confidence": event.coord_confidence,
                "is_storm_day": event_result["is_storm_day"], "is_storm_onset": event_result["is_storm_onset"],
                "n_stations_tested": n_stations, "n_tests": len(summary_rows),
                "n_significant_p_lt_05": sum(1 for p in p_values if p < 0.05),
                "min_p": min(p_values) if p_values else None,
                "data_statuses_seen": data_statuses,
            }
            run_summary["events"].append(item)
            print(f"[{group_id} {event.date}] {item['n_tests']} tests, min_p={item['min_p']}, "
                  f"statuses={data_statuses}", file=sys.stderr)

    run_summary["n_events_total"] = len(run_summary["events"])
    run_summary["n_total_tests"] = len(all_summary_rows)
    run_summary["bonferroni_alpha_at_p05"] = (
        round(0.05 / len(all_summary_rows), 8) if all_summary_rows else None
    )
    run_summary["storm_sensitivity"] = _storm_sensitivity(run_summary["events"])

    (OUT_DIR / "all_events_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    pd.DataFrame(all_summary_rows).to_csv(OUT_DIR / "summary.csv", index=False)
    if run_injection:
        (OUT_DIR / "injection_power_curve.json").write_text(
            json.dumps({"runs": all_injection_runs}, indent=2))

    print(f"[run_all] {run_summary['n_events_total']} events, {run_summary['n_total_tests']} total tests "
          f"-> {OUT_DIR}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run only the synthetic sanity check")
    ap.add_argument("--no-injection", action="store_true", help="skip the positive-control injection test")
    ap.add_argument("--all", action="store_true",
                     help="run every event in events.py's 117-event registry instead of just G10's anchor")
    args = ap.parse_args()

    if args.self_test:
        passed = self_test()
        sys.exit(0 if passed else 1)

    if not self_test():
        print("[main] synthetic self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    if args.all:
        if not args.no_injection:
            print("[main] --all always runs with injection testing disabled (already validated once on "
                  "G10, see injection_power_curve.json / plan file) -- ignoring lack of --no-injection",
                  file=sys.stderr)
        run_all(run_injection=False)
    else:
        run(run_injection=not args.no_injection)
