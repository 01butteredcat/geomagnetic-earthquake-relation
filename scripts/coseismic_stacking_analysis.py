"""同震階躍／突波統計量的跨事件疊加時間
疊加——`coseismic_step_analysis.py` 自己的 docstring
（見它「擴展到全部 117 起事件」那一節）明確點名為
自然的下一步、但刻意沒有實作的工作：

    「跨事件疊加（用 stat_utils.mad_zscore 以每個事件自己的事件外
    統計量值標準化它的窗口，再像
    superposed_epoch_analysis.py 疊加日尺度序列那樣跨事件平均，
    並以每個事件的隨機參考時間建立虛無帶）
    是自然的下一步，如果 117 起事件的執行顯示
    許多個別檢定力不足、但有提示性的事件；但這一輪
    刻意不實作。」

`coseismic_step_analysis.py` 發現 31 起事件中有 19 起在發震秒附近
個別顯著，但它只報告每個事件在 ±120 秒搜尋窗口內的*極值*
統計量——從不報告完整的窗口形狀，所以
沒有辦法看出這些事件是否真的長得一樣（真實的共同機制
會產生一致的階躍／突波形狀），還是
只是各自以不相關的方式剛好跨過 p<0.05 門檻。這支
腳本的解答方式是：取出每個事件完整的窗口*剖面*（不
只是極值），用該事件自己的事件外
雜訊底做 z 標準化（讓絕對雜訊底差很多的事件——例如
安靜的純量 F 站 vs 雜訊大的近場 XYZ 站——在平均前
處在可比較的基礎上），再跨事件疊加——完全
仿照 `superposed_epoch_analysis.py::run_band` 日尺度的 bootstrap 信賴區間
+ 虛無帶設計，只是解析度從 1 天變成 1 秒。

## 為什麼是新檔案，而不是擴充 coseismic_step_analysis.py

那支腳本自己的輸出（`data/interim/coseismic_step_analysis/events/*.json`、
`summary.csv`、`all_events_run_summary.json`）已經被一份
已發布的報告引用。這支腳本直接匯入它的去趨勢／統計量／虛無抽樣
輔助函式（不是複製一份——不像那支腳本自己
刻意複製 `ulf_analysis.py::_detrend`，那是因為
兩支腳本的目的確實不同；這支則是
coseismic_step_analysis.py 本身同一作者的直接延續，重用
已在 117 起事件規模的真實資料上驗證過的機制），但寫到
自己新的輸出目錄樹，所以既有產物完全不會改變。

## 疊加什麼

四個統計量，分開保留（不合併）——step10/step30/step90/spike，
仿照 superposed_epoch_analysis.py 把 pc3/pc4 當成獨立
頻帶而不合併的做法。每個測站池的通道化約成一個標準、
可跨站比較的量：XYZ 測站池用 H=sqrt(X^2+Y^2)，
只有 F 的測站池用 F（X/Y 個別取決於測站
方位，不能直接跨站比較）。測站選擇
是逐事件的（不是該組相對錨點的 cfg.stations 距離，那對
非錨點事件是錯的——見 coseismic_step_analysis.py 自己的
_rank_stations_for_event 和 events.py 的 G14 note），分兩層：
對該特定事件最近和第二近的測站（"near"/"far"，
沿用 common.py 既有的 near/far 用語——不是宣稱 "far"
真的很遠，只是對這個事件排第二的測站）。兩層
都疊加，以便做距離敏感度比較：如果
第二近測站的疊加訊號弱很多，那支持
空間上局部（與近場震動耦合）的來源；如果差不多，那
不支持訊號是乾淨地局部化的。

疊加半窗口是 STACK_HALF_SEC=300（比 coseismic_step_
analysis.py 的 SCAN_HALF_SEC=120 寬），專門讓疊加形狀能顯示
平均訊號在偵測窗口之後很久是*回復*（G10 的「雜訊爆發後恢復」
型態）還是*持續升高*（G9 csg 的「持續、不恢復的偏移」
型態）——回復 vs 持續的
問題正是這支腳本的動機，而 120 秒沒留下多少
空間可以看到。

## 標準化（逐事件 x 測站 x 通道 x 統計量）

1. 在整個已載入緩衝區上去趨勢 + 計算統計量序列 `d`
   （原封不動地重用 coseismic_step_analysis.py 的 _detrend/_step_statistic/
   _spike_statistic）。
2. 事件外基準：`d` 的中位數／MAD（stat_utils.mad_zscore），限制在
   距離該事件所屬組*每一個*真實事件 ±EXCLUSION_BUFFER_SEC 以外的
   樣本（不只是被檢定的事件——有些組的
   事件在時間上夠近，例如 G9 的兩起事件相隔 17 小時，
   只排除單一事件可能讓另一起真實地震汙染
   「事件外」樣本池）。
3. 取出 [-STACK_HALF_SEC, +STACK_HALF_SEC] 窗口剖面，
   用那組（中位數, MAD）做 z 標準化——疊加的就是這條
   z 剖面。

## Bootstrap 信賴區間 + 虛無帶

直接把 superposed_epoch_analysis.py::run_band 移植到 1 秒尺度：
- Bootstrap CI90：對*事件*做 N_BOOTSTRAP=2000 次有放回重抽，
  每個延遲取重抽疊加平均的第 5／95 百分位數。
- 虛無帶：N_NULL_STACK=1000 次實現；每次為每個參與事件抽**一個**隨機
  參考時間（coseismic_step_analysis.py 自己的
  `_draw_null_centers`，遵守同樣的排除緩衝），在那裡重新取出
  該事件的窗口剖面，用*同一組*
  預先算好的（中位數, MAD）重新標準化（不重算——虛無帶問的是
  「在這個事件自己的雜訊底上隨機一秒看起來像不像這樣」，所以
  基準必須固定），再疊加。每個延遲在
  1000 次實現中取第 5/50/95 百分位數，得到用於視覺化的逐點虛無帶。
  真實疊加峰值的顯著性**不是**對照那條逐點
  虛無帶判斷（self_test() 顯示這會產生偽陽性——真實峰值
  是在 601 個延遲上取 argmax 選出的，所以拿它和單一延遲的逐點
  虛無帶比較，是有搜尋 vs 沒搜尋的不對等），而是對照
  每次虛無實現在所有延遲上**自己的**峰值 |z| 分布——和
  coseismic_step_
  analysis.py 自己的 p 值所用的「觀測和虛無用完全相同的搜尋程序」設計一樣。

固定 SEED=20260805（沿用這個程式庫其他地方的種子）以確保
可重現。

用法：
  coseismic_stacking_analysis.py --self-test          # 只做合成資料健全性檢查
  coseismic_stacking_analysis.py --group G9 --group G10   # 冒煙測試，2 組
  coseismic_stacking_analysis.py --all                # events.py 中的每一組
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
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
    MISSING_FRACTION_THRESHOLD,
    SEED,
    _build_channels,
    _detrend,
    _draw_null_centers,
    _load_station_days,
    keyed_rng,
    _pos,
    _rank_stations_for_event,
    _spike_statistic,
    _step_statistic,
    _storm_status,
)

STACK_HALF_SEC = 300
N_BOOTSTRAP = 2000
N_NULL_STACK = 1000
STATION_TIERS = ("near", "far")  # 依逐事件 haversine 距離排第 0 / 第 1——見模組 docstring
STAT_NAMES = ("step10", "step30", "step90", "spike")
REQUIRED_MARGIN_SEC = STACK_HALF_SEC + EXCLUSION_BUFFER_SEC  # 900 秒：兩端各需要的最少資料支撐
BUFFER_SEC = REQUIRED_MARGIN_SEC + DETREND_WINDOW_SEC // 2 + 60  # 2760 秒——仿照 coseismic_step_analysis.py 自己的 buffer_sec 公式
MIN_OFF_EVENT_SAMPLES = 200  # 低於這個，事件外中位數／MAD 的估計太吵，不能拿來標準化

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_stacking_analysis"


# ---------------------------------------------------------------------------
# 統計量／窗口輔助函式
# ---------------------------------------------------------------------------

def _stat_array(d: np.ndarray, stat_name: str) -> np.ndarray:
    if stat_name == "spike":
        return _spike_statistic(d)
    return _step_statistic(d, int(stat_name.replace("step", "")))


def _window_profile(arr: np.ndarray, idx: pd.DatetimeIndex, center: pd.Timestamp, half_sec: int) -> np.ndarray:
    """以 `center` 為中心的 `arr` 位置切片（和 coseismic_step_analysis.py::
    _extremum_in_window 相同的 O(window) 位置索引做法，
    只是保留整條剖面而不只是它的
    極值）。如果要求的窗口超出陣列
    邊界則補 NaN——依 BUFFER_SEC 的餘裕應該不會發生，但這裡做
    防禦性處理而不是假設。"""
    p = _pos(idx, center)
    n = len(arr)
    out = np.full(2 * half_sec + 1, np.nan)
    lo, hi = p - half_sec, p + half_sec + 1
    src_lo, src_hi = max(0, lo), min(n, hi)
    if src_hi <= src_lo:
        return out
    dst_lo = src_lo - lo
    out[dst_lo:dst_lo + (src_hi - src_lo)] = arr[src_lo:src_hi]
    return out


def _off_event_baseline(arr: np.ndarray, idx: pd.DatetimeIndex,
                         exclude_centers: list[pd.Timestamp], buffer_sec: int) -> tuple[float, float] | None:
    """`arr` 的中位數／MAD，排除距離 exclude_centers 中任何
    真實事件 buffer_sec 以內的每個樣本——也就是每個事件的
    事件窗口用來標準化的事件外雜訊底。如果排除後剩下的
    乾淨樣本太少、不足以信任估計，則為 None。"""
    mask = np.ones(len(arr), dtype=bool)
    for c in exclude_centers:
        p = _pos(idx, c)
        lo, hi = max(0, p - buffer_sec), min(len(arr), p + buffer_sec + 1)
        mask[lo:hi] = False
    vals = arr[mask]
    vals = vals[~np.isnan(vals)]
    if len(vals) < MIN_OFF_EVENT_SAMPLES:
        return None
    _, med, mad = mad_zscore(vals)
    if not np.isfinite(mad) or mad <= 1e-9:
        return None
    return med, mad


# ---------------------------------------------------------------------------
# EventSeries：把疊加核心和真實資料載入分開，讓
# self_test() 能在合成資料上執行完全相同的 stack_series() 程式碼，
# 不碰任何 .sec 檔。
# ---------------------------------------------------------------------------

@dataclass
class EventSeries:
    group_id: str
    event_date: str
    anchor: bool
    magnitude: str
    station: str
    distance_km: float
    idx: pd.DatetimeIndex
    stat_arrays: dict           # stat_name -> 完整的去趨勢統計量序列（np.ndarray）
    baselines: dict             # stat_name -> 事件外的 (median, mad)
    event_utc: pd.Timestamp
    valid_lo: pd.Timestamp
    valid_hi: pd.Timestamp
    exclude_centers: list       # list[pd.Timestamp]，這個 EventSeries 所屬組的每一個真實事件
    is_storm_day: bool | None = None    # storm_days.csv 的旗標（見 coseismic_step_analysis._storm_status）；
    is_storm_onset: bool | None = None  # None = 未知，只有磁暴敏感度子集會用到


def _build_event_series(cfg: "common.GroupConfig", group, event, station: str, distance_km: float,
                         raw: pd.Series) -> EventSeries | None:
    idx = raw.index
    event_utc = pd.Timestamp(event.time_utc)
    exclude_centers = [pd.Timestamp(e.time_utc) for e in folder_events(group.group_id)]

    if idx.min() + pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) > event_utc or \
       idx.max() - pd.Timedelta(seconds=REQUIRED_MARGIN_SEC) < event_utc:
        return None
    nan_frac = float(
        raw[(idx >= event_utc - pd.Timedelta(seconds=STACK_HALF_SEC)) &
            (idx <= event_utc + pd.Timedelta(seconds=STACK_HALF_SEC))].isna().mean()
    )
    if nan_frac > MISSING_FRACTION_THRESHOLD:
        return None

    storm = _storm_status(cfg, event_utc)
    d = _detrend(raw)
    stat_arrays, baselines = {}, {}
    for stat_name in STAT_NAMES:
        arr = _stat_array(d, stat_name)
        base = _off_event_baseline(arr, idx, exclude_centers, EXCLUSION_BUFFER_SEC)
        if base is None:
            continue
        stat_arrays[stat_name] = arr
        baselines[stat_name] = base
    if not stat_arrays:
        return None

    return EventSeries(
        group_id=cfg.group_id, event_date=event.date, anchor=event.anchor,
        magnitude=f"{event.magnitude_type}{event.magnitude}",
        station=station, distance_km=round(distance_km, 1),
        idx=idx, stat_arrays=stat_arrays, baselines=baselines,
        event_utc=event_utc,
        valid_lo=idx.min() + pd.Timedelta(seconds=STACK_HALF_SEC),
        valid_hi=idx.max() - pd.Timedelta(seconds=STACK_HALF_SEC),
        exclude_centers=exclude_centers,
        is_storm_day=storm["is_storm_day"], is_storm_onset=storm["is_storm_onset"],
    )


def load_all_event_series(group_ids: tuple[str, ...]) -> dict[tuple[str, str, str], EventSeries]:
    """回傳 {(channel_type_pool, tier, "<group_id>__<event.date>"): EventSeries}。
    每個所要求組別的每個事件，都會在兩種通道類型
    測站池（"F"/"XYZ"，看該組實際有哪一種測站——和
    coseismic_step_analysis.py::process_event 相同的雙測站池迴圈）和兩個
    測站層級上嘗試。全部 20 組的 117 起事件都符合資格——這
    只需要既有的 .sec 地磁資料，不依賴地震儀。"""
    out: dict[tuple[str, str, str], EventSeries] = {}
    for group_id in group_ids:
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)
        for event in group.events:
            event_utc = pd.Timestamp(event.time_utc)
            for channel_type in ("F", "XYZ"):
                stations = _rank_stations_for_event(cfg, event, channel_type, len(STATION_TIERS))
                for rank, tier in enumerate(STATION_TIERS):
                    if rank >= len(stations):
                        continue
                    station, distance_km = stations[rank]
                    df, missing, wanted = _load_station_days(cfg.gdms_dir, station, event_utc, BUFFER_SEC)
                    if df is None:
                        continue
                    channels = _build_channels(df)
                    ch_label = "H" if channel_type == "XYZ" else "F"
                    if ch_label not in channels:
                        continue
                    es = _build_event_series(cfg, group, event, station, distance_km, channels[ch_label])
                    if es is None:
                        continue
                    out[(channel_type, tier, f"{group_id}__{event.date}")] = es
    return out


# ---------------------------------------------------------------------------
# 疊加核心（完全在 EventSeries 物件上運作——真實資料
# 透過 run_group_ids()、合成資料透過 self_test() 都會執行到）。
# ---------------------------------------------------------------------------

def stack_series(events: list[EventSeries], stat_name: str, rng: np.random.Generator) -> dict:
    lags = np.arange(-STACK_HALF_SEC, STACK_HALF_SEC + 1)
    matrix: list[np.ndarray] = []
    usable: list[EventSeries] = []
    used: list[dict] = []
    excluded: list[dict] = []

    for es in events:
        if stat_name not in es.stat_arrays:
            excluded.append({"group": es.group_id, "date": es.event_date, "reason": "no_baseline_for_stat"})
            continue
        med, mad = es.baselines[stat_name]
        profile = _window_profile(es.stat_arrays[stat_name], es.idx, es.event_utc, STACK_HALF_SEC)
        if np.all(np.isnan(profile)):
            excluded.append({"group": es.group_id, "date": es.event_date, "reason": "empty_window"})
            continue
        matrix.append((profile - med) / mad)
        usable.append(es)
        used.append({"group": es.group_id, "date": es.event_date, "anchor": es.anchor,
                      "magnitude": es.magnitude, "station": es.station, "distance_km": es.distance_km})

    if not matrix:
        return {"stat": stat_name, "error": "no events with usable data", "events_excluded": excluded}

    M = np.array(matrix)
    n_events = M.shape[0]
    stack_mean = np.nanmean(M, axis=0)
    stack_median = np.nanmedian(M, axis=0)
    n_contributing = np.sum(~np.isnan(M), axis=0)

    # 真實疊加的 bootstrap 信賴區間（對事件做有放回重抽）
    boot = np.empty((N_BOOTSTRAP, M.shape[1]))
    for b in range(N_BOOTSTRAP):
        sample = rng.integers(0, n_events, size=n_events)
        boot[b] = np.nanmean(M[sample], axis=0)
    ci_lo = np.nanpercentile(boot, 5, axis=0)
    ci_hi = np.nanpercentile(boot, 95, axis=0)

    # 虛無帶：用每個事件的隨機參考時間取代真實發震時間，
    # 重複整個疊加程序。
    null_stacks = np.empty((N_NULL_STACK, M.shape[1]))
    for r in range(N_NULL_STACK):
        null_matrix = []
        for es in usable:
            centers = _draw_null_centers(rng, es.valid_lo, es.valid_hi, es.exclude_centers,
                                          EXCLUSION_BUFFER_SEC, 1)
            if not centers:
                continue
            med, mad = es.baselines[stat_name]
            profile = _window_profile(es.stat_arrays[stat_name], es.idx, centers[0], STACK_HALF_SEC)
            if np.all(np.isnan(profile)):
                continue
            null_matrix.append((profile - med) / mad)
        null_stacks[r] = np.nanmean(np.array(null_matrix), axis=0) if null_matrix else np.nan

    null_p5 = np.nanpercentile(null_stacks, 5, axis=0)
    null_p50 = np.nanpercentile(null_stacks, 50, axis=0)
    null_p95 = np.nanpercentile(null_stacks, 95, axis=0)

    peak_i = int(np.nanargmax(np.abs(stack_mean)))
    peak_val = float(stack_mean[peak_i])

    # 峰值的顯著性：**不是**在 peak_i 對照
    # null_p5/p95 做逐點比較——真實峰值是在全部
    # 601 個延遲上取 argmax 選出的，所以拿它和單一延遲的逐點虛無帶比較，
    # 是蘋果比橘子的「有搜尋 vs 沒搜尋」比較，會灌大偽
    # 陽性（已由 self_test() 的純雜訊檢查確認：純雜訊
    # 疊加的峰值常常光是因為有 601 次機會
    # 就超出逐點虛無帶）。正確、對等的比較——仿照
    # coseismic_step_analysis.py 自己「觀測和虛無用相同搜尋程序」
    # 的設計——是對照每次**虛無**
    # 實現在所有延遲上*自己的*峰值 |z| 分布。
    null_peak_abs = np.nanmax(np.abs(null_stacks), axis=1)
    null_peak_abs = null_peak_abs[~np.isnan(null_peak_abs)]
    if len(null_peak_abs):
        p_value = (1 + int(np.sum(null_peak_abs >= abs(peak_val)))) / (len(null_peak_abs) + 1)
        null_peak_p95 = float(np.percentile(null_peak_abs, 95))
    else:
        p_value = None
        null_peak_p95 = None
    outside = bool(p_value is not None and p_value < 0.05)

    return {
        "stat": stat_name,
        "half_sec": STACK_HALF_SEC,
        "seed": SEED, "n_bootstrap": N_BOOTSTRAP, "n_null_stacks": N_NULL_STACK,
        "n_events_used": n_events,
        "events_used": used,
        "events_excluded": excluded,
        "lags_sec": lags.tolist(),
        "stack_mean": [None if np.isnan(v) else round(float(v), 4) for v in stack_mean],
        "stack_median": [None if np.isnan(v) else round(float(v), 4) for v in stack_median],
        "n_contributing_per_lag": n_contributing.tolist(),
        "bootstrap_ci90_lo": [round(float(v), 4) for v in ci_lo],
        "bootstrap_ci90_hi": [round(float(v), 4) for v in ci_hi],
        "null_band_p5": [round(float(v), 4) for v in null_p5],
        "null_band_p50": [round(float(v), 4) for v in null_p50],
        "null_band_p95": [round(float(v), 4) for v in null_p95],
        "peak_abs_z": round(abs(peak_val), 4),
        "peak_lag_sec": int(lags[peak_i]),
        "null_peak_abs_z_p95": None if null_peak_p95 is None else round(null_peak_p95, 4),
        "p_value": p_value,
        "outside_null_band_at_peak": outside,
    }


# ---------------------------------------------------------------------------
# 合成資料自我測試：直接執行 stack_series()（不用真實資料），
# 確認 (a) 注入已知階躍的事件疊加後，在延遲=0 附近得到
# 落在自己虛無帶之外的峰值，以及 (b) 只有純雜訊的事件
# 疊加後，**不會**在虛無帶之外出現假峰值。
# ---------------------------------------------------------------------------

def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n_synth = 20
    n_samples = 2 * BUFFER_SEC + 1
    idx = pd.date_range("2024-01-01", periods=n_samples, freq="s")
    amplitude = 15.0
    center_i = BUFFER_SEC

    def make_series(inject: bool) -> EventSeries:
        noise = rng.normal(0, 1.0, size=n_samples)
        if inject:
            jitter = int(rng.integers(-5, 6))
            noise[center_i + jitter:center_i + jitter + 20] += amplitude
        event_utc = idx[center_i]
        arr = _step_statistic(noise, 30)
        base = _off_event_baseline(arr, idx, [event_utc], EXCLUSION_BUFFER_SEC)
        assert base is not None, "self-test off-event baseline computation failed"
        return EventSeries(
            group_id="SYN", event_date=f"synthetic-{int(rng.integers(0, 10**6))}", anchor=True, magnitude="M0",
            station="SYN", distance_km=0.0, idx=idx, stat_arrays={"step30": arr}, baselines={"step30": base},
            event_utc=event_utc,
            valid_lo=idx[0] + pd.Timedelta(seconds=STACK_HALF_SEC),
            valid_hi=idx[-1] - pd.Timedelta(seconds=STACK_HALF_SEC),
            exclude_centers=[event_utc],
        )

    signal_events = [make_series(inject=True) for _ in range(n_synth)]
    noise_events = [make_series(inject=False) for _ in range(n_synth)]

    rng2 = np.random.default_rng(SEED + 1)
    signal_result = stack_series(signal_events, "step30", rng2)
    noise_result = stack_series(noise_events, "step30", rng2)

    ok = True
    near_zero = abs(signal_result["peak_lag_sec"]) <= 10
    outside = signal_result["outside_null_band_at_peak"]
    status1 = "PASS" if (near_zero and outside) else "FAIL"
    print(f"[self-test] injected-signal stack: peak_lag={signal_result['peak_lag_sec']}s "
          f"peak_abs_z={signal_result['peak_abs_z']} outside_null={outside}  {status1}")
    ok = ok and near_zero and outside

    noise_outside = noise_result["outside_null_band_at_peak"]
    status2 = "PASS" if not noise_outside else "FAIL"
    print(f"[self-test] noise-only stack (false-positive check): peak_abs_z={noise_result['peak_abs_z']} "
          f"outside_null={noise_outside}  {status2}")
    ok = ok and not noise_outside

    return ok


# ---------------------------------------------------------------------------
# 真實資料的統籌
# ---------------------------------------------------------------------------

# 磁暴敏感度子集：同樣的疊加，去掉被標記為磁暴的事件後重跑
# （另外也只去掉較嚴格的磁暴開始日）。每個（組合,
# 子集）都從自己的 keyed_rng 亂數流抽樣，所以兩輪互不干擾。
# 旗標未知的事件保留——去掉它們會把「沒有
# 太空天氣資料」和「磁暴」混為一談。
STORM_SUBSETS = {
    "exclude_storm_or_recovery": lambda es: es.is_storm_day is not True,
    "exclude_storm_onset": lambda es: es.is_storm_onset is not True,
}


def _storm_sensitivity_stacks(all_series: dict) -> list[dict]:
    rows: list[dict] = []
    for channel_type, ch_label in (("XYZ", "H"), ("F", "F")):
        for tier in STATION_TIERS:
            events = [es for (ct, t, _key), es in all_series.items() if ct == channel_type and t == tier]
            if not events:
                continue
            for subset_name, keep in STORM_SUBSETS.items():
                kept = [es for es in events if keep(es)]
                for stat_name in STAT_NAMES:
                    combo_id = f"{tier}__{ch_label}__{stat_name}"
                    row = {"combo_id": combo_id, "subset": subset_name,
                           "n_events_before": len(events), "n_events_kept": len(kept)}
                    rng = keyed_rng("stack", combo_id, subset_name)
                    result = stack_series(kept, stat_name, rng) if kept else {"error": "no events left"}
                    if "error" in result:
                        row["error"] = result["error"]
                    else:
                        row.update({k: result[k] for k in ("n_events_used", "peak_abs_z", "peak_lag_sec",
                                                            "p_value", "outside_null_band_at_peak")})
                    rows.append(row)
                    detail = row.get("error") or f"peak_abs_z={row['peak_abs_z']} p={row['p_value']}"
                    print(f"[storm:{subset_name}] {combo_id} kept={len(kept)}/{len(events)} {detail}",
                          file=sys.stderr)
    pd.DataFrame(rows).to_csv(OUT_DIR / "stack_storm_sensitivity.csv", index=False)
    return rows


def run_group_ids(group_ids: tuple[str, ...]) -> dict:
    all_series = load_all_event_series(group_ids)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "combos").mkdir(exist_ok=True)

    channel_labels = {"XYZ": "H", "F": "F"}
    summary_rows: list[dict] = []
    run_summary = {"seed": SEED, "group_ids": list(group_ids), "combos": []}

    for channel_type, ch_label in channel_labels.items():
        for tier in STATION_TIERS:
            events = [es for (ct, t, _key), es in all_series.items() if ct == channel_type and t == tier]
            if not events:
                continue
            for stat_name in STAT_NAMES:
                combo_id = f"{tier}__{ch_label}__{stat_name}"
                result = stack_series(events, stat_name, keyed_rng("stack", combo_id))
                result["combo_id"] = combo_id
                result["station_tier"] = tier
                result["channel_type_pool"] = channel_type
                result["channel"] = ch_label

                out_path = OUT_DIR / "combos" / f"stack__{combo_id}.json"
                out_path.write_text(json.dumps(result, indent=2))

                if "error" in result:
                    print(f"[{combo_id}] {result['error']}", file=sys.stderr)
                    run_summary["combos"].append({"combo_id": combo_id, "status": "error", "error": result["error"]})
                    continue

                print(f"[{combo_id}] n_events={result['n_events_used']} "
                      f"peak_abs_z={result['peak_abs_z']} lag={result['peak_lag_sec']}s "
                      f"outside_null={result['outside_null_band_at_peak']}", file=sys.stderr)
                run_summary["combos"].append({
                    "combo_id": combo_id, "status": "ok",
                    "n_events_used": result["n_events_used"],
                    "peak_abs_z": result["peak_abs_z"], "peak_lag_sec": result["peak_lag_sec"],
                    "outside_null_band_at_peak": result["outside_null_band_at_peak"],
                })
                summary_rows.append({
                    "combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                    "channel": ch_label, "statistic": stat_name,
                    "n_events_used": result["n_events_used"],
                    "peak_abs_z": result["peak_abs_z"], "peak_lag_sec": result["peak_lag_sec"],
                    "outside_null_band_at_peak": result["outside_null_band_at_peak"],
                })

    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "stack_summary.csv", index=False)
    run_summary["storm_sensitivity"] = _storm_sensitivity_stacks(all_series)
    (OUT_DIR / "all_stacks_run_summary.json").write_text(json.dumps(run_summary, indent=2))
    print(f"[run] {len(summary_rows)} combos -> {OUT_DIR}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="run only the synthetic sanity check")
    ap.add_argument("--group", action="append", dest="groups",
                     help="restrict to this group ID (repeatable); default is every group in events.py")
    ap.add_argument("--all", action="store_true", help="explicit alias for the default (no --group filter)")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] synthetic self-test FAILED -- aborting before touching real data", file=sys.stderr)
        sys.exit(1)

    group_ids = tuple(args.groups) if args.groups else tuple(GROUPS.keys())
    run_group_ids(group_ids)
