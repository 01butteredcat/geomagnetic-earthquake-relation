"""起始時間走時差檢定：同震地磁異常在離震源越遠的測站
是否越晚開始？

這個檢定既不需要地震儀資料，也不需要雜訊／訊號標籤。三種
候選來源預測的全網起始型態各不相同：

  - 外部擾動（SSC、亞暴開始、Pi2）：各處同時發生，
    起始時間和距離無關（斜率 ~ 0 s/km）；
  - 磁力儀被地動震到：起始時間跟著震波走，
    起始 ~ R / V，V ~ 6 km/s（P）或 ~ 3.5 km/s（S），也就是斜率約
    0.17-0.29 s/km。20 -> 200 km 的差距約 30-50 秒，在 1 Hz 下可以清楚分辨；
  - 震源處的岩石圈來源：震央附近最早，
    可能比 P 波還早，但沒有理由以地震波速傳播。

## 起始偵測器

在每個測站原始的 1 Hz 一階差分上（向量站用 H = sqrt(X^2+Y^2)，
純量站用 F——和 coseismic_step_analysis.py 的
`_build_channels` 相同的通道）：

  sigma  = 1.4826 * NOISE_PRE_LAGS 內一階差分的 MAD
           （-660..-60 秒，和
           seismometer_comparison.py::geomag_noise_ratio 相同的震前參考）
  exceed = |diff| > K_SIGMA * sigma
  onset  = SEARCH_LAGS（-60..+180 秒，每個事件再由
           coseismic_step_analysis._effective_half_sec 設上限）中第一個 exceed[t] 成立、
           且 t..t+4 這 MIN_EXCEED_WINDOW 秒中至少 MIN_EXCEED 秒 exceed 的秒 t。

這是 G10 試行時看到的震動雜訊特徵（一階差分雜訊
約是震前水準的 11.6 倍），計時取它的第一個樣本，而不是
step30 峰值，後者可能落在爆發段中的任何地方。

## 虛無校準

同一個偵測器在同一份已載入資料中的 N_NULL 個隨機參考時間上執行
（原始資料夾中每個真實事件都保持 NULL_EXCLUSION_SEC 的距離，讓
整個 -660..+180 秒窗口都乾淨）。2026-09-28 起，參考
時間由**同一事件的所有測站共用**（每個事件抽一次，只保留
發震時被檢定的每個測站在該時間也都得到有效結果的；
這樣的時間少於 MIN_SHARED_NULL 個時，該事件退回把
無效的抽樣算成「未觸發」，並標記 `shared_null_fallback`）。在
那之前每個測站各自抽時間，這掩蓋了外部
擾動會同時打到整個測站網的事實。觸發的比例就是該
測站在這個偵測器下的誤觸發率；`trigger_p` = (1 + #虛無觸發
至少一樣早的次數) / (N_NULL + 1)，下面沒有拿來做任何事，但
可以讓一次觸發對照它自己測站的雜訊來解讀。

## 走時差統計量

震源距離 R = sqrt(震央距^2 + 深度^2)。

  - 逐事件（>= MIN_STATIONS_FOR_FIT 個觸發測站）：起始延遲對 R 的
    Theil-Sen 斜率，附 90% 信賴區間，以及視速度 1/slope。
    信賴區間包含 0 且低於 P 波慢度時為 `simultaneous`；
    信賴區間高於 0 且和 [1/V_P_KMS, 1/V_S_KMS] 重疊時為 `seismic_moveout`；
    其他為 `indeterminate`。
  - 跨事件合併：扣掉每個事件自己的延遲中位數和 R 中位數後，
    （起始延遲 vs. R）的 Theil-Sen 斜率（讓發震偏移
    不同的事件不會假裝成走時差），p 值來自在*同一事件*的
    測站間打亂起始時間（單尾，H1：slope > 0）。M >= 6 和 M < 6
    （補登的 68 起 M5 事件）分開報告，也分成
    全部觸發 vs. 只看自己誤觸發率低於
    CLEAN_FALSE_RATE 的測站。
  - 到時窗檢定（主要檢定，2026-09-28 起用事件層級虛無分布）：T =
    起始時間落在 [R/V_P - 5 s, R/V_S + 30 s] 內的測站數。
    虛無分布在一個共用的隨機參考時間上評估整個事件：每次
    模擬每個事件抽一個共用索引，計算該事件
    虛無觸發落在自己窗口內的測站數，再對事件加總。
    這保留了同一事件測站之間的相關性（一次亞暴
    會讓它們同時觸發），而較舊的逐觸發 Poisson-二項檢定
    （觀測計數 vs. 每次觸發虛無落窗比例之和，保留
    為 `arrival_window_all_triggers` 供參考）把它們當成獨立。
    這是「和震波綁在一起」的直接檢定；斜率擬合
    很容易被單一個雜訊大的測站拉走。
  - 和 seismometer_comparison.py 的同址地震起始時間交叉檢查
    （`shaking_onset_lag_sec`）：地磁起始減地震起始，針對
    那支腳本每個事件比對的那一個測站。

## 前一事件的震動（只做敏感度分析）

在同一原始資料夾中、緊接在另一個已登錄事件之後
<= PRIOR_EVENT_EXCLUSION_SEC（900 秒：-660 秒的雜訊參考加上
前一事件幾分鐘的震動）發生的事件，它的雜訊參考和搜尋
窗口都在前一次震動之中，所以它的起始時間不是對照
平靜基準計時的。這類事件會得到 `prior_event_shaking` = True。主要子集
（`m6`、`m5`）保留它們——主要檢定是在注意到這件事之前定下的——
`m6_no_prior_shaking` / `m5_no_prior_shaking` 則在排除它們後重算每一個統計量，
各自用自己的逐鍵 rng 亂數流，所以主要數字不會
變動。只能檢查已登錄的事件：M5 事件只從
2024-09 起有登錄，GDMS json 匯出檔只有 M>=6，所以較小的餘震
（例如 G10 2024-04-03c 內、M7.2 後 2 小時的）無法排除。

用法：
  coseismic_onset_moveout.py --self-test
  coseismic_onset_moveout.py --all [--min-mag 6] [--group G10 ...]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, theilslopes

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from coseismic_step_analysis import (SEED, _build_channels, _effective_half_sec,  # noqa: E402
                                     _load_station_days, keyed_rng)
from events import GROUPS, folder_events  # noqa: E402
from stat_utils import MAD_SCALE  # noqa: E402

NOISE_PRE_LAGS = (-660, -60)
SEARCH_LAGS = (-60, 180)
K_SIGMA = 4.0
MIN_EXCEED_WINDOW = 5
MIN_EXCEED = 3
MIN_PRE_SAMPLES = 300
MAX_SEARCH_MISSING = 0.2
LOAD_BUFFER_SEC = 6 * 3600
N_NULL = 200
MIN_SHARED_NULL = 50  # 全部有效的共用虛無時間少於這個數就退回（見 docstring）
NULL_EXCLUSION_SEC = -NOISE_PRE_LAGS[0] + 600
V_P_KMS, V_S_KMS = 6.0, 3.5
MIN_STATIONS_FOR_FIT = 3
CI_ALPHA = 0.90
CLEAN_FALSE_RATE = 0.05
CROSSCHECK_FALSE_RATE = 0.2
PRIOR_EVENT_EXCLUSION_SEC = 900  # -660 秒的雜訊參考 + 前一事件的震動
ARRIVAL_PAD_SEC = (-5, 30)  # 「到時窗」= 發震後 [R/V_P - 5, R/V_S + 30] 秒
N_SIM = 20000
N_PERM = 2000
N_WORKERS = 8

COMPARISON_CSV = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison" / "comparison_summary.csv"
OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_onset_moveout"


# ---------------------------------------------------------------------------
# 偵測器
# ---------------------------------------------------------------------------

def detect_onset(values: np.ndarray, p0: int, search_hi: int) -> dict:
    """values：無缺口的 1 Hz 網格（NaN = 缺值）；p0：參考時間在陣列中的
    位置；搜尋窗口是 SEARCH_LAGS[0]..min(SEARCH_LAGS[1], search_hi)。"""
    lo_pre, hi_pre = p0 + NOISE_PRE_LAGS[0], p0 + NOISE_PRE_LAGS[1]
    lo_s, hi_s = p0 + SEARCH_LAGS[0], p0 + min(SEARCH_LAGS[1], search_hi)
    if lo_pre < 1 or hi_s + MIN_EXCEED_WINDOW >= len(values):
        return {"status": "insufficient_window"}
    diff = values[lo_pre:hi_s + MIN_EXCEED_WINDOW + 1] - values[lo_pre - 1:hi_s + MIN_EXCEED_WINDOW]
    pre = diff[:hi_pre - lo_pre]
    pre = pre[~np.isnan(pre)]
    if len(pre) < MIN_PRE_SAMPLES:
        return {"status": "insufficient_pre"}
    sigma = float(np.median(np.abs(pre - np.median(pre)))) * MAD_SCALE
    if sigma <= 1e-9:
        return {"status": "flat_pre"}
    tail = diff[lo_s - lo_pre:]
    search = tail[:hi_s - lo_s + 1]
    if np.mean(np.isnan(search)) > MAX_SEARCH_MISSING:
        return {"status": "gappy_search", "sigma_nt": sigma}
    exceed = np.nan_to_num(np.abs(tail), nan=0.0) > K_SIGMA * sigma
    counts = np.convolve(exceed.astype(int), np.ones(MIN_EXCEED_WINDOW, dtype=int), "full")[
        MIN_EXCEED_WINDOW - 1:]  # counts[i] = exceed[i..i+4] 的個數
    hits = np.flatnonzero(exceed[:len(search)] & (counts[:len(search)] >= MIN_EXCEED))
    peak = float(np.nanmax(np.abs(search))) / sigma if not np.all(np.isnan(search)) else None
    if len(hits) == 0:
        return {"status": "no_onset", "sigma_nt": sigma, "peak_ratio": peak}
    i = int(hits[0])
    burst = np.abs(search[i:i + 60])
    return {"status": "onset", "sigma_nt": sigma, "peak_ratio": peak,
            "onset_lag_sec": int(i + SEARCH_LAGS[0]),
            "burst_peak_ratio": float(np.nanmax(burst)) / sigma}


def arrival_window(hypocentral_km: float) -> tuple[float, float]:
    return hypocentral_km / V_P_KMS + ARRIVAL_PAD_SEC[0], hypocentral_km / V_S_KMS + ARRIVAL_PAD_SEC[1]


def prior_event_gap(event_utc: pd.Timestamp, others: list[pd.Timestamp]) -> float | None:
    """距離最近一個較早的已登錄事件的秒數（沒有則為 None）。"""
    gaps = [(event_utc - t).total_seconds() for t in others]
    gaps = [g for g in gaps if g > 0]
    return min(gaps) if gaps else None


def _grid(series: pd.Series) -> pd.Series:
    return series.reindex(pd.date_range(series.index[0], series.index[-1], freq="s"))


def shared_null_outcomes(stations: dict[str, tuple[np.ndarray, pd.Timestamp]], event_utc: pd.Timestamp,
                         exclude: list[pd.Timestamp], rng: np.random.Generator,
                         search_hi: int) -> tuple[dict[str, list[int | None]], bool]:
    """每個測站在同樣 N_NULL 個隨機參考時間的偵測結果
    （None = 未觸發）。stations：名稱 -> (網格化的值, 第一個
    時間戳)。只保留每個測站都得到有效結果的時間；
    這樣的時間少於 MIN_SHARED_NULL 個時，改用前 N_NULL 個候選時間，
    無效結果算成未觸發（fallback = True）。"""
    lo = -LOAD_BUFFER_SEC - NOISE_PRE_LAGS[0] + 1
    hi = LOAD_BUFFER_SEC - SEARCH_LAGS[1] - MIN_EXCEED_WINDOW - 2
    ex = np.array([(e - event_utc).total_seconds() for e in exclude])
    kept: dict[str, list] = {st: [] for st in stations}
    loose: dict[str, list] = {st: [] for st in stations}
    n_kept = n_loose = attempts = 0
    while n_kept < N_NULL and attempts < N_NULL * 5:
        attempts += 1
        off = int(rng.integers(lo, hi))
        if len(ex) and np.min(np.abs(ex - off)) < NULL_EXCLUSION_SEC:
            continue
        t = event_utc + pd.Timedelta(seconds=off)
        res = {}
        for st, (values, t0) in stations.items():
            r = detect_onset(values, int(round((t - t0).total_seconds())), search_hi)
            res[st] = r["onset_lag_sec"] if r["status"] == "onset" else (None if r["status"] == "no_onset" else "invalid")
        if n_loose < N_NULL:
            n_loose += 1
            for st, v in res.items():
                loose[st].append(None if v == "invalid" else v)
        if all(v != "invalid" for v in res.values()):
            n_kept += 1
            for st, v in res.items():
                kept[st].append(v)
    return (kept, False) if n_kept >= MIN_SHARED_NULL else (loose, True)


# ---------------------------------------------------------------------------
# 逐事件／逐組
# ---------------------------------------------------------------------------

def process_group(group_id: str, min_mag: float | None = None) -> tuple[list[dict], dict]:
    """資料列（每個事件 x 測站一列），加上每個事件中每個被檢定測站
    在共用參考時間上的虛無落窗指標（給
    事件層級到時窗檢定用）。"""
    cfg = common.load_group_config(group_id)
    group = GROUPS[group_id]
    exclude = [pd.Timestamp(e.time_utc) for e in folder_events(group_id)]
    rows: list[dict] = []
    nulls: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for event in group.events:
        if min_mag is not None and event.magnitude < min_mag:
            continue
        event_utc = pd.Timestamp(event.time_utc)
        search_hi = _effective_half_sec(SEARCH_LAGS[1], event_utc, exclude)
        prior = prior_event_gap(event_utc, exclude)
        ev_rows: list[dict] = []
        loaded: dict[str, tuple[np.ndarray, pd.Timestamp]] = {}
        for station, meta in sorted(cfg.stations.items()):
            if meta["reported"] not in ("F", "XYZF"):
                continue
            epi = common.haversine_km(meta["lat"], meta["lon"], event.lat, event.lon)
            base = {"group": group_id, "event_date": event.date, "event_time_utc": event.time_utc,
                    "magnitude": event.magnitude, "anchor": event.anchor, "depth_km": event.depth_km,
                    "prior_event_sec": prior,
                    "prior_event_shaking": prior is not None and prior <= PRIOR_EVENT_EXCLUSION_SEC,
                    "station": station, "channel": "F" if meta["reported"] == "F" else "H",
                    "epicentral_km": round(epi, 1),
                    "hypocentral_km": round(float(np.hypot(epi, event.depth_km)), 1)}
            df, _, _ = _load_station_days(cfg.gdms_dir, station, event_utc, LOAD_BUFFER_SEC)
            if df is None:
                ev_rows.append({**base, "status": "no_data"})
                continue
            channels = _build_channels(df)
            if base["channel"] not in channels:  # 檔頭說是 F 但日檔帶的是 XYZ，或反過來
                base["channel"] = next(iter(channels)) if len(channels) == 1 else "H"
            series = _grid(channels[base["channel"]])
            if not (series.index[0] <= event_utc <= series.index[-1]):
                ev_rows.append({**base, "status": "no_data"})
                continue
            values = series.to_numpy(dtype=float)
            p0 = int(round((event_utc - series.index[0]).total_seconds()))
            res = detect_onset(values, p0, search_hi)
            if res["status"] in ("onset", "no_onset"):
                loaded[station] = (values, series.index[0])
            ev_rows.append({**base, **res})

        outcomes, fallback = shared_null_outcomes(loaded, event_utc, exclude,
                                                  keyed_rng("onset", group_id, event.date, "shared"), search_hi) \
            if loaded else ({}, False)
        ev_null: dict[str, np.ndarray] = {}
        for row in ev_rows:
            null = outcomes.get(row["station"])
            if null is None:
                continue
            n_trig = sum(v is not None for v in null)
            lo_a, hi_a = arrival_window(row["hypocentral_km"])
            row.update({"n_null": len(null), "shared_null_fallback": fallback,
                        "null_false_rate": round(n_trig / len(null), 4) if null else None})
            ev_null[row["station"]] = np.array([v is not None and lo_a <= v <= hi_a for v in null])
            if row["status"] == "onset":
                null_on = [v for v in null if v is not None]
                row["in_arrival_window"] = bool(lo_a <= row["onset_lag_sec"] <= hi_a)
                # 以這個測站自己的雜訊觸發時，落在同一窗口的機率；
                # 如果虛無分布從未觸發，就用搜尋窗口上的均勻分布
                row["null_q_arrival"] = round(float(np.mean([lo_a <= v <= hi_a for v in null_on])) if null_on else
                                              max(0.0, min(hi_a, search_hi) - max(lo_a, SEARCH_LAGS[0]) + 1)
                                              / (search_hi - SEARCH_LAGS[0] + 1), 4)
                if null:
                    row["trigger_p"] = round((1 + sum(v is not None and v <= row["onset_lag_sec"] for v in null))
                                             / (len(null) + 1), 4)
        for row in ev_rows:
            for k in ("sigma_nt", "peak_ratio", "burst_peak_ratio"):
                if row.get(k) is not None:
                    row[k] = round(row[k], 4)
        rows.extend(ev_rows)
        nulls[(group_id, event.date)] = ev_null
        print(f"[{group_id}] {event.date} M{event.magnitude}: "
              f"{sum(r.get('status') == 'onset' for r in ev_rows)} onsets"
              f"{'（共用虛無分布退回）' if fallback else ''}", file=sys.stderr)
    return rows, nulls


# ---------------------------------------------------------------------------
# 走時差統計量
# ---------------------------------------------------------------------------

def classify_slope(lo: float, hi: float) -> str:
    if lo <= 0 <= hi and hi < 1 / V_P_KMS:
        return "simultaneous"
    if lo > 0 and lo <= 1 / V_S_KMS and hi >= 1 / V_P_KMS:
        return "seismic_moveout"
    return "indeterminate"


def per_event_moveout(onsets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (g, d), ev in onsets.groupby(["group", "event_date"], sort=False):
        trig = ev[ev.status == "onset"]
        row = {"group": g, "event_date": d, "magnitude": ev.magnitude.iloc[0],
               "prior_event_shaking": bool(ev.prior_event_shaking.iloc[0]),
               "n_stations_tested": int(ev.status.isin(["onset", "no_onset"]).sum()),
               "n_triggered": int(len(trig)),
               "nearest_triggered_km": None if trig.empty else float(trig.hypocentral_km.min()),
               "earliest_onset_sec": None if trig.empty else int(trig.onset_lag_sec.min()),
               "median_resid_p_sec": None if trig.empty else
               float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_P_KMS)),
               "median_resid_s_sec": None if trig.empty else
               float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_S_KMS))}
        if len(trig) >= MIN_STATIONS_FOR_FIT and trig.hypocentral_km.nunique() > 1:
            s, _, lo, hi = theilslopes(trig.onset_lag_sec, trig.hypocentral_km, alpha=CI_ALPHA)
            row.update({"slope_s_per_km": round(s, 4), "slope_ci_lo": round(lo, 4), "slope_ci_hi": round(hi, 4),
                        "apparent_velocity_kms": round(1 / s, 2) if s > 0 else None,
                        "verdict": classify_slope(lo, hi)})
        else:
            row["verdict"] = "too_few_stations"
        rows.append(row)
    return pd.DataFrame(rows)


def pooled_moveout(trig: pd.DataFrame, rng: np.random.Generator, n_perm: int = N_PERM) -> dict:
    """在事件內置中的 (R, onset) 上做 Theil-Sen；p 值來自在同一事件的
    測站間打亂起始時間，單尾 H1：slope > 0。"""
    trig = trig.groupby(["group", "event_date"]).filter(lambda e: len(e) >= 2)
    if len(trig) < 6:
        return {"n_onsets": int(len(trig)), "note": "too few multi-station events"}
    key = trig.group + "|" + trig.event_date
    x = (trig.hypocentral_km - trig.groupby(key).hypocentral_km.transform("median")).to_numpy()
    y = trig.onset_lag_sec.to_numpy(dtype=float)
    ym = trig.groupby(key).onset_lag_sec.transform("median").to_numpy(dtype=float)
    s, _, lo, hi = theilslopes(y - ym, x, alpha=CI_ALPHA)
    idx = [np.flatnonzero((key == k).to_numpy()) for k in key.unique()]
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = y.copy()
        for ix in idx:
            yp[ix] = y[rng.permutation(ix)]
        null[i] = theilslopes(yp - ym, x)[0]
    p = float((1 + np.sum(null >= s)) / (n_perm + 1))
    return {"n_onsets": int(len(trig)), "n_events": len(idx),
            "slope_s_per_km": round(float(s), 4), "slope_ci90": [round(float(lo), 4), round(float(hi), 4)],
            "apparent_velocity_kms": round(1 / s, 2) if s > 0 else None,
            "p_one_sided_within_event": round(p, 5), "verdict": classify_slope(lo, hi)}


def arrival_window_test(trig: pd.DataFrame, rng: np.random.Generator) -> dict:
    """起始時間落在地震到時窗的頻率，是否比同一批
    測站自己的雜訊觸發更高？觀測計數 vs. 以逐觸發機率
    null_q_arrival 建立的 Poisson-二項虛無分布（單尾）。"""
    if trig.empty:
        return {"n_onsets": 0}
    q = trig.null_q_arrival.to_numpy(dtype=float)
    obs = int(trig.in_arrival_window.astype(bool).sum())
    sims = (rng.random((N_SIM, len(q))) < q).sum(axis=1)
    return {"n_onsets": int(len(trig)), "n_in_window": obs, "expected_by_chance": round(float(q.sum()), 2),
            "p_one_sided": round(float((1 + np.sum(sims >= obs)) / (N_SIM + 1)), 5),
            "median_resid_s_sec": round(float(np.median(trig.onset_lag_sec - trig.hypocentral_km / V_S_KMS)), 2)}


def arrival_window_event_level(tested: pd.DataFrame, nulls: dict, rng: np.random.Generator) -> dict:
    """主要的到時窗檢定。tested：status 為 onset/no_onset 的資料列。
    T = 落在自己到時窗內的起始時間數。虛無分布：每次模擬每個事件
    一個共用參考時間；計算該事件虛無觸發
    落在自己窗口內的測站數，再對事件加總（單尾）。"""
    if tested.empty:
        return {"n_events": 0}
    obs = int((tested.status.eq("onset") & tested.get("in_arrival_window", pd.Series(False, index=tested.index))
               .fillna(False).astype(bool)).sum())
    counts, indep_var = [], 0.0
    for (g, d), ev in tested.groupby(["group", "event_date"], sort=False):
        nl = nulls.get((g, d), {})
        arr = [nl[st] for st in ev.station if st in nl and len(nl[st])]
        if arr and len({len(a) for a in arr}) == 1:
            counts.append(np.sum(arr, axis=0))
            q = np.mean(arr, axis=1)
            indep_var += float(np.sum(q * (1 - q)))
    if not counts:
        return {"n_events": 0}
    sims = np.zeros(N_SIM)
    for c in counts:
        sims += c[rng.integers(0, len(c), N_SIM)]
    event_var = float(sum(c.var() for c in counts))
    return {"n_events": len(counts), "n_stations_tested": int(len(tested)), "n_in_window": obs,
            "expected_by_chance": round(float(sum(c.mean() for c in counts)), 2),
            "null_sd": round(float(sims.std()), 2),
            # 事件加總的虛無變異數相對於測站獨立時應有的值：
            # > 1 代表同一事件的測站會一起觸發（外部擾動）
            "dispersion_ratio": round(event_var / indep_var, 3) if indep_var > 0 else None,
            "p_one_sided": round(float((1 + np.sum(sims >= obs)) / (N_SIM + 1)), 5)}


def timing_test_dispersion_adjusted(trig: pd.DataFrame, dispersion: float | None) -> dict:
    """逐觸發的時間問題（已知一個測站觸發了，它的
    起始時間落在到時窗的頻率是否比它自己的虛無
    觸發更高？），Poisson-二項變異數乘上事件層級的
    離散比，讓同一事件一起觸發的測站不會被
    當成獨立證據。常態近似，單尾。"""
    if trig.empty or not dispersion:
        return {"n_onsets": int(len(trig))}
    q = trig.null_q_arrival.to_numpy(dtype=float)
    obs = int(trig.in_arrival_window.astype(bool).sum())
    sd = float(np.sqrt(max(dispersion, 1.0) * np.sum(q * (1 - q))))
    z = (obs - q.sum()) / sd if sd > 0 else float("nan")
    return {"n_onsets": int(len(trig)), "n_in_window": obs, "expected_by_chance": round(float(q.sum()), 2),
            "variance_inflation": round(max(dispersion, 1.0), 3), "z": round(z, 2),
            "p_one_sided_normal": float(f"{norm.sf(z):.3g}") if z == z else None}


def trigger_rate_by_distance(onsets: pd.DataFrame) -> list[dict]:
    tested = onsets[onsets.status.isin(["onset", "no_onset"])]
    out = []
    for b, d in tested.groupby(pd.cut(tested.hypocentral_km, [0, 30, 60, 100, 150, 250, 1000]), observed=False):
        if d.empty:
            continue
        out.append({"hypocentral_km": f"{b.left:g}-{b.right:g}", "n": int(len(d)),
                    "trigger_rate": round(float((d.status == "onset").mean()), 4),
                    "median_null_false_rate": round(float(d.null_false_rate.median()), 4)})
    return out


def seismometer_crosscheck(onsets: pd.DataFrame) -> dict:
    if not COMPARISON_CSV.exists():
        return {"note": "comparison_summary.csv missing"}
    cmp_ = pd.read_csv(COMPARISON_CSV)
    cmp_ = cmp_[(cmp_.status == "ok") & cmp_.shaking_onset_lag_sec.notna()]
    m = onsets[onsets.status == "onset"].merge(
        cmp_[["group", "date", "geomag_station", "seismic_station", "seismic_colocation_km", "shaking_onset_lag_sec"]],
        left_on=["group", "event_date", "station"], right_on=["group", "date", "geomag_station"])
    # 地震儀波形從發震前 60 秒開始；起始時間落在波形開頭是被截斷的挑選，不是到時
    m = m[m.shaking_onset_lag_sec > -55]
    return {"all_triggers": _crosscheck_stats(m),
            # xcg 在約 80% 的隨機時間都會觸發，所以它的「起始」大多只是窗口中第一個雜訊樣本
            f"null_false_rate_lt_{CROSSCHECK_FALSE_RATE:g}": _crosscheck_stats(m[m.null_false_rate < CROSSCHECK_FALSE_RATE])}


def _crosscheck_stats(m: pd.DataFrame) -> dict:
    if m.empty:
        return {"n": 0}
    dt = m.onset_lag_sec - m.shaking_onset_lag_sec
    return {"n": int(len(m)), "n_m6": int((m.magnitude >= 6).sum()),
            "median_geomag_minus_seismic_sec": round(float(dt.median()), 2),
            "iqr_sec": [round(float(dt.quantile(0.25)), 2), round(float(dt.quantile(0.75)), 2)],
            "frac_within_5s": round(float((dt.abs() <= 5).mean()), 4),
            "frac_geomag_earlier_by_gt5s": round(float((dt < -5).mean()), 4)}


def summarize(onsets: pd.DataFrame, per_event: pd.DataFrame, nulls: dict) -> dict:
    trig = onsets[onsets.status == "onset"]
    m6, clean_prior = onsets.magnitude >= 6, ~onsets.prior_event_shaking.astype(bool)
    subsets = {"m6": m6, "m5": ~m6,
               "m6_no_prior_shaking": m6 & clean_prior, "m5_no_prior_shaking": ~m6 & clean_prior}
    excluded = onsets[~clean_prior].drop_duplicates(["group", "event_date"])
    out = {"seed": SEED, "k_sigma": K_SIGMA, "search_lags_sec": list(SEARCH_LAGS),
           "v_p_kms": V_P_KMS, "v_s_kms": V_S_KMS, "n_perm": N_PERM,
           "prior_event_exclusion_sec": PRIOR_EVENT_EXCLUSION_SEC,
           "prior_event_shaking_events": [f"{r.group} {r.event_date} M{r.magnitude:g} ({r.prior_event_sec:.0f} s)"
                                          for r in excluded.itertuples()],
           "n_station_tests": int(onsets.status.isin(["onset", "no_onset"]).sum()),
           "n_onsets": int(len(trig)),
           "overall_null_false_rate": round(float(onsets.null_false_rate.mean()), 4),
           "subsets": {}}
    for name, mask in subsets.items():
        t = trig[mask[trig.index]]
        tested = onsets[mask & onsets.status.isin(["onset", "no_onset"])]
        ev_all = arrival_window_event_level(tested, nulls, keyed_rng("onset", name, "window_event_all"))
        ev_clean = arrival_window_event_level(tested[tested.null_false_rate < CLEAN_FALSE_RATE], nulls,
                                              keyed_rng("onset", name, "window_event_clean"))
        pe = per_event[(per_event.magnitude >= 6) if name.startswith("m6") else (per_event.magnitude < 6)]
        if name.endswith("no_prior_shaking"):
            pe = pe[~pe.prior_event_shaking]
        out["subsets"][name] = {
            "n_events": int(len(pe)),
            "event_verdicts": pe.verdict.value_counts().to_dict(),
            # 主要：事件層級虛無分布（同一事件的測站共用參考時間）
            "arrival_window_event_level_all": ev_all,
            "arrival_window_event_level_clean": ev_clean,
            # 只看時間（以有觸發為條件），變異數乘上事件層級離散比
            "timing_dispersion_adjusted_all": timing_test_dispersion_adjusted(t, ev_all.get("dispersion_ratio")),
            "timing_dispersion_adjusted_clean": timing_test_dispersion_adjusted(
                t[t.null_false_rate < CLEAN_FALSE_RATE], ev_clean.get("dispersion_ratio")),
            "n_events_shared_null_fallback": int(onsets[mask & onsets.shared_null_fallback.eq(True)]
                                                 .drop_duplicates(["group", "event_date"]).shape[0]),
            # 參考：逐觸發 Poisson-二項（把測站當成獨立）
            "arrival_window_all_triggers": arrival_window_test(t, keyed_rng("onset", name, "window_all")),
            "arrival_window_clean_stations": arrival_window_test(t[t.null_false_rate < CLEAN_FALSE_RATE],
                                                                 keyed_rng("onset", name, "window_clean")),
            "pooled_all_triggers": pooled_moveout(t, keyed_rng("onset", name, "pooled_all")),
            "pooled_clean_stations": pooled_moveout(t[t.null_false_rate < CLEAN_FALSE_RATE],
                                                    keyed_rng("onset", name, "pooled_clean")),
            "trigger_rate_by_distance": trigger_rate_by_distance(onsets[mask]),
        }
    out["seismometer_crosscheck"] = seismometer_crosscheck(onsets)
    return out


def plot(onsets: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    trig = onsets[onsets.status == "onset"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    r = np.linspace(0, max(300, trig.hypocentral_km.max() if len(trig) else 300), 50)
    for ax, (label, mask) in zip(axes, [("M ≥ 6", trig.magnitude >= 6), ("M < 6", trig.magnitude < 6)]):
        d = trig[mask]
        clean = d.null_false_rate < CLEAN_FALSE_RATE
        ax.scatter(d.hypocentral_km[clean], d.onset_lag_sec[clean], s=14, label="station false rate < 5%")
        ax.scatter(d.hypocentral_km[~clean], d.onset_lag_sec[~clean], s=14, marker="x", alpha=0.6,
                   label="noisier station")
        ax.plot(r, r / V_P_KMS, "--", lw=1, label=f"P ({V_P_KMS:g} km/s)")
        ax.plot(r, r / V_S_KMS, ":", lw=1.5, label=f"S ({V_S_KMS:g} km/s)")
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_title(f"{label} (n={len(d)})")
        ax.set_xlabel("hypocentral distance (km)")
    axes[0].set_ylabel("geomagnetic onset lag after origin (s)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def run_all(group_ids: list[str] | None, min_mag: float | None) -> dict:
    group_ids = group_ids or list(GROUPS)
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        results = list(ex.map(process_group, group_ids, [min_mag] * len(group_ids)))
    onsets = pd.DataFrame([r for rs, _ in results for r in rs])
    nulls = {k: v for _, ns in results for k, v in ns.items()}
    per_event = per_event_moveout(onsets)
    summary = summarize(onsets, per_event, nulls)
    summary.update({"groups": group_ids, "min_mag": min_mag})
    out_dir = OUT_DIR if group_ids == list(GROUPS) and min_mag is None else \
        OUT_DIR / ("subset_" + "_".join(group_ids) + (f"_m{min_mag:g}" if min_mag is not None else ""))
    out_dir.mkdir(parents=True, exist_ok=True)
    onsets.to_csv(out_dir / "onsets.csv", index=False)
    per_event.to_csv(out_dir / "moveout_per_event.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str))
    plot(onsets, out_dir / "travel_time.png")
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=str), file=sys.stderr)
    return summary


# ---------------------------------------------------------------------------
# 自我測試
# ---------------------------------------------------------------------------

def _synthetic(delays: list[float], rng: np.random.Generator) -> list[int | None]:
    out = []
    for dly in delays:
        v = np.cumsum(rng.normal(0, 0.02, 3600))
        p0 = 2000
        s = p0 + int(round(dly))
        v[s:s + 60] += np.cumsum(rng.normal(0, 0.3, 60))
        r = detect_onset(v, p0, SEARCH_LAGS[1])
        out.append(r.get("onset_lag_sec"))
    return out


def _correlated_null_trial(rng: np.random.Generator, n_ev: int = 30, n_st: int = 6, k: int = N_NULL):
    """一份 H0 資料集，所有參考時間中有一半帶有全網
    擾動（每個測站在同一延遲觸發），其他時候每個
    測站有 10% 的機率單獨觸發——發震時間只是另一個
    參考時間。回傳形狀和 process_group 相同的 (rows, nulls)。"""
    rows, nulls = [], {}
    for e in range(n_ev):
        R = rng.uniform(20, 300, n_st)
        lo, hi = R / V_P_KMS + ARRIVAL_PAD_SEC[0], R / V_S_KMS + ARRIVAL_PAD_SEC[1]
        lag = np.full((n_st, k + 1), np.nan)
        common_t = rng.random(k + 1) < 0.5
        lag[:, common_t] = rng.uniform(SEARCH_LAGS[0], SEARCH_LAGS[1], common_t.sum())
        own = (rng.random((n_st, k + 1)) < 0.1) & ~common_t
        lag[own] = rng.uniform(SEARCH_LAGS[0], SEARCH_LAGS[1], own.sum())
        inw = (lag >= lo[:, None]) & (lag <= hi[:, None])
        nulls[("S", str(e))] = {f"s{i}": inw[i, 1:] for i in range(n_st)}
        for i in range(n_st):
            on = lag[i, 1:][~np.isnan(lag[i, 1:])]
            trig = not np.isnan(lag[i, 0])
            rows.append({"group": "S", "event_date": str(e), "station": f"s{i}",
                         "status": "onset" if trig else "no_onset", "hypocentral_km": R[i],
                         "onset_lag_sec": lag[i, 0], "in_arrival_window": bool(inw[i, 0]) if trig else None,
                         "null_q_arrival": float(np.mean((on >= lo[i]) & (on <= hi[i]))) if len(on) else 0.0})
    return pd.DataFrame(rows), nulls


def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    R = np.array([20, 45, 80, 120, 160, 210, 260], dtype=float)
    ok = True
    for name, delays, want in [("simultaneous", np.full(len(R), 12.0), "simultaneous"),
                               ("S-wave moveout", R / V_S_KMS, "seismic_moveout")]:
        on = _synthetic(list(delays), rng)
        if any(o is None for o in on):
            print(f"[self-test] {name}：漏掉的起始 {on}  FAIL", file=sys.stderr)
            ok = False
            continue
        s, _, lo, hi = theilslopes(on, R, alpha=CI_ALPHA)
        verdict = classify_slope(lo, hi)
        err = max(abs(o - d) for o, d in zip(on, delays))
        good = verdict == want and err <= 6  # 逐漸增強的爆發：最初幾個樣本可能低於 K_SIGMA
        print(f"[self-test] {name}: slope={s:.3f} s/km CI=[{lo:.3f},{hi:.3f}] -> {verdict}, "
              f"最大起始誤差 {err:.1f}s  {'PASS' if good else 'FAIL'}", file=sys.stderr)
        ok &= good
    quiet = [detect_onset(np.cumsum(rng.normal(0, 0.02, 3600)), 2000, SEARCH_LAGS[1]) for _ in range(200)]
    fr = np.mean([q["status"] == "onset" for q in quiet])
    print(f"[self-test] 純雜訊誤觸發率 {fr:.3f}  {'PASS' if fr < 0.02 else 'FAIL'}", file=sys.stderr)
    ok &= fr < 0.02
    fp_event = fp_station = 0
    n_rep = 200
    for _ in range(n_rep):
        df, nulls = _correlated_null_trial(rng)
        fp_event += arrival_window_event_level(df, nulls, rng)["p_one_sided"] < 0.05
        fp_station += arrival_window_test(df[df.status == "onset"], rng)["p_one_sided"] < 0.05
    good = fp_event / n_rep <= 0.08 and fp_station > fp_event
    print(f"[self-test] 相關的 H0（全網擾動），5% 水準下的偽陽性："
          f"事件層級 {fp_event / n_rep:.3f}，逐測站 {fp_station / n_rep:.3f}  "
          f"{'PASS' if good else 'FAIL'}", file=sys.stderr)
    ok &= good
    t0 = pd.Timestamp("2024-04-03 00:00:00")
    gaps = [prior_event_gap(t0 + pd.Timedelta(seconds=d), [t0]) for d in (PRIOR_EVENT_EXCLUSION_SEC,
                                                                        PRIOR_EVENT_EXCLUSION_SEC + 1, 0)]
    flags = [g is not None and g <= PRIOR_EVENT_EXCLUSION_SEC for g in gaps]
    good = flags == [True, False, False]
    print(f"[self-test] 900/901/0 秒時的前一事件旗標：{flags}  {'PASS' if good else 'FAIL'}", file=sys.stderr)
    ok &= good
    return bool(ok)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--all", action="store_true", help="在真實資料上執行（預設動作）")
    ap.add_argument("--group", action="append", help="限定在這一組（可重複）")
    ap.add_argument("--min-mag", type=float, default=None)
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not self_test():
        print("[main] 自我測試失敗——在碰真實資料之前中止", file=sys.stderr)
        sys.exit(1)
    run_all(args.group, args.min_mag)
