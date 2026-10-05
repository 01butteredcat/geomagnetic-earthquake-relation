"""聯合統計檢定，結合兩條到目前為止互相獨立的
同震證據線：

  - `coseismic_stacking_analysis.py`：跨事件疊加時，地磁 1Hz 資料
    是否在發震秒顯示一致的階躍／突波特徵？
    （完全沒說明這個特徵是真實的磁場
    變化，還是震動造成的儀器雜訊）
  - `seismometer_comparison.py`：在有獨立
    地震儀／加速度儀資料的事件子集中，每個事件特徵的時間
    是和真實地動對齊（`aligned_with_shaking`，和
    儀器雜訊一致），還是領先／比震動持續更久
    （`leads_shaking` / `persists_after_shaking_ends`，和
    真實機制一致）？

## 這支腳本問的問題

只限兩條線**都**有資料的事件，*疊加後的*
地磁特徵是否真的跟著地震儀推出的雜訊 vs 訊號
區分走？如果同震訊號大多是震動雜訊，被獨立
分類為 `aligned_with_shaking` 的事件疊加出來的特徵，
應該至少和 `leads_shaking`/`persists_after_shaking_ends`
事件一樣強（甚至可以說更強，因為雜訊會隨
地動振幅放大）。如果也存在真實的地球物理機制，
領先／持續組疊加出來的應該明顯更強及／或
比對齊組持續更久。

## 資料快照（2026-08-16，直接讀自
`data/interim/seismometer_comparison/comparison_summary.csv`——如果那個檔案
重新產生，請重新推導，不要寫死）

27 起事件有地震資料；其中 23 起也有可用的地磁資料
（`status == "ok"`；另外 4 起——G2_G3 的兩起、G4、G16——失敗並回報
`no_geomag_data`，因為它們近站自己那個日曆天的 `.sec` 檔不存在，
是真實的歷史缺口，不是 bug）。這 23 起中：
`aligned_with_shaking`=12、`leads_shaking`=9、`persists_after_shaking_ends`=2。
光是 `persists_after_shaking_ends`（n=2）太少，無法單獨疊加，
所以這支腳本用**兩組**拆分，而不是三組：
  - `noise_arm`  = alignment_verdict == "aligned_with_shaking"        （12）
  - `signal_arm` = alignment_verdict in ("leads_shaking",
                                          "persists_after_shaking_ends") （11）
這直接沿用 `seismometer_comparison.py::alignment_verdict` 自己的
框架（它的 docstring 已經把領先／持續歸在一起，視為「傾向
真實的地球物理機制」）。

（上面的數字來自最初 27 起事件那一輪。截至 2026-09-24，兩組
是 14 雜訊／12 訊號 = 26 起。2026-09-23 把登錄表補到 117
起事件——68 起來自 CWA GDMS 匯出檔的 M5.0-5.9 非錨點事件——沒有增加
任何組員，因為新事件都沒有地震儀資料；唯一新
可配對的 G11 2025-01-21b 是 `insufficient_data`。它仍會稍微影響
結果：新事件擴大了 G11/G12/G13/G19/G20 中 6 個組員
事件周圍的事件外排除範圍，改變了它們的基準。2026-09-24 重跑：
遠站 H 的差值在小數第三位變動，最小 p_tail 維持 0.069，最小 p_peak
從 0.095 -> 0.093。）

（2026-09-25：80 起事件的地震儀批次到齊後，全部判定共給出 94 起分組
事件（48 雜訊／46 訊號），其中 38 起 M>=6（18／20）。只跑 M>=6
用 --min-mag 6；M5 的標籤接近擲硬幣——見 run_all()。
全部事件：32 個 p 值中有 3 個 < 0.05（遠站 H 突波峰值 0.010 和尾段 0.048，
兩者都是**雜訊組**較強；近站 H step30 尾段 0.024，訊號組
較強），沒有一個撐過 Bonferroni 校正。只看 M>=6：沒有低於 0.05 的，最小 0.062。）

## 方法

1. 建立 `EventSeries` 物件（`coseismic_stacking_analysis.py` 自己的
   dataclass + 載入機制，未修改），只針對 23 起事件的
   允許清單，透過新的載入器（`load_event_series_for_events`），它仿照
   `load_all_event_series`，但限制在明確的 (group_id, event_date)
   配對，而不是一組中的每個事件——加在這裡而不是
   `coseismic_stacking_analysis.py` 本身，讓那支已經發布的
   腳本保持不動。
2. 對每個 (channel_type_pool, station_tier, stat_name) 組合（仿照
   `coseismic_stacking_analysis.py::run_group_ids` 自己的組合迴圈）：
   觀測 Delta = 兩組在**兩個**統計量上的差（不只看峰值，
   因為 `persists` 型訊號的區別特徵是
   震動後的**尾段**，不一定是峰值）：
     - `delta_peak` = peak_abs_z(訊號組) - peak_abs_z(雜訊組)
     - `delta_tail` = mean(|stack_mean|，lag>0)(訊號組) - 同上(雜訊組)
3. 虛無分布：標籤置換檢定——在同一批事件間打亂雜訊／訊號
   組標籤（兩組大小固定）
   N_PERM=2000 次，每次用輕量的
   只算疊加平均的輔助函式（`_lightweight_stack_delta`）重算兩個差值，跳過
   `stack_series()` 自己的 bootstrap 信賴區間 + 虛無帶（2000+1000 次重抽——
   在 2000 次置換抽樣的每一次都重做那個，會慢上
   好幾個數量級）。完整的 `stack_series()`（含
   bootstrap 信賴區間和虛無帶）每個組合仍然剛好呼叫兩次
   ——每組一次——純粹用來報告／畫出兩組各自的疊加
   形狀，不用於置換檢定本身。
4. p 值 = (1 + #{|delta_perm| >= |delta_obs|}) / (N_PERM + 1)，和這個
   程式庫其他所有重抽檢定一樣，採用對 |.| 的單尾
   慣例。

## 誠實的注意事項（事先說明，不是事後補述）

N=94 起事件分成 48/46（雜訊／訊號，2026-09-28；M>=6 更少）是
探索性的樣本數，不是高檢定力的。這裡的結果最多只能報告成有提示性，
絕不能單獨當作確認性的發現。

用法：
  coseismic_joint_analysis.py --self-test    # 只做合成資料健全性檢查
  coseismic_joint_analysis.py --all          # 真實資料，全部 16 個組合
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from events import GROUPS, get_group  # noqa: E402
from coseismic_step_analysis import (  # noqa: E402
    EXCLUSION_BUFFER_SEC,
    SEED,
    _build_channels,
    _load_station_days,
    _rank_stations_for_event,
    _step_statistic,
)
from coseismic_stacking_analysis import (  # noqa: E402
    BUFFER_SEC,
    STACK_HALF_SEC,
    STAT_NAMES,
    STATION_TIERS,
    EventSeries,
    _build_event_series,
    _off_event_baseline,
    _window_profile,
    stack_series,
)

N_PERM = 2000
NOISE_VERDICTS = ("aligned_with_shaking",)
SIGNAL_VERDICTS = ("leads_shaking", "persists_after_shaking_ends")

OUT_DIR = common.PROJECT_DIR / "data" / "interim" / "coseismic_joint_analysis"
COMPARISON_CSV = common.PROJECT_DIR / "data" / "interim" / "seismometer_comparison" / "comparison_summary.csv"


# ---------------------------------------------------------------------------
# 分組：讀 seismometer_comparison.py 自己的輸出，而不是
# 在這裡重算 alignment_verdict——唯一的真相來源。
# ---------------------------------------------------------------------------

def load_arm_assignment(csv_path: Path = COMPARISON_CSV, verdict_col: str = "alignment_verdict") -> dict[str, str]:
    """回傳 {"<group_id>__<event_date>": "noise"|"signal"}，涵蓋每一列
    status=="ok" 且判定屬於兩個已知組之一的資料。任何其他
    status（例如 no_geomag_data）或判定（insufficient_data）的資料列
    就不會出現在回傳的 dict 中——不是錯誤，只是從
    聯合分析中排除，和模組 docstring 寫的完全一樣。"""
    df = pd.read_csv(csv_path)
    arm_of: dict[str, str] = {}
    for _, row in df.iterrows():
        if row.get("status") != "ok":
            continue
        verdict = row.get(verdict_col)
        key = f"{row['group']}__{row['date']}"
        if verdict in NOISE_VERDICTS:
            arm_of[key] = "noise"
        elif verdict in SIGNAL_VERDICTS:
            arm_of[key] = "signal"
    return arm_of


# ---------------------------------------------------------------------------
# 限制在明確事件允許清單的 EventSeries 載入器（這裡新加的——
# coseismic_stacking_analysis.py::load_all_event_series 沒有動）。
# ---------------------------------------------------------------------------

def load_event_series_for_events(event_keys: list[tuple[str, str]]) -> dict[tuple[str, str, str], EventSeries]:
    """和 coseismic_stacking_analysis.py::
    load_all_event_series 相同的載入邏輯，限制在明確的 (group_id, event_date)
    允許清單，而不是所要求組別中的每個事件——需要這樣是因為
    兩組是逐事件定義的（來自 seismometer_comparison.py 的
    逐事件 alignment_verdict），不是逐組。"""
    wanted: dict[str, set[str]] = {}
    for group_id, event_date in event_keys:
        wanted.setdefault(group_id, set()).add(event_date)

    out: dict[tuple[str, str, str], EventSeries] = {}
    for group_id, dates in wanted.items():
        cfg = common.load_group_config(group_id)
        group = get_group(group_id)
        for event in group.events:
            if event.date not in dates:
                continue
            event_utc = pd.Timestamp(event.time_utc)
            for channel_type in ("F", "XYZ"):
                stations = _rank_stations_for_event(cfg, event, channel_type, len(STATION_TIERS))
                for rank, tier in enumerate(STATION_TIERS):
                    if rank >= len(stations):
                        continue
                    station, distance_km = stations[rank]
                    df, missing, wanted_dates = _load_station_days(cfg.gdms_dir, station, event_utc, BUFFER_SEC)
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
# 置換迴圈用的輕量、只算疊加平均的路徑（跳過
# stack_series() 自己的 bootstrap 信賴區間 + 虛無帶——原因見模組 docstring：
# 2000 次置換抽樣 x 那套機制會慢得太多）。
# ---------------------------------------------------------------------------

def _lightweight_stack_delta(events: list[EventSeries], stat_name: str) -> dict | None:
    lags = np.arange(-STACK_HALF_SEC, STACK_HALF_SEC + 1)
    matrix = []
    for es in events:
        if stat_name not in es.stat_arrays:
            continue
        med, mad = es.baselines[stat_name]
        profile = _window_profile(es.stat_arrays[stat_name], es.idx, es.event_utc, STACK_HALF_SEC)
        if np.all(np.isnan(profile)):
            continue
        matrix.append((profile - med) / mad)
    if not matrix:
        return None
    M = np.array(matrix)
    stack_mean = np.nanmean(M, axis=0)
    peak_abs_z = float(np.nanmax(np.abs(stack_mean)))
    tail_vals = stack_mean[lags > 0]
    tail_vals = tail_vals[~np.isnan(tail_vals)]
    if len(tail_vals) == 0:
        return None
    return {"peak_abs_z": peak_abs_z, "tail_mean_abs_z": float(np.mean(np.abs(tail_vals))), "n_events": len(matrix)}


def _arm_deltas(events: list[EventSeries], labels: np.ndarray, stat_name: str) -> tuple[float, float] | None:
    noise = [es for es, lab in zip(events, labels) if lab == "noise"]
    signal = [es for es, lab in zip(events, labels) if lab == "signal"]
    dn = _lightweight_stack_delta(noise, stat_name)
    ds = _lightweight_stack_delta(signal, stat_name)
    if dn is None or ds is None:
        return None
    return ds["peak_abs_z"] - dn["peak_abs_z"], ds["tail_mean_abs_z"] - dn["tail_mean_abs_z"]


def permutation_test(events: list[EventSeries], arm_labels: np.ndarray, stat_name: str,
                      rng: np.random.Generator, n_perm: int = N_PERM) -> dict:
    obs = _arm_deltas(events, arm_labels, stat_name)
    if obs is None:
        return {"error": "insufficient events in one or both arms for this stat/combo"}
    obs_peak, obs_tail = obs

    perm_peak, perm_tail = [], []
    for _ in range(n_perm):
        shuffled = rng.permutation(arm_labels)
        d = _arm_deltas(events, shuffled, stat_name)
        if d is None:
            continue
        perm_peak.append(d[0])
        perm_tail.append(d[1])
    perm_peak = np.array(perm_peak)
    perm_tail = np.array(perm_tail)

    def _p(obs_val, perm_vals):
        if len(perm_vals) == 0:
            return None
        return (1 + int(np.sum(np.abs(perm_vals) >= abs(obs_val)))) / (len(perm_vals) + 1)

    return {
        "delta_peak_abs_z": round(obs_peak, 4),
        "delta_tail_mean_abs_z": round(obs_tail, 4),
        "p_value_peak": _p(obs_peak, perm_peak),
        "p_value_tail": _p(obs_tail, perm_tail),
        "n_perm_used": len(perm_peak),
    }


# ---------------------------------------------------------------------------
# 合成資料自我測試（不用真實資料）：確認置換機制
# 本身 (a) 能偵測到真實注入的組間差異，(b) 沒有差異時
# 不會產生偽陽性——和 coseismic_stacking_analysis.py::self_test
# 相同的雙向健全性檢查模式。
# ---------------------------------------------------------------------------

def _make_synth_series(rng: np.random.Generator, amplitude: float, tag: str) -> EventSeries:
    n_samples = 2 * BUFFER_SEC + 1
    idx = pd.date_range("2024-01-01", periods=n_samples, freq="s")
    center_i = BUFFER_SEC
    noise = rng.normal(0, 1.0, size=n_samples)
    if amplitude > 0:
        jitter = int(rng.integers(-5, 6))
        noise[center_i + jitter:center_i + jitter + 20] += amplitude
    event_utc = idx[center_i]
    arr = _step_statistic(noise, 30)
    base = _off_event_baseline(arr, idx, [event_utc], EXCLUSION_BUFFER_SEC)
    assert base is not None, "自我測試的事件外基準計算失敗"
    return EventSeries(
        group_id="SYN", event_date=f"{tag}-{int(rng.integers(0, 10**6))}", anchor=True, magnitude="M0",
        station="SYN", distance_km=0.0, idx=idx, stat_arrays={"step30": arr}, baselines={"step30": base},
        event_utc=event_utc,
        valid_lo=idx[0] + pd.Timedelta(seconds=STACK_HALF_SEC),
        valid_hi=idx[-1] - pd.Timedelta(seconds=STACK_HALF_SEC),
        exclude_centers=[event_utc],
    )


def self_test() -> bool:
    rng = np.random.default_rng(SEED)
    n_each = 12
    ok = True

    # 情況 (a)：真實的組間差異（訊號組注入較強的訊號）
    noise_events = [_make_synth_series(rng, 2.0, "noise") for _ in range(n_each)]
    signal_events = [_make_synth_series(rng, 15.0, "signal") for _ in range(n_each)]
    events_a = noise_events + signal_events
    labels_a = np.array(["noise"] * n_each + ["signal"] * n_each)
    result_a = permutation_test(events_a, labels_a, "step30", rng)
    detected = (result_a.get("p_value_peak") is not None and result_a["p_value_peak"] < 0.05
                and result_a["delta_peak_abs_z"] > 0)
    status_a = "PASS" if detected else "FAIL"
    print(f"[self-test] 有真實差異的情況：delta_peak_abs_z={result_a.get('delta_peak_abs_z')} "
          f"p_value_peak={result_a.get('p_value_peak')}  {status_a}")
    ok = ok and detected

    # 情況 (b)：沒有真實組間差異（兩組振幅相同）——偽陽性檢查
    rng2 = np.random.default_rng(SEED + 1)
    events_b = ([_make_synth_series(rng2, 5.0, "noise") for _ in range(n_each)]
                + [_make_synth_series(rng2, 5.0, "signal") for _ in range(n_each)])
    labels_b = np.array(["noise"] * n_each + ["signal"] * n_each)
    result_b = permutation_test(events_b, labels_b, "step30", rng2)
    no_false_positive = result_b.get("p_value_peak") is not None and result_b["p_value_peak"] >= 0.05
    status_b = "PASS" if no_false_positive else "FAIL"
    print(f"[self-test] 沒有差異的情況（偽陽性檢查）："
          f"delta_peak_abs_z={result_b.get('delta_peak_abs_z')} p_value_peak={result_b.get('p_value_peak')}  {status_b}")
    ok = ok and no_false_positive

    return ok


# ---------------------------------------------------------------------------
# 真實資料的統籌
# ---------------------------------------------------------------------------

def run_all(min_mag: float | None = None, gated: bool = False) -> dict:
    """min_mag 把兩組都限制在規模至少那麼大的 events.py 事件，
    寫到另一個目錄（coseismic_joint_analysis_m<min_mag>）。
    自 2026-09-25 的 M5 批次以來，值得和完整集合一起跑：對
    小事件來說，磁力儀記錄到的大多是雜訊，而雜訊峰值落在
    ±180 秒搜尋窗口中任何地方，約有一半的機率會在震動開始之前，
    所以 M5 的 "leads_shaking" 標籤接近擲硬幣。

    gated 使用 seismometer_comparison.py 的 alignment_verdict_gated（只有
    該事件自己的 step30 異常 p < 0.05 時才有判定），寫到
    coseismic_joint_analysis_gated[_m<mag>]/。2026-09-25 這樣只剩 7 起事件
    （5 對齊、1 領先、1 持續）——太少無法檢定；跑這個是為了把
    這件事明確呈現，而不是要當成結果解讀。"""
    name = OUT_DIR.name + ("_gated" if gated else "") + ("" if min_mag is None else f"_m{min_mag:g}")
    out_dir = OUT_DIR.with_name(name)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "combos").mkdir(exist_ok=True)

    arm_of = load_arm_assignment(verdict_col="alignment_verdict_gated" if gated else "alignment_verdict")
    if min_mag is not None:
        mag = {f"{gid}__{e.date}": e.magnitude for gid, g in GROUPS.items() for e in g.events}
        arm_of = {k: v for k, v in arm_of.items() if mag.get(k, 0) >= min_mag}
    event_keys = [tuple(k.split("__", 1)) for k in arm_of]
    print(f"[joint] {len(event_keys)} 起事件有明確的分組 "
          f"(noise={sum(1 for v in arm_of.values() if v=='noise')}, "
          f"signal={sum(1 for v in arm_of.values() if v=='signal')})", file=sys.stderr)

    all_series = load_event_series_for_events(event_keys)
    rng = np.random.default_rng(SEED)

    channel_labels = {"XYZ": "H", "F": "F"}
    summary_rows: list[dict] = []
    run_summary = {"seed": SEED, "n_perm": N_PERM,
                   "arm_assignment": arm_of, "combos": []}

    for channel_type, ch_label in channel_labels.items():
        for tier in STATION_TIERS:
            events = [es for (ct, t, _key), es in all_series.items() if ct == channel_type and t == tier]
            if not events:
                continue
            labels = np.array([arm_of[f"{es.group_id}__{es.event_date}"] for es in events])
            if "noise" not in labels or "signal" not in labels:
                continue

            for stat_name in STAT_NAMES:
                combo_id = f"{tier}__{ch_label}__{stat_name}"
                perm_result = permutation_test(events, labels, stat_name, rng)

                noise_events = [es for es, lab in zip(events, labels) if lab == "noise"]
                signal_events = [es for es, lab in zip(events, labels) if lab == "signal"]
                noise_stack = stack_series(noise_events, stat_name, rng)
                signal_stack = stack_series(signal_events, stat_name, rng)

                result = {
                    "combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                    "channel": ch_label, "stat": stat_name,
                    "n_noise_events": len(noise_events), "n_signal_events": len(signal_events),
                    "permutation_test": perm_result,
                    "noise_arm_stack": noise_stack, "signal_arm_stack": signal_stack,
                }
                (out_dir / "combos" / f"joint__{combo_id}.json").write_text(json.dumps(result, indent=2, default=str))

                if "error" in perm_result:
                    print(f"[{combo_id}] {perm_result['error']}", file=sys.stderr)
                    run_summary["combos"].append({"combo_id": combo_id, "status": "error",
                                                   "error": perm_result["error"]})
                    continue

                print(f"[{combo_id}] n_noise={len(noise_events)} n_signal={len(signal_events)} "
                      f"delta_peak={perm_result['delta_peak_abs_z']} p_peak={perm_result['p_value_peak']} "
                      f"delta_tail={perm_result['delta_tail_mean_abs_z']} p_tail={perm_result['p_value_tail']}",
                      file=sys.stderr)
                row = {"combo_id": combo_id, "station_tier": tier, "channel_type_pool": channel_type,
                       "channel": ch_label, "stat": stat_name,
                       "n_noise_events": len(noise_events), "n_signal_events": len(signal_events), **perm_result}
                summary_rows.append(row)
                run_summary["combos"].append({"combo_id": combo_id, "status": "ok", **perm_result})

    run_summary["min_mag"] = min_mag
    run_summary["gated"] = gated
    pd.DataFrame(summary_rows).to_csv(out_dir / "joint_summary.csv", index=False)
    (out_dir / "all_joint_run_summary.json").write_text(json.dumps(run_summary, indent=2, default=str))
    print(f"[run] {len(summary_rows)} combos -> {out_dir}", file=sys.stderr)
    return run_summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--self-test", action="store_true", help="只跑合成資料健全性檢查")
    ap.add_argument("--all", action="store_true", help="跑真實資料的聯合分析（預設動作）")
    ap.add_argument("--gated", action="store_true",
                    help="分組取自 alignment_verdict_gated（只有顯著異常），輸出到另外的 _gated 目錄")
    ap.add_argument("--min-mag", type=float, default=None,
                    help="只用規模至少這麼大的事件，輸出到另外的 _m<mag> 目錄")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("[main] 合成資料自我測試失敗——在碰真實資料之前中止", file=sys.stderr)
        sys.exit(1)

    run_all(min_mag=args.min_mag, gated=args.gated)
