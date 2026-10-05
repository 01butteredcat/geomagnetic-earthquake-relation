"""驗證步驟（計畫第 7 節）——在相信這個流程得出的任何前兆結論之前，
必須執行並檢視所有檢查。

用法：verify_pipeline.py --group G10
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import DATA_TIMEZONE, auto_outage_dates, list_day_refs, load_group_config, resolve_day_ref  # noqa: E402
from parser import is_scalar_only, open_raw, parse_day_file  # noqa: E402
from build_daily_features import SPIKE_THRESHOLD_NT  # noqa: E402
from compute_indices import MIN_CLEAN_POINTS, MIN_NIGHT_MINUTES, NIGHT_HOURS_UTC, TRAILING_WINDOW_DAYS  # noqa: E402
from stat_utils import bootstrap_ci  # noqa: E402

random.seed(42)
checks: dict = {}


def check_raw_vs_parsed_spotcheck(cfg):
    samples = []
    all_stations = list(cfg.stations.keys())
    for _ in range(4):
        station = random.choice(all_stations)
        refs = list_day_refs(cfg.gdms_dir, station)
        ref = random.choice(refs)
        df = parse_day_file(ref, station)
        scalar = is_scalar_only(ref)
        row_idx = random.randint(0, len(df) - 1)
        parsed_row = df.iloc[row_idx]
        ts = df.index[row_idx]

        raw_line = None
        with open_raw(ref) as f:
            for line in f:
                if line.startswith(ts.strftime("%Y-%m-%d %H:%M:%S")):
                    raw_line = line
                    break
        assert raw_line is not None, f"找不到 {ts} 的原始資料列：{ref.label}"
        parts = raw_line.split()
        raw = dict(zip(["X", "Y", "Z", "F"], (float(v) for v in parts[3:7])))

        real_cols = ["F"] if scalar else ["X", "Y", "Z"]
        junk_cols = ["X", "Y", "Z"] if scalar else ["F"]
        # 垃圾欄位（對這個檔案的通道類型來說從來不是真實資料）
        # 通常是未回報佔位值（88888），但在一些較舊的
        # （2024 年前、CWB 時代）完全沒資料的日子，檔案會改把
        # 中斷哨兵值（99999）也填進垃圾欄位——例如整個
        # G2_G3 時代的檔案，向量型測站在完全沒有資料的那天
        # 可能讀到 X=Y=Z=88888/F=99999。流程本身
        # 已經正確處理這種情況（parser.py 的防禦性 NaN 處理
        # 也把真實欄位中零星的 88888 當成無效），所以這個
        # 檢查只要求垃圾值是兩個已知哨兵值**之一**，
        # 不指定是哪一個——真正錯誤／損壞的原始
        # 值仍然會被抓到。
        junk_ok = all(raw[c] in (88888.0, 99999.0) for c in junk_cols)
        # 真實欄位在解析輸出中，遇到中斷
        # 哨兵值（99999）一律是 NaN，遇到零星的未回報佔位值
        # （88888）也是——但只限向量（X/Y/Z）欄位，那是
        # parser.py 防禦性備援套用的地方（完全沒資料的日子
        # 顯然仍可能把 88888 寫進本該是真實向量的
        # 欄位；依 common.py／NOTES.md 的慣例，實際資料
        # 不應該發生這種事，但流程防禦性地把它當成
        # 無效，而不是一個假的真實讀值）。純量（F）欄位
        # 在 parser.py 中沒有這種防禦性對應，所以那裡的 88888 不在
        # 預期／檢查範圍內。
        nan_sentinels = {"X": (88888.0, 99999.0), "Y": (88888.0, 99999.0), "Z": (88888.0, 99999.0), "F": (99999.0,)}
        real_ok = all(
            pd.isna(parsed_row[c]) if raw[c] in nan_sentinels[c] else np.isclose(parsed_row[c], raw[c])
            for c in real_cols
        )
        ok = junk_ok and real_ok
        samples.append({"station": station, "file": ref.label, "timestamp": str(ts), "ok": bool(ok)})

    all_ok = all(s["ok"] for s in samples)
    checks["raw_vs_parsed_spotcheck"] = {"pass": all_ok, "samples": samples}


def check_known_outage_ground_truth(cfg):
    """人工預先確認的三個特定 twu 中斷日的實際資料
    ——只對 G10 有意義，那是在流程最初建立時對照原始
    檔案人工檢查過的。其他組沒有對應的
    人工確認實際資料，所以這個檢查在那裡不適用
    （它們的缺口處理正確性改由下面的
    check_no_spikes_adjacent_to_missing_data 通用地涵蓋）。"""
    if cfg.group_id != "G10":
        checks["known_outage_ground_truth"] = {
            "pass": None,
            "note": "skipped_not_applicable——這份人工確認的實際資料只適用於 G10",
        }
        return

    results = {}
    df419 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240419"), "twu")
    n_nan_419 = int(df419["X"].isna().sum())
    results["twu_20240419_nan_count_expected_404"] = n_nan_419

    df424 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240424"), "twu")
    n_nan_424 = int(df424["X"].isna().sum())
    results["twu_20240424_nan_count_expected_86400"] = n_nan_424

    df425 = parse_day_file(resolve_day_ref(cfg.gdms_dir, "twu", "20240425"), "twu")
    n_nan_425 = int(df425["X"].isna().sum())
    results["twu_20240425_nan_count_expected_86400"] = n_nan_425

    ok = n_nan_419 == 404 and n_nan_424 == 86400 and n_nan_425 == 86400
    checks["known_outage_ground_truth"] = {"pass": ok, **results}


def check_timezone_conclusion(cfg):
    """common.DATA_TIMEZONE（'UTC'）是已經確立、全資料集共用的
    實際資料，在 G10 乾淨的 2024 資料上以高信心確認
    （見 common.py 的註解）——它應該對每一組都成立，因為所有
    組別都使用同一個 CWA IAGA-2002 紀錄慣例。這個檢查
    逐組重跑 timezone_check.py 是抽查，不是從頭重新
    證明：資料品質較差的測站（例如 G11 的
    twu，有很多部分中斷）合理地可能得到
    UNCERTAIN／低信心的結果，那不代表有任何問題——
    只是那個測站這次的訊號不夠乾淨。真正的問題是
    一個**有信心**、卻和
    common.DATA_TIMEZONE 矛盾的結論；那才是真正的警訊。"""
    tz = json.loads((cfg.interim_dir / "timezone_check.json").read_text())
    conclusion, confidence = tz["conclusion_timezone"], tz["confidence"]
    expected = "UTC" if DATA_TIMEZONE == "UTC" else "LOCAL (UTC+8)"
    if confidence != "high":
        result_pass = None
    else:
        result_pass = conclusion == expected
    checks["timezone_conclusion"] = {
        "pass": result_pass,
        "conclusion": conclusion,
        "confidence": confidence,
        "expected_per_common_py": expected,
        "note": (
            "檔頭沒有時區欄位；common.DATA_TIMEZONE 是全資料集共用的實際資料"
            "（已在 G10 的乾淨資料上確認）——這個檢查只有在出現「有信心」的矛盾時才失敗，"
            "雜訊較大的測站／組別得到低信心或無法判定的抽查結果時不會失敗"
        ),
    }


# 磁暴消除檢定（多日）。遠站->近站迴歸只在
# 非磁暴日上擬合，所以每個磁暴日都是外推，而只看
# 單一天的判定（舊的「Kp 最高那天」版本）會被那天的
# 雜訊主導：G2/G3 的 F 指標在那天只有約 -2 z，G6/G13 的近站和遠站
# 那天就是不一致。改成合併**每一個**磁暴開始日
# （Kp>=KP_STORM_THRESHOLD 或 Dst<=DST_STORM_THRESHOLD，取自 storm_days.csv）
# 中指標可計算的日子，比較該組典型磁暴日的
# |near_index| 和迴歸後典型的 |local_anomaly_index|。
# 磁暴日排除在迴歸擬合之外，所以這個檢定是樣本外的：
# 擬合天數少（G13：18 天）不會讓強的合併結果失效，只會
# 得到一個警告（WARN_FIT_DAYS_BELOW）——判定取決於磁暴日。
WARN_FIT_DAYS_BELOW = 30
MIN_STORM_DAYS_FOR_TEST = 5  # 可用磁暴日少於這個：不判定
STORM_SIGNAL_FLOOR = 0.5  # |near_index| 中位數（z）低於這個：沒有可量測的磁暴訊號可以消除
STORM_CANCELLATION_RATIO = 0.6  # 中位數|local| < 這個 * 中位數|near| 則通過（和舊的單日規則同樣是 0.6；判斷值）
# 只看點估計判定，會讓運氣好的小樣本比值通過，
# 即使重抽的不確定性顯示它並沒有可靠地低於
# 門檻——現在通過需要 bootstrap 信賴區間的**上界**低於
# STORM_CANCELLATION_RATIO，而不只是點估計（嚴格收緊：
# ci_hi >= point_estimate 永遠成立，所以這只會把通過變成失敗，
# 不會反過來）。LOW_BOOTSTRAP_POWER_THRESHOLD 標記（不擋）天數
# 少到 bootstrap 沒有多少組合空間的情況：n 個項目的
# 多重集合重抽只有 C(2n-1, n) 種不同結果（n=5 時 126，
# n=10 時約 92k，n=15 時約 7.76e7），所以低於約 15 時，信賴區間比
# N_BOOTSTRAP=2000 次抽樣所暗示的更粗。
CI_LEVEL = 0.90
N_BOOTSTRAP = 2000  # 和這個 repo 中其他所有 bootstrap 腳本相同（下面的 SEED 也是）——這裡的成本反正可以忽略
SEED = 20260805
LOW_BOOTSTRAP_POWER_THRESHOLD = 15


def _median_ratio(near_abs: np.ndarray, local_abs: np.ndarray) -> float:
    """給 bootstrap_ci 的 statistic()：一次配對重抽上的
    中位數|local| / 中位數|near|。當某次重抽的近站日中位數
    塌成 0 時回傳 NaN（而不是拋出例外），讓 nanpercentile 丟掉那個罕見的
    退化重抽，而不是讓整個 bootstrap 當掉。"""
    med_near = np.median(near_abs)
    return float(np.median(local_abs) / med_near) if med_near > 0 else float("nan")


def check_storm_cancellation(cfg):
    idx_field = "H" if cfg.xyz_pool.sufficient else ("F" if cfg.f_pool.sufficient else None)
    storm_path = cfg.interim_dir / "storm_days.csv"
    if not storm_path.exists() or idx_field is None:
        checks["storm_cancellation_test"] = {
            "pass": None,
            "note": "無法判定——沒有磁暴日清單，或沒有足以計算指標的測站池",
        }
        return

    idx = pd.read_csv(cfg.interim_dir / "local_anomaly_index.csv", dtype={"date": str})
    near_col, local_col = f"near_index_{idx_field}", f"local_anomaly_index_{idx_field}"
    if near_col not in idx.columns:
        checks["storm_cancellation_test"] = {"pass": None, "note": f"無法判定——{near_col} 沒有算出來"}
        return

    storm = pd.read_csv(storm_path, dtype={"date": str})
    storm_onsets = set(storm.loc[storm["is_storm_onset"], "date"])
    rows = idx[idx.date.isin(storm_onsets)].dropna(subset=[near_col, local_col])
    near_abs = rows[near_col].abs().to_numpy()
    local_abs = rows[local_col].abs().to_numpy()

    cand_path = cfg.interim_dir / "candidate_windows.json"
    cand = json.loads(cand_path.read_text()) if cand_path.exists() else {}
    n_fit = (cand.get(f"fit_{idx_field}") or {}).get("n_fit_days")
    summary_path = cfg.interim_dir / "storm_days_summary.json"
    storm_conf = json.loads(summary_path.read_text()).get("confidence") if summary_path.exists() else None

    result = {
        "pass": None,
        "field": idx_field,
        "n_storm_days_total": len(storm_onsets),
        "n_storm_days_used": int(len(rows)),
        "n_fit_days": n_fit,
        "storm_days_confidence": storm_conf,
        "note": "",
    }
    if len(rows):
        med_near, med_local = float(np.median(near_abs)), float(np.median(local_abs))
        near_v, local_v = rows[near_col].to_numpy(), rows[local_col].to_numpy()
        result.update(
            median_abs_near_index=med_near,
            median_abs_local_anomaly_index=med_local,
            ratio_local_over_near=(med_local / med_near) if med_near > 0 else None,
            frac_days_reduced=float(np.mean(local_abs < near_abs)),
            frac_days_sign_flipped=float(np.mean(np.sign(near_v) != np.sign(local_v))),
        )
        # 參考用：跨磁暴日 |local| < |near| 的單尾 Wilcoxon 符號等級檢定
        if len(rows) >= MIN_STORM_DAYS_FOR_TEST and np.any(near_abs != local_abs):
            from scipy.stats import wilcoxon

            result["wilcoxon_p_local_smaller"] = float(wilcoxon(near_abs, local_abs, alternative="greater").pvalue)

        if len(rows) >= MIN_STORM_DAYS_FOR_TEST:
            rng = np.random.default_rng(SEED)
            boot = bootstrap_ci((near_abs, local_abs), _median_ratio, N_BOOTSTRAP, rng, ci=CI_LEVEL)
            result.update(
                bootstrap_ci90_lo=boot["ci_lo"],
                bootstrap_ci90_hi=boot["ci_hi"],
                n_bootstrap=N_BOOTSTRAP,
                seed=SEED,
            )

    result_warnings = []
    if len(rows) < MIN_STORM_DAYS_FOR_TEST:
        result["note"] = (
            f"無法判定——只有 {len(rows)} 個磁暴日可計算指標（< {MIN_STORM_DAYS_FOR_TEST}）；"
            "太少，無法判斷"
        )
    elif result["median_abs_near_index"] < STORM_SIGNAL_FLOOR:
        result["note"] = (
            f"無法判定——磁暴日 |near_index| 中位數只有 {result['median_abs_near_index']:.2f} z"
            f"（< {STORM_SIGNAL_FLOOR}）；這組沒有可量測的磁暴訊號可以消除"
        )
    elif result.get("bootstrap_ci90_hi") is None or np.isnan(result["bootstrap_ci90_hi"]):
        result["note"] = (
            "無法判定——比值的 bootstrap 信賴區間沒有定義（有一部分退化的重抽"
            "中位數|near_index| == 0）；無法和比值門檻比較"
        )
    else:
        ci_hi = result["bootstrap_ci90_hi"]
        result["pass"] = bool(ci_hi < STORM_CANCELLATION_RATIO)
        result["note"] = (
            f"合併 {len(rows)} 個磁暴日：點估計 = {result['ratio_local_over_near']:.2f}，"
            f"{int(CI_LEVEL * 100)}% bootstrap 信賴區間 = [{result['bootstrap_ci90_lo']:.2f}, {ci_hi:.2f}]，"
            f"重抽 {N_BOOTSTRAP} 次（seed={SEED}）；通過需要信賴區間上界 < {STORM_CANCELLATION_RATIO}"
            "（比舊的點估計規則嚴格——點估計勉強通過的情況，如果信賴區間很寬，現在可能會失敗）；"
            "預期遠站迴歸會吸收共模的磁暴訊號"
        )
        if len(rows) < LOW_BOOTSTRAP_POWER_THRESHOLD:
            result_warnings.append(
                f"low_bootstrap_power——只用了 {len(rows)} 個磁暴日（< {LOW_BOOTSTRAP_POWER_THRESHOLD}）；"
                "這個 n 下重抽空間的組合數很少，所以信賴區間的邊界比較粗、比較不可信"
                "——這不會降低判定結果，只是對信賴區間精度的提醒"
            )
    if n_fit is None or n_fit < WARN_FIT_DAYS_BELOW:
        result_warnings.append(f"只有 {n_fit} 個迴歸擬合日（< {WARN_FIT_DAYS_BELOW}）；遠站->近站的擬合很薄")
    if result_warnings:
        result["warnings"] = result_warnings
    checks["storm_cancellation_test"] = result


def check_candidates_not_in_outage_windows(cfg):
    """候選日絕不應該落在 compute_indices.py 自己的
    is_clean_day 認為不乾淨的日子——candidate_flag 的定義要求
    is_clean_day==True，所以這是確認這道關卡真的
    守住的回歸測試，不是同義反覆：它獨立重建
    compute_indices.py 的 clean_overall 現在使用的同一個「中斷」事件
    （一天只有在測站池近站組的**每一個**測站、或遠站組的每一個
    測站那天都中斷時，才算該測站池的中斷——和近站中位數／
    遠站中位數指標可以容忍單一測站缺漏的方式一致），而不是
    較早、過度嚴格的「任何一個相關測站不乾淨」聯集，
    那在 auto_outage_dates 開始找到長期的
    單一測站中斷後產生了假失敗，而那些中斷近站／遠站中位數早就
    能毫無問題地繞過。測站–日是否「中斷」依夜間
    窗口判斷（見下面），和 compute_indices.py 的夜間特徵
    （MIN_NIGHT_MINUTES）實際把關的方式一致。"""
    cand = json.loads((cfg.interim_dir / "candidate_windows.json").read_text())
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    all_days = sorted(daily["date"].unique())

    # 這裡的「中斷」依實際建立指標用的**夜間**窗口判斷
    # （NIGHT_HOURS_UTC，該測站＋通道 >= MIN_NIGHT_MINUTES 個有效分鐘），
    # 而不是整天的 pct_missing：一個測站可能 24 小時內缺漏 >5%，
    # 夜間窗口卻完整，反之亦然。正是這種
    # 不一致讓這個檢查標記出 G14/G18 的候選日，而它們的近站或
    # 遠站中位數其實是用完整的夜間資料算出來的。
    night_ok_cache: dict[tuple[str, str], set[str]] = {}

    def night_ok_dates(station, channel):
        key = (station, channel)
        if key not in night_ok_cache:
            df = pd.read_parquet(cfg.interim_dir / f"minute_series_{station}.parquet", columns=[channel])
            night = df[df.index.hour.isin(NIGHT_HOURS_UTC)][channel].dropna()
            n = night.groupby(night.index.strftime("%Y%m%d")).size()
            night_ok_cache[key] = set(n[n >= MIN_NIGHT_MINUTES].index)
        return night_ok_cache[key]

    def pool_outage_dates(pool, channel):
        near_out = {d for d in all_days if pool.near and not any(d in night_ok_dates(s, channel) for s in pool.near)}
        far_out = {d for d in all_days if pool.far and not any(d in night_ok_dates(s, channel) for s in pool.far)}
        return near_out | far_out

    overlap = set()
    if cand.get("candidate_dates_H") or cand.get("candidate_dates_Z"):
        overlap |= set(cand.get("candidate_dates_H", [])) & pool_outage_dates(cfg.xyz_pool, "H")
        overlap |= set(cand.get("candidate_dates_Z", [])) & pool_outage_dates(cfg.xyz_pool, "Z")
    if cand.get("candidate_dates_F"):
        overlap |= set(cand.get("candidate_dates_F", [])) & pool_outage_dates(cfg.f_pool, "F")

    all_candidates = set(cand.get("candidate_dates_H", [])) | set(cand.get("candidate_dates_Z", [])) | set(cand.get("candidate_dates_F", []))
    checks["candidates_not_in_outage_windows"] = {
        "pass": len(overlap) == 0,
        "candidate_dates": sorted(all_candidates),
        "overlap_with_known_outages": sorted(overlap),
    }


def check_quiet_day_smoothness(cfg):
    """目視健全性檢查的代理：在時區檢查用的
    平靜日上，解析出的曲線不應該有**解析**造成的假象。
    我們不要求每一秒都接近平坦——短暫（1-2 個樣本）
    同時跳動、最高約 O(100nT) 在物理上是真實的（例如
    急始型脈衝），即使在
    其他方面平靜的日子也偶爾會出現；n_spikes==0 已經確認沒有任何東西超過
    流程其他地方使用的 300nT 結構性故障門檻（見
    build_daily_features.SPIKE_THRESHOLD_NT）。這個檢查只是
    確認沒有普遍的損壞（許多大跳動，像 twu 那樣）。"""
    tz = json.loads((cfg.interim_dir / "timezone_check.json").read_text())
    station = tz.get("station")
    quiet_days = tz["quiet_days_used"][:5]
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    sub = daily[(daily.station == station) & (daily.date.isin(quiet_days))]
    ok = bool((sub["n_spikes"] == 0).all()) if len(sub) else None
    checks["quiet_day_smoothness"] = {
        "pass": ok,
        "station": station,
        "days_checked": quiet_days,
        "max_abs_jump_values": sub["max_abs_jump"].tolist(),
        "note": "偶爾出現兩位數到約 O(100nT) 的單秒跳動是預期中的真實瞬變，不是解析 bug；n_spikes==0 確認沒有任何跳動超過 300nT 的結構性故障門檻",
    }


def check_no_spikes_adjacent_to_missing_data(cfg):
    """歷史 bug 的回歸測試：NaN 差分（來自
    中斷哨兵值缺口）曾被誤算成突波，那會扭曲
    n_spikes 並灌大表面上的異常數。build_daily_features.py 的
    _despike() 現在依靠「和 NaN 比較為 False」來排除缺口
    邊界（見它的 docstring）——這個測試**不**呼叫 _despike()
    本身（那會是同義反覆），而是獨立重新解析每個
    有任何缺漏樣本的測站／日，從頭重算差分／
    門檻檢查，再斷言沒有任何被標記的突波緊鄰
    NaN 樣本。"""
    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    partial = daily[(daily.n_valid > 0) & (daily.n_valid < daily.n_total)]

    bad = []
    for _, row in partial.iterrows():
        station, date = row["station"], row["date"]
        ref = resolve_day_ref(cfg.gdms_dir, station, date)
        df = parse_day_file(ref, station)
        for col in df.columns:  # 純量檔是 ["F"]，向量檔是 ["X","Y","Z"]
            s = df[col]
            d = s.diff().abs()
            flagged = np.where((d > SPIKE_THRESHOLD_NT).values)[0]
            for i in flagged:
                if pd.isna(s.iloc[i]) or pd.isna(s.iloc[i - 1]):
                    bad.append({"station": station, "date": date, "column": col, "row": int(i)})

    checks["no_spikes_adjacent_to_missing_data"] = {
        "pass": len(bad) == 0,
        "station_days_checked": len(partial),
        "bad_flags": bad[:20],
        "note": "這裡只要有任何條目，就代表有突波被標記在緊鄰 NaN 樣本的位置——正是 NaN 被誤算成突波的失效模式",
    }


def check_baseline_window_excludes_storms(cfg):
    """歷史 bug 的回歸測試（在 G10/2024 資料上發現）：
    太短的滑動基準窗口讓磁暴
    汙染（或完全餓死）在錨點事件前關鍵窗口中
    判斷候選異常日所用的參考。不依賴
    compute_indices.mad_zscore，獨立為關鍵窗口中的每一天重建
    TRAILING_WINDOW_DAYS 個日曆天的窗口（錨點事件日期
    -8..+1 天——這是對 G10 2024-04-03 主震
    和 2024-03-30 候選日有影響的偏移，因為
    失效模式——磁暴剛好落在窗口之前——是通用的，所以當成全資料集共用的常數重用）
    並檢查 (a) 那些窗口中至少有一個真的和實際磁暴
    日期重疊——確認這個測試確實觸及失效情境，
    而不是空洞的——以及 (b) 那些窗口中的每一個都仍
    有足夠的乾淨（非磁暴、非中斷）日來支撐真正的基準
    （>= MIN_CLEAN_POINTS，compute_indices.py 自己的門檻，低於它
    mad_zscore 會回傳 NaN，也就是完全沒有基準）。

    這裡的「中斷」指那天的 near_index（或 far_index）會是
    NaN——也就是**每一個**近站（或每一個遠站）那天都中斷
    ——和 compute_indices.py 的 clean_overall 定義一致（3 個近站
    取中位數可以容忍其中一個停擺；只有全部都停擺時
    才真的變成 NaN）。三個近站中只有一個
    停擺的日子，在這裡不算中斷。資料範圍以外的日子
    （第一個檔案之前）永遠不算乾淨。這個檢查不
    呼叫 compute_indices.py 本身（那會是同義反覆），而是
    從 daily_features.csv 獨立重建同樣的「整個測站池停擺」
    條件。"""
    storm_path = cfg.interim_dir / "storm_days.csv"
    if not storm_path.exists():
        checks["baseline_window_excludes_storms"] = {"pass": None, "note": "無法判定——缺少 storm_days.csv"}
        return
    storm_dates = set(pd.read_csv(storm_path, dtype={"date": str})["date"])

    daily = pd.read_csv(cfg.interim_dir / "daily_features.csv", dtype={"date": str})
    outage_by_station = auto_outage_dates(daily)

    pool = cfg.xyz_pool if cfg.xyz_pool.sufficient else cfg.f_pool
    all_days = sorted(daily["date"].unique())
    data_days = set(all_days)
    near_out = {d: all(d in outage_by_station.get(s, set()) for s in pool.near) for d in all_days}
    far_out = {d: all(d in outage_by_station.get(s, set()) for s in pool.far) for d in all_days}
    outage_dates = {d for d in all_days if near_out[d] or far_out[d]}

    anchor_date = pd.to_datetime(cfg.anchor_event.date)
    critical_window = pd.date_range(
        anchor_date - pd.Timedelta(days=8), anchor_date + pd.Timedelta(days=1)
    ).strftime("%Y%m%d").tolist()

    results = []
    for date_str in critical_window:
        d = pd.to_datetime(date_str, format="%Y%m%d")
        window = pd.date_range(
            d - pd.Timedelta(days=TRAILING_WINDOW_DAYS), d - pd.Timedelta(days=1)
        ).strftime("%Y%m%d").tolist()
        storm_in_window = sorted(set(window) & storm_dates)
        # 只有資料中存在的日子才能當基準日（G12 的資料只在錨點前 16 天
        # 開始；那之前的日子以前被算成乾淨日）
        clean_days = [w for w in window if w in data_days and w not in storm_dates and w not in outage_dates]
        results.append(
            {
                "probe_date": date_str,
                "trailing_window": [window[0], window[-1]],
                "storm_dates_in_window": storm_in_window,
                "n_clean_days_in_window": len(clean_days),
            }
        )

    overlap_exercised = any(r["storm_dates_in_window"] for r in results)
    all_sufficient = all(r["n_clean_days_in_window"] >= MIN_CLEAN_POINTS for r in results)
    checks["baseline_window_excludes_storms"] = {
        "pass": all_sufficient if overlap_exercised else None,
        "critical_window_checked": [critical_window[0], critical_window[-1]],
        "trailing_window_days": TRAILING_WINDOW_DAYS,
        "min_clean_points_threshold": MIN_CLEAN_POINTS,
        "probe_days": results,
        "note": (
            "檢查錨點事件前關鍵窗口中的每一天，在排除重疊的磁暴日之後，仍有 >= MIN_CLEAN_POINTS 個乾淨"
            "基準日" + ("" if overlap_exercised else
            "——無法判定：沒有任何探測日的窗口真的和磁暴日重疊，所以這個測試"
            "沒有觸及失效情境")
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    check_raw_vs_parsed_spotcheck(cfg)
    check_known_outage_ground_truth(cfg)
    check_timezone_conclusion(cfg)
    check_storm_cancellation(cfg)
    check_candidates_not_in_outage_windows(cfg)
    check_quiet_day_smoothness(cfg)
    check_no_spikes_adjacent_to_missing_data(cfg)
    check_baseline_window_excludes_storms(cfg)

    n_pass = sum(1 for c in checks.values() if c.get("pass") is True)
    n_fail = sum(1 for c in checks.values() if c.get("pass") is False)
    n_inconclusive = sum(1 for c in checks.values() if c.get("pass") is None)

    report = {
        "group": args.group,
        "summary": {"pass": n_pass, "fail": n_fail, "inconclusive": n_inconclusive},
        "checks": checks,
    }
    (cfg.interim_dir / "verification_report.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))

    if n_fail > 0:
        print(f"\n{n_fail} 項檢查失敗——修好之前不要相信任何前兆結論", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
