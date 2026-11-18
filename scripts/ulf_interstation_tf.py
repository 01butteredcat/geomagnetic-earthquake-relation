"""站間轉換函數（interstation transfer function, ISTF）差分：磁暴日能不能放回來？（探索性診斷）

回答使用者 2026-11-13 的問題：「如果直接在原始資料這一層把近站和遠站相減，磁暴日能不能
加回分析、讓樣本變多？」11-02 試過的兩種做法都是先算出每日指標再相減（Z/H 近減遠、
Han／Warden 日能量比值），磁暴都消不乾淨。這裡改在時間序列（頻域）這一層，用遠站的
X/Y/Z 預測近站的 Z，取殘差：全球性的外部源場（包括磁暴的 Pc3/Pc4）兩站都收得到，會被
預測掉；只有近站才有的局部訊號會留在殘差裡。

文獻：
  Harada, M., Hattori, K., & Isezaki, N. (2004). Transfer function approach to signal
    discrimination of ULF geomagnetic data. Physics and Chemistry of the Earth 29(4–9),
    409–417. doi:10.1016/j.pce.2004.03.002
  概念（以遠端參考站三分量的站間轉換函數預測本站，殘差 = 局部訊號）取自這篇；論文用
  小波轉換，本專案拿不到全文，**估計細節（FFT 分段、頻率格合併、最小平方）是本專案自訂**，
  採用大地電磁法常見的頻域最小平方寫法。

每站每日（夜間，同 ulf_analysis 的 UTC 17–19 時 = 當地 01:00–03:59）：
  1. 取 UTC 16:00–21:00（前後各多 1 小時給去趨勢），任一分量有缺值或 INTERMAGNET 尖峰
     （build_daily_features._despike）就整晚不用。
  2. X/Y/Z 減掉 1 小時置中滑動平均（同 ulf_analysis._detrend），截出 UTC 17:00–20:00。
  3. 切成 1200 秒、重疊 50% 的子段（17 段），加 Hann 窗做 FFT，只留 Pc4（45–150 秒）到
     Pc3（10–45 秒）範圍內的頻率格。頻譜快取在 data/interim/ulf_istf/spectra/<資料夾>/。

每一對（近站 n、遠站 f），每個頻率格 k（合併 k±2 共 5 格的訓練資料）：
  複數最小平方 Z_n = t_x X_f + t_y Y_f + t_z Z_f，只用訓練日估計。
  訓練日 = 兩站都有資料的非磁暴日，扣掉錨點震前 30 天窗口和緊接在前的安慰劑區塊
  （排名檢定會檢定的兩個區塊），避免轉換函數把可能的前兆一起預測掉。
  每日能量（頻帶內頻率格加總，子段平均）：
    res    = |Z_n − 預測|²    ← 本腳本的指標（ISTF 殘差）
    direct = |Z_n − Z_f|²     ← 逐秒直接相減（使用者原本問的做法，對照）
    raw    = |Z_n|²           ← 不相減（對照）
  各配對（raw 是各近站）先除以自己非磁暴日的中位數，每日再取配對中位數。

診斷（只列描述數字，不設門檻、不算前兆 p 值；判準由使用者看完再定）：
  能量類指標取 log10 後算（能量是乘法尺度）；現行 Z/H 近減遠（ulf_near_far_index.csv 的
  <band>_diff_zh）照 11-02 用原尺度。每組 × 頻帶 × 版本：
  - 磁暴／平靜 MAD 比、磁暴日中位數偏移（以平靜日 MAD 為單位），定義同
    surrogate_test.storm_kept_one；
  - 指標與每日最大 Kp 的 Spearman ρ；
  - 錨點震前窗口內磁暴日比例 vs 整條序列、窗口最高／最低日是不是磁暴日；
  - 樣本增益：震前窗口與安慰劑區塊可用天數，排除磁暴日 vs 保留；
  - ISTF 預測得多好：平靜日殘差能量 ÷ 原始能量的中位數（越接近 1 越沒預測到）。

輸出：
  data/interim/ulf_istf/spectra/<資料夾>/<測站>.npz  每站每晚的頻譜（快取）
  data/interim/<group>/ulf_istf_index.csv            每日指標（res／direct／raw × pc3／pc4）
  data/interim/<group>/ulf_istf_fits.json            每配對訓練日數、平靜日殘差占比
  data/interim/ulf_istf/summary.{json,md}            跨組診斷

用法：
  ulf_interstation_tf.py --self-test
  ulf_interstation_tf.py --group G10 [--from-cache]
  ulf_interstation_tf.py --all [--from-cache]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
import method_comparison as mc  # noqa: E402
from build_daily_features import _despike  # noqa: E402
from common import PROJECT_DIR, list_day_refs, load_group_config  # noqa: E402
from cross_group_analysis import pre_event_anchor, pre_event_window  # noqa: E402
from parser import parse_day_file  # noqa: E402
from surrogate_test import ULF_GROUPS  # noqa: E402
from ulf_analysis import BANDS, NIGHT_HOURS_UTC  # noqa: E402

FS = 1.0
PAD_START = pd.Timedelta(hours=min(NIGHT_HOURS_UTC) - 1)   # UTC 16:00
NIGHT_START = pd.Timedelta(hours=min(NIGHT_HOURS_UTC))     # UTC 17:00
NIGHT_END = pd.Timedelta(hours=max(NIGHT_HOURS_UTC) + 1)   # UTC 20:00（不含）
PAD_END = NIGHT_END + pd.Timedelta(hours=1)                # UTC 21:00（不含）
N_PAD = int((PAD_END - PAD_START).total_seconds())         # 18000
NIGHT = slice(int((NIGHT_START - PAD_START).total_seconds()), int((NIGHT_END - PAD_START).total_seconds()))
N_NIGHT = NIGHT.stop - NIGHT.start                          # 10800
SEG_LEN = 1200
SEG_STEP = 600
N_SEG = (N_NIGHT - SEG_LEN) // SEG_STEP + 1                 # 17
FREQS = np.fft.rfftfreq(SEG_LEN, 1 / FS)
KEEP = (FREQS >= BANDS["pc4"][0]) & (FREQS <= BANDS["pc3"][1])
KEPT_FREQS = FREQS[KEEP]
BAND_BINS = {b: (KEPT_FREQS >= lo) & (KEPT_FREQS < hi if b == "pc4" else KEPT_FREQS <= hi)
             for b, (lo, hi) in BANDS.items()}
SMOOTH_HALF = 2          # 每個頻率格合併 k±2 格估計轉換函數
MIN_TRAIN_DAYS = 10      # 配對的訓練日少於這個就不估
VERSIONS = ("res", "direct", "raw")
COMPONENTS = ("X", "Y", "Z")
N_WORKERS = 10
OUT_DIR = PROJECT_DIR / "data" / "interim" / "ulf_istf"
SPEC_DIR = OUT_DIR / "spectra"
_HANN = np.hanning(SEG_LEN)


# ---------------------------------------------------------------------------
# 每站每晚頻譜

def segment_spectra(night: np.ndarray) -> np.ndarray:
    """night：(N_NIGHT, 3) 已去趨勢的 X/Y/Z。回傳 (N_SEG, 頻率格數, 3) complex64。"""
    out = np.empty((N_SEG, KEEP.sum(), 3), dtype=np.complex64)
    for s in range(N_SEG):
        seg = night[s * SEG_STEP: s * SEG_STEP + SEG_LEN] * _HANN[:, None]
        out[s] = np.fft.rfft(seg, axis=0)[KEEP]
    return out


def night_block(df: pd.DataFrame, date_str: str) -> np.ndarray | str:
    """回傳 (N_NIGHT, 3) 去趨勢後的夜間 X/Y/Z，或略過原因。"""
    if "X" not in df.columns:
        return "scalar_only"
    day = pd.Timestamp(date_str)
    seg = df.loc[(df.index >= day + PAD_START) & (df.index < day + PAD_END), list(COMPONENTS)].astype(float)
    if len(seg) != N_PAD:
        return f"window_has_{len(seg)}_samples"
    if seg.isna().to_numpy().any():
        return "gap_in_window"
    if any(_despike(seg[c])[1].any() for c in COMPONENTS):
        return "spike_in_window"
    cols = []
    for c in COMPONENTS:  # 同 ulf_analysis._detrend：1 小時置中滑動平均
        s = seg[c]
        cols.append((s - s.rolling(window=3601, center=True, min_periods=1).mean()).to_numpy()[NIGHT])
    return np.column_stack(cols)


def _process_day(ref):
    try:
        df = parse_day_file(ref, ref.station)
    except Exception as exc:  # noqa: BLE001
        return ref.station, ref.date_str, str(exc)
    blk = night_block(df, ref.date_str)
    if isinstance(blk, str):
        return ref.station, ref.date_str, blk
    return ref.station, ref.date_str, segment_spectra(blk)


def _spec_path(cfg, station: str) -> Path:
    return SPEC_DIR / cfg.gdms_dir.name / f"{station}.npz"


def load_spectra(cfg, station: str) -> dict[str, np.ndarray]:
    path = _spec_path(cfg, station)
    if not path.exists():
        return {}
    z = np.load(path)
    return dict(zip(z["dates"].tolist(), z["spec"]))


def build_spectra(cfg, stations, from_cache: bool) -> dict[str, dict[str, np.ndarray]]:
    """每站 {date: (N_SEG, K, 3)}。同一資料夾的兄弟組（G6/G7/G8、G2/G3、G23/G24）共用快取，
    只補算快取裡沒有、也不在略過清單中的測站–日。"""
    specs = {st: load_spectra(cfg, st) for st in stations}
    if from_cache:
        return specs
    skipped_path = SPEC_DIR / cfg.gdms_dir.name / "skipped.csv"
    skipped = pd.read_csv(skipped_path, dtype={"date": str}) if skipped_path.exists() else \
        pd.DataFrame(columns=["station", "date", "reason"])
    skip_keys = set(zip(skipped["station"], skipped["date"]))
    refs = [r for st in stations for r in list_day_refs(cfg.gdms_dir, st)
            if r.date_str not in specs[st] and (st, r.date_str) not in skip_keys]
    if not refs:
        return specs
    print(f"[{cfg.group_id}] {len(refs)} 個測站–日要算夜間頻譜", file=sys.stderr)
    new_skips = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for i, fut in enumerate(as_completed([pool.submit(_process_day, r) for r in refs]), 1):
            st, d, res = fut.result()
            if isinstance(res, str):
                new_skips.append({"station": st, "date": d, "reason": res})
            else:
                specs[st][d] = res
            if i % 200 == 0:
                print(f"  {i}/{len(refs)}", file=sys.stderr)
    for st in stations:
        if specs[st]:
            path = _spec_path(cfg, st)
            path.parent.mkdir(parents=True, exist_ok=True)
            ds = sorted(specs[st])
            np.savez(path, dates=np.array(ds), spec=np.stack([specs[st][d] for d in ds]))
    pd.concat([skipped, pd.DataFrame(new_skips, columns=["station", "date", "reason"])]) \
        .to_csv(skipped_path, index=False)
    print(f"[{cfg.group_id}] 新增略過 {len(new_skips)} 個測站–日（缺值／尖峰／純量）", file=sys.stderr)
    return specs


# ---------------------------------------------------------------------------
# 轉換函數與每日能量

def fit_transfer(zn: np.ndarray, bf: np.ndarray) -> np.ndarray:
    """zn：(D, S, K) 近站 Z；bf：(D, S, K, 3) 遠站 X/Y/Z（都只含訓練日）。
    回傳 (K, 3) 複數轉換函數：每個頻率格 k 用 k±SMOOTH_HALF 格的所有子段做最小平方。"""
    K = zn.shape[2]
    t = np.full((K, 3), np.nan + 0j, dtype=np.complex128)
    for k in range(K):
        lo, hi = max(0, k - SMOOTH_HALF), min(K, k + SMOOTH_HALF + 1)
        a = bf[:, :, lo:hi, :].reshape(-1, 3).astype(np.complex128)
        b = zn[:, :, lo:hi].reshape(-1).astype(np.complex128)
        t[k] = np.linalg.lstsq(a, b, rcond=None)[0]
    return t


def band_energy(spec: np.ndarray, band: str) -> np.ndarray:
    """spec：(D, S, K) 複數 → (D,) 頻帶內能量（頻率格加總、子段平均）。"""
    return (np.abs(spec[:, :, BAND_BINS[band]]) ** 2).sum(axis=2).mean(axis=1)


def pair_energies(sn: dict, sf: dict, dates: list[str], train: set[str]) -> tuple[dict, dict]:
    """一對近站／遠站的每日 res／direct 能量（兩站都有資料的日子），與擬合資訊。"""
    both = [d for d in dates if d in sn and d in sf]
    tr = [d for d in both if d in train]
    info = {"n_days_both": len(both), "n_train_days": len(tr)}
    if len(tr) < MIN_TRAIN_DAYS:
        return {}, info
    zn_all = np.stack([sn[d][:, :, 2] for d in both])
    bf_all = np.stack([sf[d] for d in both])
    idx_tr = np.array([d in train for d in both])
    t = fit_transfer(zn_all[idx_tr], bf_all[idx_tr])
    resid = zn_all - np.einsum("dskc,kc->dsk", bf_all.astype(np.complex128), t)
    direct = zn_all - bf_all[..., 2]
    out = {}
    for band in BANDS:
        e_res, e_dir, e_raw = band_energy(resid, band), band_energy(direct, band), band_energy(zn_all, band)
        out[f"{band}_res"] = pd.Series(e_res, index=both)
        out[f"{band}_direct"] = pd.Series(e_dir, index=both)
        # 預測得多好：訓練日（平靜）殘差 ÷ 原始能量的中位數
        info[f"{band}_quiet_res_over_raw"] = round(float(np.median(e_res[idx_tr] / e_raw[idx_tr])), 4)
    return out, info


def _normalize(s: pd.Series, quiet: set[str]) -> pd.Series:
    q = s[s.index.isin(quiet)]
    return s / q.median() if len(q) and q.median() > 0 else s * np.nan


def index_table(specs: dict, near, far, dates: list[str], storm: set[str], exclude: set[str]) -> tuple[pd.DataFrame, dict]:
    """每日指標：res／direct 是配對中位數，raw 是近站中位數；都先除以自己非磁暴日中位數。"""
    quiet = {d for d in dates if d not in storm}
    train = quiet - exclude
    out = pd.DataFrame(index=pd.Index(dates, name="date"))
    out["is_storm"] = [d in storm for d in dates]
    fits: dict = {"n_train_candidates": len(train)}
    pair_series: dict[str, dict[str, pd.Series]] = {f"{b}_{v}": {} for b in BANDS for v in ("res", "direct")}
    for n in near:
        for f in far:
            if n not in specs or f not in specs:
                continue
            ser, info = pair_energies(specs[n], specs[f], dates, train)
            fits[f"{n}-{f}"] = info
            for key, s in ser.items():
                pair_series[key][f"{n}-{f}"] = _normalize(s, quiet)
    for key, d in pair_series.items():
        wide = pd.DataFrame(d).reindex(dates) if d else pd.DataFrame(index=dates)
        out[key] = wide.median(axis=1) if d else np.nan
        out[f"{key.split('_')[0]}_n_pairs"] = wide.notna().sum(axis=1) if d else 0
    for band in BANDS:
        raw = {}
        for n in near:
            sn = specs.get(n, {})
            ds = [d for d in dates if d in sn]
            if ds:
                e = band_energy(np.stack([sn[d][:, :, 2] for d in ds]), band)
                raw[n] = _normalize(pd.Series(e, index=ds), quiet)
        out[f"{band}_raw"] = pd.DataFrame(raw).reindex(dates).median(axis=1) if raw else np.nan
    return out.reset_index(), fits


def tested_blocks(cfg) -> tuple[set[str], set[str]]:
    """排名檢定會檢定的兩個區塊：錨點震前窗口、緊接在前的安慰劑區塊（同 surrogate_test）。"""
    end = pre_event_anchor(cfg.anchor_event)
    window = set(pre_event_window(end, mc.PRE_WINDOW_DAYS))
    prev = set(pre_event_window(end - pd.Timedelta(days=mc.PRE_WINDOW_DAYS), mc.PRE_WINDOW_DAYS))
    return window, prev


def run_group(group_id: str, from_cache: bool) -> pd.DataFrame | None:
    cfg = load_group_config(group_id)
    if not cfg.xyz_pool.sufficient:
        print(f"[{group_id}] XYZ 測站池不足，略過", file=sys.stderr)
        return None
    near, far = cfg.xyz_pool.near, cfg.xyz_pool.far
    specs = build_spectra(cfg, list(near) + list(far), from_cache)
    storm = set(pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})["date"])
    dates = sorted({d for st in specs.values() for d in st})
    window, prev = tested_blocks(cfg)
    out, fits = index_table(specs, near, far, dates, storm, window | prev)
    out.to_csv(cfg.interim_dir / "ulf_istf_index.csv", index=False)
    (cfg.interim_dir / "ulf_istf_fits.json").write_text(json.dumps(fits, indent=2, ensure_ascii=False))
    print(f"[{group_id}] 已寫入 ulf_istf_index.csv（{len(out)} 列）", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------
# 診斷

def _mad(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return 1.4826 * float(np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


def daily_kp_max(cfg) -> pd.Series:
    path = cfg.external_dir / "kp.csv"
    if not path.exists():  # 共用資料夾的兄弟組（例如 G24 用 G23 的）：Kp 是全球指數，只看日期
        sibs = [g for g in ULF_GROUPS if g != cfg.group_id
                and (PROJECT_DIR / "data" / "external" / g / "kp.csv").exists()
                and load_group_config(g).gdms_dir == cfg.gdms_dir]
        if not sibs:
            return pd.Series(dtype=float)
        path = PROJECT_DIR / "data" / "external" / sibs[0] / "kp.csv"
    kp = pd.read_csv(path, dtype={"date": str})
    return kp.groupby("date")["kp"].max()


def series_diagnostics(s: pd.Series, storm: set[str], window: set[str], prev: set[str], kp: pd.Series,
                       log: bool) -> dict:
    s = s.dropna()
    if log:
        s = np.log10(s[s > 0])
    is_storm = s.index.isin(list(storm))
    out = {"n_days": int(len(s)), "n_storm_days": int(is_storm.sum()),
           "storm_fraction_all": round(float(is_storm.mean()), 3) if len(s) else None}
    v = s.to_numpy()
    if is_storm.sum() >= 5 and (~is_storm).sum() >= 5:
        mq = _mad(v[~is_storm])
        out["mad_storm_over_quiet"] = round(_mad(v[is_storm]) / mq, 3)
        out["storm_median_shift_mad"] = round((np.median(v[is_storm]) - np.median(v[~is_storm])) / mq, 3)
    k = kp.reindex(s.index)
    ok = k.notna().to_numpy()
    if ok.sum() >= 20:
        out["spearman_kp"] = round(float(spearmanr(v[ok], k.to_numpy()[ok]).statistic), 3)
    in_win = s[s.index.isin(list(window))]
    out["window_days_storm_excluded"] = int((~in_win.index.isin(list(storm))).sum())
    out["window_days_storm_kept"] = int(len(in_win))
    in_prev = s[s.index.isin(list(prev))]
    out["prev_block_days_storm_excluded"] = int((~in_prev.index.isin(list(storm))).sum())
    out["prev_block_days_storm_kept"] = int(len(in_prev))
    if len(in_win):
        out["window_storm_fraction"] = round(float(in_win.index.isin(list(storm)).mean()), 3)
        out["window_max_date"], out["window_min_date"] = str(in_win.idxmax()), str(in_win.idxmin())
        out["window_max_is_storm"] = out["window_max_date"] in storm
        out["window_min_is_storm"] = out["window_min_date"] in storm
    return out


def group_diagnostics(group_id: str) -> dict | None:
    cfg = load_group_config(group_id)
    path = cfg.interim_dir / "ulf_istf_index.csv"
    if not path.exists():
        return None
    idx = pd.read_csv(path, dtype={"date": str}).set_index("date")
    storm = set(pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})["date"])
    window, prev = tested_blocks(cfg)
    kp = daily_kp_max(cfg)
    zh_path = cfg.interim_dir / "ulf_near_far_index.csv"
    zh = pd.read_csv(zh_path, dtype={"date": str}).set_index("date") if zh_path.exists() else None
    fits = json.loads((cfg.interim_dir / "ulf_istf_fits.json").read_text())
    res = {}
    for band in BANDS:
        r = {v: series_diagnostics(idx[f"{band}_{v}"], storm, window, prev, kp, log=True) for v in VERSIONS}
        if zh is not None and f"{band}_diff_zh" in zh:
            r["zh_diff"] = series_diagnostics(zh[f"{band}_diff_zh"], storm, window, prev, kp, log=False)
        ratios = [v[f"{band}_quiet_res_over_raw"] for k, v in fits.items()
                  if isinstance(v, dict) and f"{band}_quiet_res_over_raw" in v]
        r["quiet_res_over_raw_median"] = round(float(np.median(ratios)), 4) if ratios else None
        r["n_pairs_fitted"] = len(ratios)
        res[band] = r
    return res


def summarize(groups: list[str]) -> dict:
    per_group = {g: d for g in groups if (d := group_diagnostics(g)) is not None}
    return {"question": "在原始資料層用遠站以站間轉換函數預測近站 Z、取殘差，磁暴日能不能放回分析？",
            "status": "探索性診斷（2026-11-13）；只列描述數字，不設通過門檻、不算前兆 p 值",
            "reference": "Harada, Hattori & Isezaki 2004, doi:10.1016/j.pce.2004.03.002（概念；估計細節為本專案自訂）",
            "per_group": per_group}


VERSION_LABEL = {"res": "ISTF 殘差", "direct": "直接相減", "raw": "不相減", "zh_diff": "現行 Z/H 近減遠"}


def render(s: dict) -> str:
    L = ["# 站間轉換函數差分：磁暴日能不能放回來？（探索性診斷）", "",
         "用遠站 X/Y/Z 在頻域（1200 秒 FFT 子段，頻率格 ±2 合併）以最小平方預測近站 Z，取殘差的夜間能量。"
         "轉換函數只用非磁暴日、且不在錨點震前 30 天窗口與安慰劑前一區塊內的日子估計。"
         "概念取自 Harada, Hattori & Isezaki (2004)，論文用小波；這裡的估計細節是本專案自訂。", "",
         "**只列描述數字**：不設通過門檻、不算前兆 p 值；能不能用由使用者看完再定。", "",
         "四個版本：ISTF 殘差｜直接相減（Z 近 − Z 遠）｜不相減（近站 Z）｜現行 Z/H 近減遠。"
         "前三個是能量，取 log10 後計算；Z/H 照 11-02 用原尺度。", "",
         "- **MAD 比** = 磁暴日 MAD ÷ 平靜日 MAD；約 1 表示磁暴日和平靜日一樣分散。",
         "- **中位偏移** = (磁暴日中位數 − 平靜日中位數) ÷ 平靜日 MAD；約 0 表示磁暴日沒有整體偏高或偏低。",
         "- **ρ(Kp)** = 指標與當日最大 Kp 的 Spearman 相關；約 0 表示指標不跟著地磁活動走。",
         "- **平靜日殘差占比** = 訓練日殘差能量 ÷ 原始能量的中位數；越小表示遠站預測得越好。", ""]
    for band in BANDS:
        L += [f"## {band.upper()}", "",
              "| 組 | 平靜日殘差占比 | 版本 | 磁暴日／總日數 | MAD 比 | 中位偏移 | ρ(Kp) | 窗口磁暴比例／全序列 | 窗口最高日／最低日是磁暴？ |",
              "|---|---|---|---|---|---|---|---|---|"]
        for g, r in s["per_group"].items():
            b = r[band]
            for i, v in enumerate(("res", "direct", "raw", "zh_diff")):
                d = b.get(v)
                if d is None:
                    continue
                head = f"| {g} | {b['quiet_res_over_raw_median']} " if i == 0 else "| | "
                L.append(head + f"| {VERSION_LABEL[v]} | {d['n_storm_days']}／{d['n_days']} "
                         f"| {d.get('mad_storm_over_quiet')} | {d.get('storm_median_shift_mad')} | {d.get('spearman_kp')} "
                         f"| {d.get('window_storm_fraction')}／{d.get('storm_fraction_all')} "
                         f"| {d.get('window_max_is_storm')}／{d.get('window_min_is_storm')} |")
        L.append("")
    L += ["## 跨組彙整（中位數，各組 × 頻帶）", "",
          "| 頻帶 | 版本 | MAD 比 | 中位偏移 | ρ(Kp) | MAD 比 ≤ 1.2 的組數* |", "|---|---|---|---|---|---|"]
    for band in BANDS:
        for v in ("res", "direct", "raw", "zh_diff"):
            rows = [r[band][v] for r in s["per_group"].values() if v in r[band]]
            def med(key):
                xs = [x[key] for x in rows if x.get(key) is not None]
                return f"{np.median(xs):+.3f}" if xs else "—"
            mads = [x["mad_storm_over_quiet"] for x in rows if x.get("mad_storm_over_quiet") is not None]
            L.append(f"| {band} | {VERSION_LABEL[v]} | {med('mad_storm_over_quiet')} | {med('storm_median_shift_mad')} "
                     f"| {med('spearman_kp')} | {sum(m <= 1.2 for m in mads)}／{len(mads)} |")
    L += ["", "\\* 1.2 只是方便閱讀的參考線，不是判準。", "",
          "## 樣本增益（ISTF 殘差版本）", "",
          "震前 30 天窗口與安慰劑前一區塊的可用天數：排除磁暴日 → 保留磁暴日。"
          f"排名檢定要求區塊內至少 {mc.MIN_WINDOW_DAYS} 個可用日。", "",
          "| 組 | 頻帶 | 震前窗口 | 安慰劑區塊 | 整條序列 |", "|---|---|---|---|---|"]
    for g, r in s["per_group"].items():
        for band in BANDS:
            d = r[band]["res"]
            L.append(f"| {g} | {band} | {d['window_days_storm_excluded']} → {d['window_days_storm_kept']} "
                     f"| {d['prev_block_days_storm_excluded']} → {d['prev_block_days_storm_kept']} "
                     f"| {d['n_days'] - d['n_storm_days']} → {d['n_days']} |")
    L += ["", "## 和文獻的差異", "",
          "- Harada et al. (2004) 用小波轉換估計站間轉換函數；這裡用 FFT 分段＋頻率格合併的最小平方（自訂）。",
          "- 文獻用一個遠端參考站（柿岡）；這裡每組多個近站 × 遠站，取配對中位數，站距也短得多。",
          "- 轉換函數在訓練日（平靜日）上是樣本內估計；每個頻率格 3 個參數、數千個子段，過擬合可忽略。", ""]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# 自我測試（合成資料）

def synthetic_specs(n_days: int = 160, storm_days=range(30, 70), local_days=range(140, 160), seed: int = 0,
                    drop_n2_from: int | None = None):
    """共同源場（Pc3/Pc4 頻帶雜訊）× 各站不同的頻率相依轉換（1–2 秒延遲）；源場和近站局部
    雜訊每天各有對數常態起伏（像真實資料），磁暴日源場再 ×10。local_days 加只有 n1 才有的
    20 秒正弦。局部雜訊約占近站 Z 能量的一成；若局部雜訊極小（< 1%）而磁暴極大，
    FFT 分段對延遲的近似誤差（能量約 1e-4）也會被磁暴放大到看得見。"""
    from scipy.signal import butter, sosfiltfilt
    rng = np.random.default_rng(seed)
    sos = butter(4, [BANDS["pc4"][0] / 0.5, BANDS["pc3"][1] / 0.5], btype="band", output="sos")
    dates = pd.date_range("2024-01-01", periods=n_days).strftime("%Y%m%d").tolist()
    storm = {dates[i] for i in storm_days}
    specs = {st: {} for st in ("n1", "n2", "f1", "f2")}
    t = np.arange(N_NIGHT)
    for i, d in enumerate(dates):
        amp = np.exp(rng.normal(0, 0.5)) * (10.0 if d in storm else 1.0)
        loc_amp = 0.2 * np.exp(rng.normal(0, 0.2))
        sx, sy = (amp * sosfiltfilt(sos, rng.normal(size=N_NIGHT)) for _ in range(2))
        noise = lambda a: a * sosfiltfilt(sos, rng.normal(size=N_NIGHT))  # noqa: E731
        stations = {
            # 近站 Z：對源場有延遲（頻率相依相位）與不同增益；遠站 Z 的響應完全不同
            "n1": (1.1 * sx, 0.9 * sy, 0.6 * np.roll(sx, 2) - 0.3 * np.roll(sy, 1) + noise(loc_amp)),
            "n2": (1.0 * sx, 1.0 * sy, 0.4 * np.roll(sx, 1) + 0.5 * np.roll(sy, 2) + noise(loc_amp)),
            "f1": (sx + noise(0.01), sy + noise(0.01), 0.15 * sx + noise(0.01)),
            "f2": (0.9 * sx + noise(0.01), 1.1 * sy + noise(0.01), -0.1 * sy + noise(0.01)),
        }
        if i in local_days:
            x, y, z = stations["n1"]
            stations["n1"] = (x, y, z + 1.0 * np.sin(2 * np.pi * t / 20))
        for st, (x, y, z) in stations.items():
            if st == "n2" and drop_n2_from is not None and i >= drop_n2_from:
                continue
            specs[st][d] = segment_spectra(np.column_stack([x, y, z]))
    return specs, dates, storm


def self_test() -> bool:
    ok = True
    t = np.arange(N_NIGHT, dtype=float)
    blk = np.column_stack([np.sin(2 * np.pi * t / 20), np.sin(2 * np.pi * t / 100), np.sin(2 * np.pi * t / 5)])
    sp = segment_spectra(blk)[None]
    e = {c: {b: float(band_energy(sp[..., i], b)[0]) for b in BANDS} for i, c in enumerate("XYZ")}
    c1 = e["X"]["pc3"] > 100 * e["X"]["pc4"] and e["Y"]["pc4"] > 100 * e["Y"]["pc3"] \
        and e["Z"]["pc3"] < 1e-3 * e["X"]["pc3"]
    print(f"[self-test] 頻帶：20 秒落在 Pc3、100 秒落在 Pc4、5 秒被濾掉  {'PASS' if c1 else 'FAIL'}", file=sys.stderr)
    ok &= c1

    specs, dates, storm = synthetic_specs()
    exclude = set(dates[120:140])
    tab, _ = index_table(specs, ["n1", "n2"], ["f1", "f2"], dates, storm, exclude)
    tab = tab.set_index("date")
    calm = [d for d in dates[:140] if d not in storm]
    ratio = {v: tab["pc3_" + v][list(storm)].median() / tab["pc3_" + v][calm].median() for v in VERSIONS}
    shift = {}
    for v in VERSIONS:
        s = np.log10(tab["pc3_" + v])
        shift[v] = (s[list(storm)].median() - s[calm].median()) / _mad(s[calm].to_numpy())
    c2 = ratio["res"] < 1.5 and ratio["direct"] > 20 and ratio["raw"] > 20 and abs(shift["res"]) < 1
    print(f"[self-test] 磁暴源場 ×10：磁暴／平靜中位數比 ISTF {ratio['res']:.2f}、直接相減 {ratio['direct']:.0f}、"
          f"不相減 {ratio['raw']:.0f}；ISTF 中位偏移 {shift['res']:+.2f} MAD  {'PASS' if c2 else 'FAIL'}", file=sys.stderr)
    ok &= c2
    loc = tab["pc3_res"].iloc[140:].median() / tab["pc3_res"].iloc[100:140].median()
    c3 = loc > 2
    print(f"[self-test] 只有近站有的局部 20 秒訊號：殘差上升 {loc:.1f} 倍  {'PASS' if c3 else 'FAIL'}", file=sys.stderr)
    ok &= c3

    specs2, dates2, storm2 = synthetic_specs(drop_n2_from=100, local_days=range(0))
    tab2, _ = index_table(specs2, ["n1", "n2"], ["f1", "f2"], dates2, storm2, set())
    p = tab2.set_index("date")["pc3_res"]
    jump = p.iloc[100:120].median() / p.iloc[80:100].median()
    c4 = abs(np.log(jump)) < np.log(1.5)
    print(f"[self-test] n2 中途離開前後殘差比 {jump:.2f}（≈1）  {'PASS' if c4 else 'FAIL'}", file=sys.stderr)
    ok &= c4
    return bool(ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", action="append")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--from-cache", action="store_true", help="只用既有頻譜快取，不讀原始資料")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if not self_test():
        sys.exit("自我測試失敗")
    if args.self_test and not (args.group or args.all):
        return
    groups = list(ULF_GROUPS) if args.all else (args.group or [])
    for g in groups:
        run_group(g, args.from_cache)
    s = summarize(list(ULF_GROUPS))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "summary.json").write_text(json.dumps(s, indent=2, ensure_ascii=False))
    (OUT_DIR / "summary.md").write_text(render(s))
    print(f"已寫入 {OUT_DIR / 'summary.md'}（{len(s['per_group'])} 組）", file=sys.stderr)


if __name__ == "__main__":
    main()
