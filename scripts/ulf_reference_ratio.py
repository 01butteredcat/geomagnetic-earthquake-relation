"""遠站參考比值（Han et al. 2014 的方法，照 Warden et al. 2020 的獨立重做寫法）。

回答「既然近站減遠站，能不能一開始就用遠站把磁暴背景消掉」：不看 Z/H，而是用
遠站的每日 ULF 能量預測近站的能量，取「實測 ÷ 預測」比值 P。全球性的變化
（包括磁暴）兩站都有，會被預測掉；P 偏高代表近站有遠站沒有的局部能量。

文獻：
  Han, P., Hattori, K., Hirokawa, M., Zhuang, J., et al. (2014). Statistical analysis of
    ULF seismomagnetic phenomena at Kakioka, Japan, during 2001–2010. JGR Space Physics
    119(6), 4998–5011. doi:10.1002/2014JA019789
  Warden, S., MacLean, L., Lemon, J., & Schneider, D. (2020). Statistical analysis of
    pre-earthquake electromagnetic anomalies in the ULF range. JGR Space Physics.
    doi:10.1029/2020JA027955（獨立重做 Han 2014，方法細節以這篇為準）

每站每日（照文獻）：
  1. Z 分量，1 Hz；取 130 分鐘（UTC 18:10–20:20），段內有缺值或 INTERMAGNET 尖峰就略過。
  2. db5 小波 DWT 6 層，只留第 6 層細節重建（中心頻率約 0.01 Hz，週期約 100 秒）。
  3. 截中間 90 分鐘（UTC 18:30–20:00 = 當地 02:30–04:00），能量 = 變異數。
  H（X、Y 能量相加）一併算，只供比對；文獻主分量是 Z。

比值（照文獻，最小平方線性迴歸，直接對能量、不取對數）：
  每一對（近站 n、遠站 f）：E_n* = beta * E_f + c，P = E_n / E_n*（E_n* <= 0 時為 NaN）。
  文獻只有一對測站；這裡每組有多個近站、遠站，每日取所有可用配對 P 的中位數。
  每一對的 P 都以 1 為中心，所以測站加入／離開不會讓指標跳動。

三個版本（回答磁暴問題）：
  all        —— 磁暴日留在迴歸和序列中（文獻做法；Han 只剔除少數高能量日，
                Warden 發現那些大多是磁暴日，但沒有系統性排除）← 主版本
  fit_quiet  —— 只用平靜日估 beta、c，序列保留所有日子
  drop_storm —— 同 fit_quiet，序列也去掉磁暴／恢復期日（和現行主要檢定一致）

統計（主要設定在看結果前寫死：PRIMARY）：
  - 和主要檢定相同的排名檢定（不重疊 30 天區塊、Fisher 合併、同一套敏感度），
    方向是上尾（文獻報告震前能量比值增加）。
  - 描述性：P > 中位數 + 1.5 IQR 的異常日（文獻門檻，用整條序列算），窗口內 vs 窗口外的比例。
  - 磁暴日 vs 平靜日的 P 的 MAD 比：約 1 代表遠站確實把磁暴消掉了。

輸出：
  data/interim/<group>/ulf_ref_daily.csv        每站每日能量
  data/interim/<group>/ulf_ref_ratio_index.csv  每日 P（三個版本 × Z/H）
  data/interim/ulf_ref_ratio/summary.{json,md}  跨組結果

用法：
  ulf_reference_ratio.py --self-test
  ulf_reference_ratio.py --group G10 [--from-daily]
  ulf_reference_ratio.py --all [--from-daily]
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pywt

sys.path.insert(0, str(Path(__file__).parent))
import method_comparison as mc  # noqa: E402
from build_daily_features import _despike  # noqa: E402
from common import PROJECT_DIR, list_day_refs, load_group_config  # noqa: E402
from cross_group_analysis import pre_event_anchor, pre_event_window  # noqa: E402
from parser import parse_day_file  # noqa: E402
from surrogate_test import FAMILY_DROP, SEED, ULF_GROUPS, rank_sensitivity  # noqa: E402

WAVELET = "db5"
DWT_LEVEL = 6
# 130 分鐘做 DWT（避開邊緣效應），截中間 90 分鐘：Warden et al. 2020 第 2.2.2 節
DWT_START = pd.Timedelta(hours=18, minutes=10)   # UTC
CROP_START = pd.Timedelta(hours=18, minutes=30)  # UTC 18:30 = 當地 02:30
CROP_END = pd.Timedelta(hours=20)                # UTC 20:00 = 當地 04:00（不含）
DWT_END = pd.Timedelta(hours=20, minutes=20)
N_DWT = int((DWT_END - DWT_START).total_seconds())     # 7800
CROP = slice(int((CROP_START - DWT_START).total_seconds()), int((CROP_END - DWT_START).total_seconds()))
N_WORKERS = 10
VERSIONS = ("all", "fit_quiet", "drop_storm")
COMPONENTS = ("Z", "H")
# 在看任何結果之前就定下來（2026-11-02）：文獻版本、文獻分量、文獻方向
PRIMARY = {"version": "all", "component": "Z", "tail": "upper"}
IQR_K = 1.5  # 中位數 + 1.5 IQR（Han 2014，引用 Zwillinger & Kokoska 2000）


# ---------------------------------------------------------------------------
# 每站每日能量

def dwt_level_detail(x: np.ndarray, wavelet: str = WAVELET, level: int = DWT_LEVEL) -> np.ndarray:
    """只用第 `level` 層細節係數重建訊號（其他係數歸零）。"""
    coeffs = pywt.wavedec(x, wavelet, level=level)
    kept = [np.zeros_like(c) for c in coeffs]
    kept[1] = coeffs[1]  # coeffs = [cA_L, cD_L, cD_{L-1}, ..., cD_1]
    return pywt.waverec(kept, wavelet)[: len(x)]


def band_energy(x: np.ndarray) -> float:
    return float(np.var(dwt_level_detail(x)[CROP]))


def _process_day(ref):
    station, date_str = ref.station, ref.date_str
    try:
        df = parse_day_file(ref, station)
    except Exception as exc:  # noqa: BLE001
        return {"station": station, "date": date_str, "error": str(exc)}
    if "Z" not in df.columns:
        return {"station": station, "date": date_str, "error": "scalar_only"}
    day = pd.Timestamp(date_str)
    seg = df.loc[(df.index >= day + DWT_START) & (df.index < day + DWT_END)]
    if len(seg) != N_DWT:
        return {"station": station, "date": date_str, "error": f"window_has_{len(seg)}_samples"}
    row = {"station": station, "date": date_str}
    for comp, cols in (("Z", ("Z",)), ("H", ("X", "Y"))):
        energy = 0.0
        for c in cols:
            s = seg[c].astype(float)
            if s.isna().any() or _despike(s)[1].any():
                energy = np.nan  # 缺值或尖峰：DWT 會把它攤到前後，整段不用
                break
            energy += band_energy(s.to_numpy())
        row[f"{comp}_energy"] = energy
    if np.isnan(row["Z_energy"]) and np.isnan(row["H_energy"]):
        return {"station": station, "date": date_str, "error": "gap_or_spike_in_window"}
    return row


def build_daily_energy(cfg, stations) -> pd.DataFrame:
    refs = [r for st in stations for r in list_day_refs(cfg.gdms_dir, st)]
    print(f"[{cfg.group_id}] {len(refs)} 個測站–日要算 DWT 能量", file=sys.stderr)
    rows, errors = [], []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for fut in as_completed([pool.submit(_process_day, r) for r in refs]):
            r = fut.result()
            (errors if "error" in r else rows).append(r)
    daily = pd.DataFrame(rows, columns=["station", "date", "Z_energy", "H_energy"]).sort_values(["station", "date"])
    daily.to_csv(cfg.interim_dir / "ulf_ref_daily.csv", index=False)
    if errors:
        pd.DataFrame(errors).to_csv(cfg.interim_dir / "ulf_ref_errors.csv", index=False)
    print(f"[{cfg.group_id}] 已寫入 ulf_ref_daily.csv（{len(daily)} 列，略過 {len(errors)} 個）", file=sys.stderr)
    return daily


# ---------------------------------------------------------------------------
# 比值

def pair_ratio(e_near: pd.Series, e_far: pd.Series, fit_mask: pd.Series) -> tuple[pd.Series, dict]:
    """最小平方 E_near = beta * E_far + c（只用 fit_mask 的日子估），回傳 P = E_near / 預測值。"""
    both = e_near.notna() & e_far.notna()
    fit = both & fit_mask
    if fit.sum() < 10:
        return pd.Series(np.nan, index=e_near.index), {"n_fit": int(fit.sum())}
    beta, c = np.polyfit(e_far[fit].to_numpy(), e_near[fit].to_numpy(), 1)
    pred = beta * e_far + c
    p = (e_near / pred).where(both & (pred > 0))
    r = float(np.corrcoef(e_far[fit], e_near[fit])[0, 1])
    return p, {"beta": float(beta), "c": float(c), "pearson_r": round(r, 4), "n_fit": int(fit.sum()),
               "n_pred_le_0": int((both & (pred <= 0)).sum())}


def ratio_table(daily: pd.DataFrame, near, far, storm_dates: set[str]) -> tuple[pd.DataFrame, dict]:
    dates = sorted(daily["date"].astype(str).unique())
    is_storm = pd.Series([d in storm_dates for d in dates], index=dates)
    out = pd.DataFrame(index=pd.Index(dates, name="date"))
    out["is_storm"] = is_storm
    fits: dict = {}
    for comp in COMPONENTS:
        wide = daily.pivot_table(index="date", columns="station", values=f"{comp}_energy").reindex(dates)
        for version in VERSIONS:
            fit_mask = pd.Series(True, index=dates) if version == "all" else ~is_storm
            ps = {}
            for n in near:
                for f in far:
                    if n in wide and f in wide:
                        ps[f"{n}-{f}"], fits[f"{comp}/{version}/{n}-{f}"] = pair_ratio(wide[n], wide[f], fit_mask)
            p = pd.DataFrame(ps).median(axis=1) if ps else pd.Series(np.nan, index=dates)
            if version == "drop_storm":
                p = p.where(~is_storm)
            out[f"{comp}_P_{version}"] = p
            out[f"{comp}_n_pairs"] = pd.DataFrame(ps).notna().sum(axis=1) if ps else 0
    return out.reset_index(), fits


def build_ratio_index(cfg, daily: pd.DataFrame) -> pd.DataFrame:
    storm = set(pd.read_csv(cfg.interim_dir / "storm_days.csv", dtype={"date": str})["date"])
    out, fits = ratio_table(daily.assign(date=daily["date"].astype(str)), cfg.xyz_pool.near, cfg.xyz_pool.far, storm)
    out.to_csv(cfg.interim_dir / "ulf_ref_ratio_index.csv", index=False)
    (cfg.interim_dir / "ulf_ref_ratio_fits.json").write_text(json.dumps(fits, indent=2))
    print(f"[{cfg.group_id}] 已寫入 ulf_ref_ratio_index.csv（{len(out)} 列）", file=sys.stderr)
    return out


# ---------------------------------------------------------------------------
# 統計

def _mad(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return 1.4826 * float(np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


def group_stats(cfg, idx: pd.DataFrame) -> dict:
    window_end = pre_event_anchor(cfg.anchor_event)
    window = set(pre_event_window(window_end, mc.PRE_WINDOW_DAYS))
    prev_window = set(pre_event_window(window_end - pd.Timedelta(days=mc.PRE_WINDOW_DAYS), mc.PRE_WINDOW_DAYS))
    idx = idx.assign(date=idx["date"].astype(str)).set_index("date")
    res = {}
    for comp in COMPONENTS:
        for version in VERSIONS:
            s = idx[f"{comp}_P_{version}"].dropna()
            key = f"{comp}/{version}"
            if len(s) < 20:
                res[key] = {"error": f"只有 {len(s)} 天"}
                continue
            storm = idx.loc[s.index, "is_storm"].to_numpy(dtype=bool)
            q1, med, q3 = np.percentile(s, [25, 50, 75])
            thr = med + IQR_K * (q3 - q1)
            in_win = s.index.isin(window)
            anom = s.to_numpy() > thr
            rk = mc.rank_window_test(s, window, tail="upper")
            prev = mc.rank_window_test(s, prev_window, tail="upper")
            res[key] = {
                **{k: rk.get(k) for k in ("rank_p", "rank_n_fake_windows", "null_ps")},
                "placebo_prev_block": {k: prev.get(k) for k in ("rank_p", "rank_n_fake_windows", "null_ps")},
                "n_days": int(len(s)), "n_storm_days": int(storm.sum()),
                "mad_storm_over_quiet": round(_mad(s.to_numpy()[storm]) / _mad(s.to_numpy()[~storm]), 3)
                if storm.sum() >= 5 and (~storm).sum() >= 5 else None,
                "anomaly_threshold": round(float(thr), 4),
                "anomaly_rate_window": round(float(anom[in_win].mean()), 4) if in_win.any() else None,
                "anomaly_rate_outside": round(float(anom[~in_win].mean()), 4),
                "n_window_days": int(in_win.sum()),
                "window_max_date": str(s[in_win].idxmax()) if in_win.any() else None,
                "window_max_is_storm": bool(idx.loc[s[in_win].idxmax(), "is_storm"]) if in_win.any() else None,
            }
    return res


def combine(per_group: dict[str, dict], key: str) -> dict:
    rows = [(g, r[key]) for g, r in per_group.items() if r.get(key, {}).get("rank_p") is not None]
    rng = np.random.default_rng([SEED, sum(map(ord, key))])
    out = {}
    for name, keep in (("all_groups", lambda g: True), ("one_per_family", lambda g: g not in FAMILY_DROP)):
        sel = [(g, pe) for g, pe in rows if keep(g)]
        out[name] = {**mc.fisher_across_groups([pe["rank_p"] for _, pe in sel], [pe["null_ps"] for _, pe in sel], rng),
                     "groups": [g for g, _ in sel]} if sel else None
    out["sensitivity"] = rank_sensitivity(rows, np.random.default_rng([SEED, sum(map(ord, key)), 1])) if rows else None
    return out


def summarize(groups: list[str]) -> dict:
    per_group = {}
    for g in groups:
        cfg = load_group_config(g)
        path = cfg.interim_dir / "ulf_ref_ratio_index.csv"
        if path.exists():
            per_group[g] = group_stats(cfg, pd.read_csv(path, dtype={"date": str}))
    keys = [f"{c}/{v}" for c in COMPONENTS for v in VERSIONS]
    return {"primary": PRIMARY, "references": ["Han et al. 2014, doi:10.1002/2014JA019789",
                                               "Warden et al. 2020, doi:10.1029/2020JA027955"],
            "status": "探索性：使用者 2026-11-02 提問後新增；主要檢定仍是 Pc3 Z/H 下尾",
            "combined": {k: combine(per_group, k) for k in keys},
            "per_group": {g: {k: {kk: vv for kk, vv in v.items() if kk != "null_ps"} for k, v in r.items()}
                          for g, r in per_group.items()}}


def render(s: dict) -> str:
    def cell(c):
        return "—" if not c else f"{c['p']}（{c['n_groups']} 組）"
    L = ["# 遠站參考比值（Han 2014／Warden 2020 方法）", "",
         "用遠站的每日 ULF 能量（db5 DWT 第 6 層，約 0.01 Hz，當地 02:30–04:00）以最小平方預測近站，"
         "比值 P = 實測 ÷ 預測。P 高 = 近站有遠站沒有的局部能量。", "",
         f"**主要（看結果前定下）**：版本 `{PRIMARY['version']}`、分量 {PRIMARY['component']}、{PRIMARY['tail']} 尾。"
         "探索性；現行主要檢定仍是 Pc3 Z/H 下尾。", "",
         "## 跨組合併（排名檢定、上尾、Fisher）", "",
         "| 分量 | 版本 | 全部組 | 每家族一組 | 不含 G10 | 前一區段安慰劑 |", "|---|---|---|---|---|---|"]
    for key, c in s["combined"].items():
        comp, ver = key.split("/")
        sv = c.get("sensitivity") or {}
        L.append(f"| {comp} | {ver} | {cell(c.get('all_groups'))} | {cell(c.get('one_per_family'))} "
                 f"| {(sv.get('drop_discovery_G10') or {}).get('p')} | {(sv.get('placebo_prev_block') or {}).get('p')} |")
    L += ["", "## 磁暴有沒有被遠站消掉？（Z 分量）", "",
          "磁暴／平靜 MAD 比 = 磁暴日 P 的 MAD ÷ 平靜日 P 的 MAD；約 1 = 遠站把磁暴消掉了。", "",
          "| 組 | 磁暴日／總日數 | MAD 比（all） | MAD 比（fit_quiet） | 窗口最高日是磁暴日？ | 排名 p：all／fit_quiet／drop_storm "
          "| 異常日比例：窗口內／窗口外（all） |", "|---|---|---|---|---|---|---|"]
    for g, r in s["per_group"].items():
        a, fq, ds = r.get("Z/all", {}), r.get("Z/fit_quiet", {}), r.get("Z/drop_storm", {})
        L.append(f"| {g} | {a.get('n_storm_days')}／{a.get('n_days')} | {a.get('mad_storm_over_quiet')} "
                 f"| {fq.get('mad_storm_over_quiet')} | {a.get('window_max_is_storm')}（{a.get('window_max_date')}） "
                 f"| {a.get('rank_p')}／{fq.get('rank_p')}／{ds.get('rank_p')} "
                 f"| {a.get('anomaly_rate_window')}／{a.get('anomaly_rate_outside')} |")
    L += ["", "## 和文獻的差異", "",
          "- 文獻只有一對測站（柿岡／鹿屋，相距 > 1000 km）；這裡每組多個近站、遠站，取配對 P 的中位數，遠站距離也短得多。",
          "- Han 2014 剔除了少數高能量日但沒寫明規則；這裡 `all` 版完全不剔除，`drop_storm` 版剔除所有磁暴／恢復期日。",
          "- 統計檢定沿用本專案的排名檢定，不是文獻的 SEA 計數；文獻門檻（中位數 + 1.5 IQR）只做描述。", ""]
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------

def self_test() -> bool:
    ok = True
    t = np.arange(N_DWT, dtype=float)
    e100 = band_energy(np.sin(2 * np.pi * t / 100))
    e10 = band_energy(np.sin(2 * np.pi * t / 10))
    c1 = e100 > 20 * e10
    print(f"[self-test] DWT 第 6 層：100 秒正弦能量 {e100:.3g}，10 秒 {e10:.3g}  {'PASS' if c1 else 'FAIL'}",
          file=sys.stderr)
    ok &= c1

    # 共同源場 × 各站增益；磁暴日源場放大 20 倍；n2 中途離開；n1 最後 15 天有局部增強
    rng = np.random.default_rng(0)
    dates = pd.date_range("2024-01-01", periods=150).strftime("%Y%m%d")
    src = np.exp(rng.normal(0, 0.5, len(dates)))
    storm = set(dates[40:55]) | set(dates[90:100])
    src = src * np.where(np.isin(dates, list(storm)), 20.0, 1.0)
    gains = {"n1": 2.0, "n2": 5.0, "f1": 1.0, "f2": 0.7}
    rows = []
    for i, d in enumerate(dates):
        for st, gn in gains.items():
            if st == "n2" and i >= 70:
                continue
            e = gn * src[i] * np.exp(rng.normal(0, 0.1))
            if st == "n1" and i >= 135:
                e *= 3.0
            rows.append({"station": st, "date": d, "Z_energy": e, "H_energy": e})
    tab, _ = ratio_table(pd.DataFrame(rows), ["n1", "n2"], ["f1", "f2"], storm)
    tab = tab.set_index("date")
    p = tab["Z_P_all"]
    st_mask = tab["is_storm"].to_numpy(dtype=bool)
    mad_ratio = _mad(p[st_mask].to_numpy()) / _mad(p[~st_mask & (np.arange(len(p)) < 135)].to_numpy())
    c2 = mad_ratio < 2.0
    print(f"[self-test] 磁暴日源場 ×20：P 的磁暴／平靜 MAD 比 {mad_ratio:.2f}（< 2）  {'PASS' if c2 else 'FAIL'}",
          file=sys.stderr)
    jump = p.iloc[60:70].median() / p.iloc[70:80].median()
    c3 = abs(jump - 1) < 0.2
    print(f"[self-test] n2 離開前後 P 比 {jump:.2f}（≈1）  {'PASS' if c3 else 'FAIL'}", file=sys.stderr)
    c4 = p.iloc[135:].median() > 1.5 * p.iloc[100:135].median()
    print(f"[self-test] 局部增強 ×3：P 上升 {p.iloc[135:].median() / p.iloc[100:135].median():.2f} 倍  "
          f"{'PASS' if c4 else 'FAIL'}", file=sys.stderr)
    return bool(ok and c2 and c3 and c4)


def run_group(group_id: str, from_daily: bool) -> None:
    cfg = load_group_config(group_id)
    if not cfg.xyz_pool.sufficient:
        print(f"[{group_id}] XYZ 測站池不足，略過", file=sys.stderr)
        return
    path = cfg.interim_dir / "ulf_ref_daily.csv"
    if from_daily and path.exists():
        daily = pd.read_csv(path, dtype={"date": str})
    else:
        daily = build_daily_energy(cfg, list(cfg.xyz_pool.near) + list(cfg.xyz_pool.far))
    build_ratio_index(cfg, daily)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", action="append")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--from-daily", action="store_true", help="用既有的 ulf_ref_daily.csv，不重讀原始資料")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        sys.exit(0 if self_test() else 1)
    if not self_test():
        sys.exit("自我測試失敗，不碰真實資料")
    groups = list(ULF_GROUPS) if args.all else (args.group or [])
    if not groups:
        ap.error("請傳入 --group <id> 或 --all")
    for g in groups:
        run_group(g, args.from_daily)
    if args.all:
        out_dir = PROJECT_DIR / "data" / "interim" / "ulf_ref_ratio"
        out_dir.mkdir(parents=True, exist_ok=True)
        s = summarize(groups)
        (out_dir / "summary.json").write_text(json.dumps(s, indent=2, ensure_ascii=False))
        (out_dir / "summary.md").write_text(render(s), encoding="utf-8")
        print(f"wrote {out_dir / 'summary.md'}", file=sys.stderr)


if __name__ == "__main__":
    main()
