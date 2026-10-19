"""平行 ETL：把每個 <station><date>dsec.sec 檔濃縮成 (a) 一列
每日摘要特徵，以及 (b) 一條重取樣成 1 分鐘的序列，記憶體中同一時間
最多只放一個檔案份量的 1Hz 資料。

用法：build_daily_features.py --group G10

輸出（放在該組自己的 data/interim/<group>/ 底下）：
  daily_features.csv
  minute_series_<station>.parquet   （每站一個檔）
"""
from __future__ import annotations

import argparse
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import list_day_refs, load_group_config  # noqa: E402
from parser import DayFileRef, _parse_day_file_from_text, parse_day_file  # noqa: E402

N_WORKERS = 16

# 尖峰門檻（nT），依 INTERMAGNET 技術手冊 5.2.5 節（Despiking and Removal of
# Artificially Disturbed Data，https://tech-man.intermagnet.org/latest/chapters/processdata/dataqualitycontrol.html）：
# 自然地磁變化達不到每秒 100 nT，所以和前一點、後一點都相差超過 100 nT 的
# 資料點可以安全移除；移除後留成缺值，不內插。（2026-10-20 前用的是專案自訂的
# 單邊 300 nT；實測的故障例子：twu 震後有跳動 1e5 nT 的故障，ttn 長期有
# 最高約 7700 nT 的突波——兩者都是儀器問題。）
SPIKE_THRESHOLD_NT = 100.0


def _despike(series: pd.Series) -> tuple[pd.Series, pd.Series, float]:
    """回傳（尖峰設為 NaN 後的清理序列，只有尖峰樣本才是 True 的
    布林 spike_mask，兩個有效樣本之間的最大 |跳動|）。

    尖峰 = 和前一點、後一點都相差超過 SPIKE_THRESHOLD_NT（手冊 5.2.5 的雙邊
    規則）。第一點、最後一點和緊鄰缺口（哨兵值 NaN）的點只有一個鄰點，和 NaN
    比較都是 False，所以不會被標記。階梯（水位跳動後停在新水位）只有單邊超過
    門檻，不算尖峰，見 _one_sided_jumps。"""
    d_prev = series.diff().abs()
    d_next = series.diff(-1).abs()
    max_jump = float(d_prev.max()) if d_prev.notna().any() else 0.0
    spike_mask = (d_prev > SPIKE_THRESHOLD_NT) & (d_next > SPIKE_THRESHOLD_NT)
    cleaned = series.copy()
    cleaned[spike_mask] = np.nan
    return cleaned, spike_mask, max_jump


def _one_sided_jumps(series: pd.Series) -> int:
    """兩個相鄰樣本都不是尖峰、卻相差超過門檻的次數（階梯／水位跳動）。
    手冊把階梯和尖峰分開處理、沒有給階梯門檻，所以這裡只記錄、不處理。"""
    spike = _despike(series)[1]
    jump = series.diff().abs() > SPIKE_THRESHOLD_NT
    return int((jump & ~spike & ~spike.shift(1, fill_value=False)).sum())


def _list_files(gdms_dir, stations: dict):
    """把該組的日檔分成零散檔（一般 .sec/.sec.gz）和
    來自 tgz 的參照。分開放而不是合成一個清單，是因為
    兩者在 main() 中的處理方式完全不同：零散檔讓工作行程自己開
    很便宜，但 .tgz 成員不是（gzip 無法
    隨機存取——見模組 docstring／建置計畫），所以那些改在
    主行程中取出一次，而不是每個工作行程各取一次。"""
    refs = [r for r in list_day_refs(gdms_dir) if r.station in stations]
    loose = [r for r in refs if r.kind == "loose"]
    tgz = [r for r in refs if r.kind == "tgz"]
    return loose, tgz


def _features_from_df(df, station, date_str):
    n_total = len(df)
    if "F" in df.columns:  # 依它自己的檔頭是純量檔（見 parser.is_scalar_only）
        f_raw = df["F"]
        f, spike_mask, max_jump = _despike(f_raw)
        n_valid = int(f.notna().sum())
        row = {
            "station": station,
            "date": date_str,
            "n_total": n_total,
            "n_valid": n_valid,
            "pct_missing": round(1 - n_valid / n_total, 4),
            "n_spikes": int(spike_mask.sum()),
            "n_one_sided_jumps": _one_sided_jumps(f_raw),
            "max_abs_jump": round(max_jump, 2),
            "F_mean": float(f.mean()) if n_valid else np.nan,
            "F_std": float(f.std()) if n_valid else np.nan,
            "F_range": float(f.max() - f.min()) if n_valid else np.nan,
        }
        minute = df[["F"]].resample("1min").mean()
    else:
        x_raw, y_raw, z_raw = df["X"], df["Y"], df["Z"]
        x, spike_x, max_jump_x = _despike(x_raw)
        y, spike_y, max_jump_y = _despike(y_raw)
        z, spike_z, max_jump_z = _despike(z_raw)
        # 任何分量有突波，那一秒的 H/D 向量就無效；
        # 這是加在原本哨兵值 NaN 缺口之上（而不是取代它）
        spike_any = spike_x | spike_y | spike_z
        x, y, z = x.mask(spike_any), y.mask(spike_any), z.mask(spike_any)
        h = np.sqrt(x**2 + y**2)
        d = np.degrees(np.arctan2(y, x))
        n_valid = int(x.notna().sum())
        row = {
            "station": station,
            "date": date_str,
            "n_total": n_total,
            "n_valid": n_valid,
            "pct_missing": round(1 - n_valid / n_total, 4),
            "n_spikes": int(spike_any.sum()),
            "n_one_sided_jumps": sum(_one_sided_jumps(c) for c in (x_raw, y_raw, z_raw)),
            "max_abs_jump": round(max(max_jump_x, max_jump_y, max_jump_z), 2),
            "H_mean": float(h.mean()) if n_valid else np.nan,
            "H_std": float(h.std()) if n_valid else np.nan,
            "H_range": float(h.max() - h.min()) if n_valid else np.nan,
            "Z_mean": float(z.mean()) if n_valid else np.nan,
            "Z_std": float(z.std()) if n_valid else np.nan,
            "Z_range": float(z.max() - z.min()) if n_valid else np.nan,
            "D_mean": float(d.mean()) if n_valid else np.nan,
        }
        minute = pd.DataFrame({"X": x, "Y": y, "Z": z, "H": h}).resample("1min").mean()

    return row, (station, minute)


def _process_one_ref(ref: DayFileRef):
    """零散檔路徑：工作行程自己做 I/O + 解析，和這個模組
    支援 .tgz 之前一樣（零散檔沒有 gzip 隨機存取的代價，
    所以沒有理由改它）。"""
    try:
        df = parse_day_file(ref, ref.station)
    except Exception as exc:  # noqa: BLE001
        return {"station": ref.station, "date": ref.date_str, "error": str(exc)}, None
    return _features_from_df(df, ref.station, ref.date_str)


def _process_one_text(text: str, ref: DayFileRef):
    """來自 .tgz 的路徑：成員的文字已在
    主行程中取出一次（見 main() 的 tgz 分派迴圈），所以這個工作行程只做
    純 CPU 的解析／去突波／重取樣，完全不碰 tarfile。"""
    try:
        df = _parse_day_file_from_text(text, ref, ref.station)
    except Exception as exc:  # noqa: BLE001
        return {"station": ref.station, "date": ref.date_str, "error": str(exc)}, None
    return _features_from_df(df, ref.station, ref.date_str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    loose_refs, tgz_refs = _list_files(cfg.gdms_dir, cfg.stations)
    n_files = len(loose_refs) + len(tgz_refs)
    print(
        f"[{args.group}] {n_files} 個檔案要處理"
        f"（{len(loose_refs)} 個零散檔，{len(tgz_refs)} 個來自 .tgz）",
        file=sys.stderr,
    )

    # 把來自 .tgz 的參照依壓縮檔分組，讓每個壓縮檔只開一次、
    # 依序走一遍（見模組 docstring：gzip 無法
    # 隨機存取，所以每個成員各自冷開——例如每個工作行程一次——
    # 比單次依序掃描慢約 600 倍，這是在真實的 GDMS
    # 批次下載檔上量到的）。
    by_tgz: dict[Path, dict[str, DayFileRef]] = {}
    for r in tgz_refs:
        by_tgz.setdefault(r.source_path, {})[r.member] = r

    rows = []
    errors = []
    minute_chunks: dict[str, list[pd.DataFrame]] = {s: [] for s in cfg.stations}

    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = [pool.submit(_process_one_ref, r) for r in loose_refs]
        for tgz_path, wanted in by_tgz.items():
            with tarfile.open(tgz_path, "r:gz") as tf:
                for member in tf:  # 對這個壓縮檔做單次依序掃描
                    if member.name not in wanted:
                        continue
                    ref = wanted[member.name]
                    extracted = tf.extractfile(member)
                    text = extracted.read().decode("ascii", errors="strict")
                    futures.append(pool.submit(_process_one_text, text, ref))
        done = 0
        for fut in as_completed(futures):
            row, minute_pair = fut.result()
            done += 1
            if "error" in row:
                errors.append(row)
                print(f"ERROR {row['station']} {row['date']}: {row['error']}", file=sys.stderr)
                continue
            rows.append(row)
            station, minute = minute_pair
            minute_chunks[station].append(minute)
            if done % 200 == 0:
                print(f"  {done}/{n_files}", file=sys.stderr)

    daily = pd.DataFrame(rows).sort_values(["station", "date"])
    daily.to_csv(cfg.interim_dir / "daily_features.csv", index=False)
    print(f"已寫入 daily_features.csv（{len(daily)} 列）", file=sys.stderr)

    if errors:
        pd.DataFrame(errors).to_csv(cfg.interim_dir / "parse_errors.csv", index=False)
        print(f"已寫入 parse_errors.csv（{len(errors)} 列）", file=sys.stderr)

    for station, chunks in minute_chunks.items():
        if not chunks:
            continue
        merged = pd.concat(chunks).sort_index()
        out_path = cfg.interim_dir / f"minute_series_{station}.parquet"
        merged.to_parquet(out_path)
        print(f"wrote {out_path.name} ({len(merged)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
