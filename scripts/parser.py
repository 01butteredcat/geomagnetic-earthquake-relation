"""GDMSdata 1 秒地磁檔的 IAGA-2002 解析器。

處理整個 G1-G13 資料集中出現的兩種測站慣例：
  - 向量站：真實資料在 X/Y/Z，F 是結構性佔位值
    88888.00（未回報）。
  - 純量站：真實資料在 F，X/Y/Z 永遠是 88888.00。

適用哪一種慣例**不是**測站代碼的固定屬性——大多數
測站一開始只有純量，在資料集歷史中途才升級成向量
（依 G1-G13 測站登錄調查，交接總是落在組別年代的
邊界上；`ttn` 是唯一一個從未升級的
長期測站）。每個檔案自己的檔頭會透過
`Reported` 那一行（`F` 或 `XYZF`）說明適用哪一種，所以 `parse_day_file` 逐檔讀取它，
而不是查靜態登錄表——見 `parse_header`/`is_scalar_only`。

兩種慣例也都用 99999.00 當另一個哨兵值，代表真的
資料中斷（和 88888.00 的「通道未回報」意思不同），它會
在該測站「真實」的那些欄位中被對應成 NaN。
"""
from __future__ import annotations

import gzip
import io
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

OUTAGE_SENTINEL = 99999.0
NOT_REPORTED_PLACEHOLDER = 88888.0

_HEADER_LINE_RE = re.compile(r"^DATE\s+TIME\s+DOY\s+")
_HEADER_FIELD_RE = re.compile(r"^\s*(.+?)\s{2,}(.*?)\s*\|\s*$")


@dataclass(frozen=True)
class DayFileRef:
    """指向一個測站–日份量的 .sec 資料，不管它
    是零散檔案，還是 GDMS 批次下載 .tgz 裡的
    成員。只有 str/Path 欄位——可以安全地 pickle 跨越
    ProcessPoolExecutor 的行程邊界（不像開啟中的 tarfile handle）。
    `member` 對零散檔是 None，對 .tgz 來源是 tar 成員名稱。
    """

    station: str
    date_str: str  # YYYYMMDD 格式
    source_path: Path  # 零散檔：.sec/.sec.gz 本身；tgz：.tgz 檔
    member: str | None = None

    @property
    def kind(self) -> str:
        return "tgz" if self.member is not None else "loose"

    @property
    def label(self) -> str:
        """用於記錄／錯誤訊息，取代只印出一個裸路徑。"""
        return f"{self.source_path}!{self.member}" if self.member else str(self.source_path)


def _as_ref(x: "str | Path | DayFileRef", station_code: str | None = None) -> DayFileRef:
    """把零散路徑（str/Path，每個呼叫者過去都是這樣傳的）
    正規化成 DayFileRef，讓這個模組其餘部分只有一條程式路徑。"""
    if isinstance(x, DayFileRef):
        return x
    path = Path(x)
    stem = path.name[:-3] if path.name.endswith(".sec.gz") else path.name
    return DayFileRef(station_code or stem[:3], stem[3:11], path, None)


def _read_text(ref: DayFileRef) -> str:
    """把這個參照的完整檔案內容讀進記憶體，只讀一次。一個 .sec
    檔解壓後約 6MB，小到一開始就整個讀進來
    （而不是每一輪都重開／重解壓）對零散檔很便宜，
    對 .tgz 成員則是**必要**的：gzip 不支援隨機
    存取，所以每次讀取都重開同一個 .tgz 會比讀一次
    昂貴得多（避免每個工作行程每個檔案都這樣做的批次取出路徑，
    見 build_daily_features.py）。"""
    if ref.member is None:  # 零散的 .sec 或 .sec.gz
        path = ref.source_path
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="ascii", errors="strict") as f:
                return f.read()
        with open(path, "r", encoding="ascii", errors="strict") as f:
            return f.read()
    # .tgz 成員
    with tarfile.open(ref.source_path, "r:gz") as tf:
        extracted = tf.extractfile(ref.member)
        if extracted is None:
            raise ValueError(f"{ref.label}：tar 成員不是一般檔案")
        raw = extracted.read()
    return raw.decode("ascii", errors="strict")


def open_raw(ref: "str | Path | DayFileRef"):
    """以文字模式開啟 .sec/.sec.gz 檔或 .tgz 成員，回傳一個
    可以當 context manager、可逐行迭代的物件——介面和以前相同
    （向後相容傳入裸 str/Path 的呼叫者，例如
    verify_pipeline.py 的抽查）。內部現在一律一開始就讀進整個
    檔案／成員（見 `_read_text`），而不是回傳一個即時的
    gzip／一般檔案 handle。"""
    return io.StringIO(_read_text(_as_ref(ref)))


def _find_data_start_from_text(text: str, label: str) -> int:
    """回傳第一列資料的行號（從 0 起算）。"""
    for i, line in enumerate(text.splitlines(keepends=True)):
        if _HEADER_LINE_RE.match(line):
            return i + 1
    raise ValueError(f"找不到 IAGA-2002 欄位標題列：{label}")


def _parse_header_from_text(text: str, label: str) -> dict:
    fields: dict[str, str] = {}
    for line in text.splitlines(keepends=True):
        if _HEADER_LINE_RE.match(line):
            break
        m = _HEADER_FIELD_RE.match(line)
        if m:
            fields[m.group(1).strip()] = m.group(2).strip()

    if "Geodetic Latitude" not in fields or "Geodetic Longitude" not in fields:
        raise ValueError(f"{label}：在 IAGA-2002 檔頭中找不到 lat/lon")

    elevation = fields.get("Elevation", "")
    return {
        "station_name": fields.get("Station Name", ""),
        "lat": float(fields["Geodetic Latitude"]),
        "lon": float(fields["Geodetic Longitude"]),
        "elevation_m": float(elevation) if elevation else None,
        "reported": fields.get("Reported", ""),
        "source": fields.get("Source of Data", ""),
    }


def parse_header(ref: "str | Path | DayFileRef") -> dict:
    """把一個 .sec 檔／tgz 成員的 IAGA-2002 檔頭解析成測站
    中繼資料。

    每個檔案都帶有自己的測站名稱／緯度／經度／高程，以及
    `Reported` 欄位（`F` 或 `XYZF`）——已在所有 13 組中對每個測站代碼
    抽查多個檔案，確認同一測站／年代在不同日期間一致。
    這代表測站中繼資料和純量 vs 向量的
    分類可以直接從資料推出，而不必
    在另一張查詢表中人工維護。

    回傳的 dict 含：station_name、lat、lon、elevation_m（如果
    檔頭欄位空白則為 None——幾個較舊／已停用的測站有這種情況）、
    reported（原始字串，例如 "F" 或 "XYZF"）、source。
    """
    ref = _as_ref(ref)
    return _parse_header_from_text(_read_text(ref), ref.label)


def is_scalar_only(ref: "str | Path | DayFileRef") -> bool:
    """如果這個檔案自己的檔頭只回報 F（沒有真實 X/Y/Z）則為 True。"""
    return parse_header(ref)["reported"].strip().upper() == "F"


def _parse_day_file_from_text(text: str, ref: DayFileRef, station_code: str) -> pd.DataFrame:
    """核心解析邏輯，處理已經讀進記憶體的文字，而不是
    從磁碟／tar 重新讀取。獨立公開（而不只是內嵌的輔助函式），
    是因為 build_daily_features.py 的批次 .tgz 路徑會在主行程中
    每個檔案取出一次成員文字，直接交給工作
    行程，省掉這個函式原本要做的逐檔 I/O。"""
    skiprows = _find_data_start_from_text(text, ref.label)
    scalar_only = _parse_header_from_text(text, ref.label)["reported"].strip().upper() == "F"

    df = pd.read_csv(
        io.StringIO(text),
        skiprows=skiprows,
        sep=r"\s+",
        header=None,
        names=["DATE", "TIME", "DOY", "X", "Y", "Z", "F"],
        dtype={"X": "float32", "Y": "float32", "Z": "float32", "F": "float32"},
        engine="c",
    )

    if len(df) != 86400:
        raise ValueError(
            f"{ref.label}：預期 86400 列資料，實際 {len(df)}"
        )

    ts = pd.to_datetime(df["DATE"] + " " + df["TIME"], format="%Y-%m-%d %H:%M:%S.%f")

    if scalar_only:
        real = df[["F"]].copy()
        real.loc[real["F"] == OUTAGE_SENTINEL, "F"] = pd.NA
    else:
        real = df[["X", "Y", "Z"]].copy()
        for col in ("X", "Y", "Z"):
            real.loc[real[col] == OUTAGE_SENTINEL, col] = pd.NA
            # 防禦性處理：向量欄位中零星出現的未回報佔位值
            # 也是無效資料，雖然依實際資料不預期會出現。
            real.loc[real[col] == NOT_REPORTED_PLACEHOLDER, col] = pd.NA

    real = real.astype("float32")
    real.index = ts
    real.index.name = "time"
    return real


def parse_day_file(ref: "str | Path | DayFileRef", station_code: str) -> pd.DataFrame:
    """把一個 <station><YYYYMMDD>dsec.sec 檔（零散檔或 .tgz 成員）
    解析成 DataFrame。

    回傳以不含時區的 UTC datetime 為索引的 DataFrame（檔案中印出的
    原始時鐘值；真正的 UTC 偏移另外由
    timezone_check.py 判定），float32 欄位取決於這個特定
    檔案的檔頭回報了什麼（見 `is_scalar_only`）：
      - 向量檔：X、Y、Z
      - 純量檔：F
    缺值／中斷樣本（哨兵值 99999.00）是 NaN。`station_code` 不會
    用來判斷這件事（只保留給呼叫者記帳／錯誤訊息用）
    ——純量 vs 向量的狀態每次都從這個檔案自己的檔頭重新讀取。
    """
    ref = _as_ref(ref, station_code=station_code)
    return _parse_day_file_from_text(_read_text(ref), ref, station_code)


if __name__ == "__main__":
    import sys

    p = Path(sys.argv[1])
    station = p.name[:3]
    out = parse_day_file(p, station)
    print(out.describe())
    print(out.head())
    print("NaN 個數：\n", out.isna().sum())
