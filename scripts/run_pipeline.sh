#!/usr/bin/env bash
# 從頭到尾統籌單一組的 geomag_precursor 流程，而且
# 關鍵在於：只要 verify_pipeline.py 有任何檢查失敗，就拒絕（重新）建置報告
# ——以前沒有任何東西讓 prepare_report_data.py/
# build_artifact.py 以驗證通過為前提，所以過時或失敗的
# verification_report.json 可能默默躺在剛建好的
# 報告底下。
#
# 用法：run_pipeline.sh --group G10 [--full-report]
#
# 預設跑步驟 1-6（到 verify_pipeline.py 為止）——這就是
# run_all_groups.sh 做 20 個資料夾／23 個事件序列批次
# 跨組分析所需要的全部。
# --full-report 另外跑 prepare_report_data.py + build_artifact.py，
# （重新）建置精修的單一事件 HTML 報告；那份報告
# 手寫的敘述文字是專門針對 G10/2024 的，所以 --full-report
# 只對 --group G10 有意義。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-$SCRIPT_DIR/../.venv/bin/python3}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

GROUP=""
FULL_REPORT=0
while [ $# -gt 0 ]; do
  case "$1" in
    --group) GROUP="$2"; shift 2 ;;
    --full-report) FULL_REPORT=1; shift ;;
    *) echo "未知的參數：$1" >&2; exit 1 ;;
  esac
done
if [ -z "$GROUP" ]; then
  echo "用法：run_pipeline.sh --group <G1|G2|...|G24> [--full-report]" >&2
  exit 1
fi

STEPS=6
[ "$FULL_REPORT" = 1 ] && STEPS=8

echo "== [$GROUP] 1/$STEPS build_daily_features.py =="
"$PYTHON" build_daily_features.py --group "$GROUP"

echo "== [$GROUP] 2/$STEPS timezone_check.py =="
"$PYTHON" timezone_check.py --group "$GROUP"

echo "== [$GROUP] 3/$STEPS fetch_space_weather.py =="
# 只有完整的 Kp+Dst 結果（"confidence": "high..."）才會重用；先前抓到一半或
# 失敗的（Kyoto/GFZ 暫時性錯誤）會重抓，而不是被
# 永遠快取——G1/G2/G3/G11 就是這樣好幾週都沒有 Dst。快取
# 也必須涵蓋資料夾目前的日期範圍（G6/G7/G8/G17 漏了最後幾週）。
if "$PYTHON" fetch_space_weather.py --group "$GROUP" --check-cache; then
  echo "  storm_days.csv 已有高信心且完整涵蓋的快取，跳過抓取（刪掉它可以強制重抓）"
else
  "$PYTHON" fetch_space_weather.py --group "$GROUP"
fi

echo "== [$GROUP] 4/$STEPS compute_indices.py =="
"$PYTHON" compute_indices.py --group "$GROUP"

echo "== [$GROUP] 5/$STEPS ulf_analysis.py =="
"$PYTHON" ulf_analysis.py --group "$GROUP"

echo "== [$GROUP] 6/$STEPS verify_pipeline.py =="
if ! "$PYTHON" verify_pipeline.py --group "$GROUP"; then
  echo
  echo "[$GROUP] verify_pipeline.py 失敗——在重新產生任何報告之前停止。" >&2
  echo "修好上面失敗的檢查後，再重跑這支腳本。" >&2
  exit 1
fi

if [ "$FULL_REPORT" = 1 ]; then
  echo "== [$GROUP] 7/$STEPS prepare_report_data.py =="
  "$PYTHON" prepare_report_data.py --group "$GROUP"

  echo "== [$GROUP] 8/$STEPS build_artifact.py =="
  "$PYTHON" build_artifact.py --group "$GROUP"
fi

echo
echo "[$GROUP] 流程完成：所有驗證檢查都通過$( [ "$FULL_REPORT" = 1 ] && echo '，報告也已重建' )。"
