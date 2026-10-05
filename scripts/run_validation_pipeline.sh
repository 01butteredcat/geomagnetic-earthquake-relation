#!/usr/bin/env bash
# 從頭到尾統籌教授建議的驗證方法：
# 擴充地震目錄 -> 替代資料顯著性檢定 -> 疊加
# 時間分析 -> z<=-4.1 規則回測 -> 跨組驗證報告。
# 假設 ULF_GROUPS 中的每一組都已經跑過 run_pipeline.sh（步驟 1-6），
# 也就是每組的 ulf_near_far_index.csv 都已存在——
# 這支腳本不會重新推導原始測站資料，只是在
# run_all_groups.sh 已產出的東西之上建立新的驗證分析。
#
# 全程比較三個規模級距（嚴格 -> 寬鬆）：M>=6.0
# （N 小，每個事件的訊雜比最好）、M>=5.5（主要）、M>=5.0（N
# 最大，每個事件的訊號最弱）——這是分層的敏感度設計，而不是
# 只挑一個門檻，因為兩個極端在樣本數
# 和每個事件訊號強度之間的取捨方向相反。
#
# 用法：run_validation_pipeline.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-$SCRIPT_DIR/../.venv/bin/python3}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

echo "== 1/8 fetch_earthquake_catalog.py（M>=6.0，嚴格級距） =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 6.0

echo "== 2/8 fetch_earthquake_catalog.py（M>=5.5，主要） =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 5.5

echo "== 3/8 fetch_earthquake_catalog.py（M>=5.0，敏感度） =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 5.0

echo "== 4/8 surrogate_test.py（所有組別） =="
"$PYTHON" surrogate_test.py --all

# SEA 只跑 M>=6.0（2026-09-30）：用完整的 CWA 目錄時，M5 級地震
# 密集到 14/14（M>=5.0）和 13/14（M>=5.5）組都沒有任何一天距離每一起
# 地震 >=30 天，也就是完全沒有和地震無關的虛無時期。回測不需要虛無
# 時期，三個級距都保留。
echo "== 5/8 superposed_epoch_analysis.py (M>=6.0) =="
"$PYTHON" superposed_epoch_analysis.py --catalog ../data/external/extended_catalog_m6.0.csv --label m6.0 --min-mag 6.0
rm -f ../data/interim/superposed_epoch/pc[34]_stack_m5.[05].json

echo "== 6/8 backtest_rule.py（M>=6.0、M>=5.5 和 M>=5.0） =="
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m6.0.csv --label m6.0 --min-mag 6.0
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m5.5.csv --label m5.5 --min-mag 5.5
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m5.0.csv --label m5.0 --min-mag 5.0

echo "== 7/8 prepare_validation_report_data.py =="
"$PYTHON" prepare_validation_report_data.py

echo "== 8/8 build_validation_report.py =="
"$PYTHON" build_validation_report.py

echo
echo "驗證流程完成：output/geomag_precursor_validation_report.html"
