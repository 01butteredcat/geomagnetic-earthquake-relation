#!/usr/bin/env bash
# Orchestrates the professor-suggested validation methodology end to end:
# extended earthquake catalog -> surrogate significance test -> superposed
# epoch analysis -> z<=-4.1 rule backtest -> cross-group validation report.
# Assumes run_pipeline.sh (steps 1-6) has already been run for every group
# in ULF_GROUPS, i.e. each group's ulf_near_far_index.csv already exists --
# this script does not re-derive raw station data, it only builds the new
# validation analyses on top of what run_all_groups.sh already produced.
#
# Three magnitude tiers are compared throughout (strict -> loose): M>=6.0
# (small N, best per-event signal-to-noise), M>=5.5 (primary), M>=5.0 (largest
# N, weakest per-event signal) -- a stratified sensitivity design rather than
# picking a single threshold, since the two extremes trade off sample size
# against per-event signal strength in opposite directions.
#
# Usage: run_validation_pipeline.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-$SCRIPT_DIR/../.venv/bin/python3}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

echo "== 1/8 fetch_earthquake_catalog.py (M>=6.0, strict tier) =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 6.0

echo "== 2/8 fetch_earthquake_catalog.py (M>=5.5, primary) =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 5.5

echo "== 3/8 fetch_earthquake_catalog.py (M>=5.0, sensitivity) =="
"$PYTHON" fetch_earthquake_catalog.py --min-mag 5.0

echo "== 4/8 surrogate_test.py (all groups) =="
"$PYTHON" surrogate_test.py --all

echo "== 5/8 superposed_epoch_analysis.py (M>=6.0, M>=5.5 and M>=5.0) =="
"$PYTHON" superposed_epoch_analysis.py --catalog ../data/external/extended_catalog_m6.0.csv --label m6.0
"$PYTHON" superposed_epoch_analysis.py --catalog ../data/external/extended_catalog_m5.5.csv --label m5.5
"$PYTHON" superposed_epoch_analysis.py --catalog ../data/external/extended_catalog_m5.0.csv --label m5.0

echo "== 6/8 backtest_rule.py (M>=6.0, M>=5.5 and M>=5.0) =="
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m6.0.csv --label m6.0
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m5.5.csv --label m5.5
"$PYTHON" backtest_rule.py --catalog ../data/external/extended_catalog_m5.0.csv --label m5.0

echo "== 7/8 prepare_validation_report_data.py =="
"$PYTHON" prepare_validation_report_data.py

echo "== 8/8 build_validation_report.py =="
"$PYTHON" build_validation_report.py

echo
echo "validation pipeline complete: output/geomag_precursor_validation_report.html"
