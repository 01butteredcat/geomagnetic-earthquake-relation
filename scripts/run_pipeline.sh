#!/usr/bin/env bash
# Orchestrates the geomag_precursor pipeline for one group end to end and,
# critically, refuses to (re)build a report if verify_pipeline.py finds any
# failing check -- nothing used to gate prepare_report_data.py/
# build_artifact.py on verification passing, so a stale or failing
# verification_report.json could silently sit underneath a freshly-built
# report.
#
# Usage: run_pipeline.sh --group G10 [--full-report]
#
# By default this runs steps 1-6 (through verify_pipeline.py) -- that's all
# run_all_groups.sh needs for the 20-group-folder/23-event-sequence batch/
# cross-group analysis.
# --full-report additionally runs prepare_report_data.py + build_artifact.py
# to (re)build the polished single-event HTML report; that report's
# hand-written narrative text is specific to G10/2024, so --full-report is
# only meaningful for --group G10.
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
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
if [ -z "$GROUP" ]; then
  echo "usage: run_pipeline.sh --group <G1|G2|...|G23> [--full-report]" >&2
  exit 1
fi

STEPS=6
[ "$FULL_REPORT" = 1 ] && STEPS=8

echo "== [$GROUP] 1/$STEPS build_daily_features.py =="
"$PYTHON" build_daily_features.py --group "$GROUP"

echo "== [$GROUP] 2/$STEPS timezone_check.py =="
"$PYTHON" timezone_check.py --group "$GROUP"

echo "== [$GROUP] 3/$STEPS fetch_space_weather.py =="
# Only a complete Kp+Dst result ("confidence": "high...") is reused; a partial or
# failed earlier fetch (transient Kyoto/GFZ error) is retried instead of being
# cached forever -- that is how G1/G2/G3/G11 stayed without Dst for weeks.
if [ -f "../data/interim/$GROUP/storm_days.csv" ] && grep -q '"confidence": "high' "../data/interim/$GROUP/storm_days_summary.json" 2>/dev/null; then
  echo "  storm_days.csv already cached with high confidence, skipping fetch (delete it to force a re-fetch)"
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
  echo "[$GROUP] verify_pipeline.py FAILED -- stopping before regenerating any report." >&2
  echo "Fix the failing check(s) above, then re-run this script." >&2
  exit 1
fi

if [ "$FULL_REPORT" = 1 ]; then
  echo "== [$GROUP] 7/$STEPS prepare_report_data.py =="
  "$PYTHON" prepare_report_data.py --group "$GROUP"

  echo "== [$GROUP] 8/$STEPS build_artifact.py =="
  "$PYTHON" build_artifact.py --group "$GROUP"
fi

echo
echo "[$GROUP] pipeline complete: all verification checks passed$( [ "$FULL_REPORT" = 1 ] && echo ' and the report was rebuilt' )."
