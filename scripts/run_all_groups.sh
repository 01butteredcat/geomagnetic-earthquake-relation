#!/usr/bin/env bash
# Runs scripts/run_pipeline.sh (through verify_pipeline.py, no report build)
# for all 23 groups. A single group failing (e.g. a real verification
# check failure) does NOT abort the batch -- it's recorded and the script
# moves on, since the whole point is to see the full 20-group picture
# including groups with structurally degraded station coverage (G1 has no
# vector stations at all, G2_G3/G4 have thin vector pools; G14-G18 predate
# the modern vector-station network and may be F-pool-only -- see events.py/
# common.py). "Skipped" methods (e.g. G1's H/Z screening) are not failures;
# real verify_pipeline.py check failures are.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-$SCRIPT_DIR/../.venv/bin/python3}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

ALL_GROUP_IDS="G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22 G23"
SUMMARY_JSON="../data/interim/all_groups_run_summary.json"

echo "[" > "$SUMMARY_JSON.tmp"
first=1
overall_status=0

for g in $ALL_GROUP_IDS; do
  echo
  echo "############################################"
  echo "# GROUP $g"
  echo "############################################"

  if ./run_pipeline.sh --group "$g"; then
    run_status="ok"
  else
    run_status="verify_failed"
    overall_status=1
  fi

  # Pull a compact per-group summary out of whatever got written, even on failure
  entry=$("$PYTHON" - "$g" "$run_status" <<PYEOF
import json, sys
from pathlib import Path

group, run_status = sys.argv[1], sys.argv[2]
interim = Path("../data/interim") / group

def load(name):
    p = interim / name
    return json.loads(p.read_text()) if p.exists() else None

sys.path.insert(0, "$SCRIPT_DIR")
from common import load_group_config  # noqa: E402

cfg = load_group_config(group)
verification = load("verification_report.json")
candidates = load("candidate_windows.json")

out = {
    "group": group,
    "run_status": run_status,
    "anchor_date": cfg.anchor_event.date,
    "anchor_magnitude": f"{cfg.anchor_event.magnitude_type}{cfg.anchor_event.magnitude}",
    "coord_confidence": cfg.anchor_event.coord_confidence,
    "xyz_pool_stations": len(cfg.xyz_pool.all_stations),
    "xyz_pool_sufficient": cfg.xyz_pool.sufficient,
    "f_pool_stations": len(cfg.f_pool.all_stations),
    "f_pool_sufficient": cfg.f_pool.sufficient,
    "methods_run": candidates.get("methods_run", []) if candidates else [],
    "verification_summary": verification.get("summary") if verification else None,
}
print(json.dumps(out))
PYEOF
)
  if [ "$first" = 1 ]; then first=0; else echo "," >> "$SUMMARY_JSON.tmp"; fi
  echo "$entry" >> "$SUMMARY_JSON.tmp"
done

echo "]" >> "$SUMMARY_JSON.tmp"
"$PYTHON" -c "import json,sys; d=json.load(open('$SUMMARY_JSON.tmp')); json.dump(d, open('$SUMMARY_JSON','w'), indent=2)"
rm -f "$SUMMARY_JSON.tmp"

echo
echo "############################################"
echo "# SUMMARY (all_groups_run_summary.json)"
echo "############################################"
"$PYTHON" -c "
import json
rows = json.load(open('$SUMMARY_JSON'))
for r in rows:
    v = r['verification_summary'] or {}
    print(f\"{r['group']:10s} status={r['run_status']:14s} methods={','.join(r['methods_run']) or '(none)':8s} \"
          f\"XYZ={r['xyz_pool_stations']:2d}({'ok' if r['xyz_pool_sufficient'] else 'insuff'}) \"
          f\"F={r['f_pool_stations']:2d}({'ok' if r['f_pool_sufficient'] else 'insuff'}) \"
          f\"verify=pass:{v.get('pass','?')}/fail:{v.get('fail','?')}/inconclusive:{v.get('inconclusive','?')}\")
"

exit $overall_status
