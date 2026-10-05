#!/usr/bin/env bash
# 對全部 24 組執行 scripts/run_pipeline.sh（到 verify_pipeline.py 為止，不建置報告）。
# 單一組失敗（例如真正的驗證
# 檢查失敗）**不會**中止整批——會記錄下來，腳本
# 繼續往下跑，因為重點就是看到完整的 20 組全貌，
# 包括測站涵蓋在結構上較差的組別（G1 完全沒有
# 向量站，G2_G3/G4 的向量測站池很薄；G14-G18 早於
# 現代向量測站網，可能只有 F 測站池——見 events.py／
# common.py）。「跳過」的方法（例如 G1 的 H/Z 篩檢）不算失敗；
# verify_pipeline.py 真正的檢查失敗才算。
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-$SCRIPT_DIR/../.venv/bin/python3}"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

ALL_GROUP_IDS="G1 G2 G3 G4 G5 G6 G7 G8 G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22 G23 G24"
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

  # 從已寫出的任何東西中取出精簡的逐組摘要，即使失敗也照做
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
echo "# 摘要（all_groups_run_summary.json）"
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
