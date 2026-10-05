"""把 data/interim/<group>/report_data.json 嵌入 scripts/report_template.html，
並把最終的單檔 artifact 寫到 output/。

用法：build_artifact.py --group G10
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import OUTPUT_DIR, PROJECT_DIR, load_group_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    args = ap.parse_args()
    cfg = load_group_config(args.group)

    template = (PROJECT_DIR / "scripts" / "report_template.html").read_text()
    data_json = (cfg.interim_dir / "report_data.json").read_text()

    token = "/*__REPORT_DATA__*/null"
    assert token in template, "placeholder token missing from template"
    html = template.replace(token, "/*__REPORT_DATA__*/" + data_json)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / f"geomag_precursor_report_{args.group}.html"
    out_path.write_text(html)
    print(f"wrote {out_path} ({out_path.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
