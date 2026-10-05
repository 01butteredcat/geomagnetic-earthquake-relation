"""把 data/interim/validation_report_data.json 嵌入 scripts/report_
template_validation.html，並把最終的單檔 artifact 寫到
output/——是 build_artifact.py 的跨組版本。

用法：build_validation_report.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import OUTPUT_DIR, PROJECT_DIR  # noqa: E402


def main():
    template = (PROJECT_DIR / "scripts" / "report_template_validation.html").read_text()
    data_json = (PROJECT_DIR / "data" / "interim" / "validation_report_data.json").read_text()

    token = "/*__REPORT_DATA__*/null"
    assert token in template, "placeholder token missing from template"
    html = template.replace(token, "/*__REPORT_DATA__*/" + data_json)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / "geomag_precursor_validation_report.html"
    out_path.write_text(html)
    print(f"wrote {out_path} ({out_path.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
