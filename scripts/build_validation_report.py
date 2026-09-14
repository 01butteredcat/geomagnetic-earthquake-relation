"""Embed data/interim/validation_report_data.json into scripts/report_
template_validation.html and write the final self-contained artifact to
output/ -- the cross-group counterpart of build_artifact.py.

Usage: build_validation_report.py
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
