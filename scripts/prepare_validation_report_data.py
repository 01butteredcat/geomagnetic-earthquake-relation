"""Gather surrogate_test.py / superposed_epoch_analysis.py / backtest_rule.py
outputs (already written under data/interim/) into one JSON blob that
report_template_validation.html embeds verbatim -- the cross-group
counterpart of prepare_report_data.py, which does the same thing for a
single group's descriptive report.

Usage: prepare_validation_report_data.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import PROJECT_DIR  # noqa: E402

ULF_GROUPS = ("G4", "G5", "G6_G7_G8", "G9", "G10", "G11", "G12", "G13", "G19", "G20", "G23")
BANDS = ("pc3", "pc4")


def load_json(path: Path):
    return json.loads(path.read_text()) if path.exists() else None


def main():
    interim = PROJECT_DIR / "data" / "interim"
    external = PROJECT_DIR / "data" / "external"

    surrogate = {}
    for g in ULF_GROUPS:
        for band in BANDS:
            p = interim / "surrogate_test" / f"{g}_{band}.json"
            data = load_json(p)
            if data is not None:
                surrogate[f"{g}_{band}"] = data

    sea = {}
    backtest = {}
    catalog_counts = {}
    for label in ("m5.5", "m5.0"):
        for band in BANDS:
            p = interim / "superposed_epoch" / f"{band}_stack_{label}.json"
            data = load_json(p)
            if data is not None:
                sea[f"{band}_{label}"] = data
        p = interim / "backtest" / f"z_rule_backtest_{label}.json"
        data = load_json(p)
        if data is not None:
            backtest[label] = data

        catalog_path = external / f"extended_catalog_{label}.csv"
        n_raw = n_decl = 0
        if catalog_path.exists():
            with catalog_path.open(newline="") as f:
                for row in csv.DictReader(f):
                    n_raw += 1
                    if row["declustered"] == "True":
                        n_decl += 1
        catalog_counts[label] = {"n_raw": n_raw, "n_declustered": n_decl}

    report_data = {
        "meta": {
            "ulf_groups": list(ULF_GROUPS),
            "catalog_counts": catalog_counts,
            "fixed_rule_threshold": -4.1,
            "windows_days": [7, 14, 30],
        },
        "surrogate": surrogate,
        "sea": sea,
        "backtest": backtest,
    }

    out_path = interim / "validation_report_data.json"
    out_path.write_text(json.dumps(report_data, indent=2))
    print(f"wrote {out_path} ({out_path.stat().st_size / 1024:.0f} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
