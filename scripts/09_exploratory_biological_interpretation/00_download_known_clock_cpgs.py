#!/usr/bin/env python3

"""Download published clock CpG lists from pyaging."""

import csv
import re
from pathlib import Path

import pyaging as pya


output_dir = Path("data/reference")
pyaging_data_dir = output_dir / "pyaging_clock_files"
known_clock_file = output_dir / "known_epigenetic_clock_cpgs.csv"
summary_file = output_dir / "known_epigenetic_clock_download_summary.csv"

clock_names = [
    "Horvath2013",
    "Hannum",
    "PhenoAge",
    "DNAmPhenoAge",
    "GrimAge",
    "SkinAndBlood",
    "ZhangEN",
]

output_dir.mkdir(parents=True, exist_ok=True)
pyaging_data_dir.mkdir(parents=True, exist_ok=True)

logger = pya.logger.Logger("known_clock_cpg_download")
known_clock_rows = []
summary_rows = []

for clock_name in clock_names:
    try:
        clock = pya.pred.load_clock(
            clock_name,
            "cpu",
            str(pyaging_data_dir),
            logger,
            indent_level=1,
        )

        features = list(clock.features)
        cpgs = sorted({
            str(feature)
            for feature in features
            if re.match(r"^cg[0-9]{8}$", str(feature))
        })

        for cpg in cpgs:
            known_clock_rows.append({
                "clock_name": clock_name,
                "CpG": cpg,
            })

        summary_rows.append({
            "clock_name": clock_name,
            "status": "downloaded",
            "features_total": len(features),
            "cpg_features": len(cpgs),
        })

    except Exception as error:
        summary_rows.append({
            "clock_name": clock_name,
            "status": "failed",
            "features_total": "",
            "cpg_features": "",
            "error": str(error),
        })

if len(known_clock_rows) == 0:
    raise SystemExit("No CpG features were downloaded from pyaging")

with known_clock_file.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["clock_name", "CpG"])
    writer.writeheader()
    writer.writerows(known_clock_rows)

with summary_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "clock_name",
            "status",
            "features_total",
            "cpg_features",
            "error",
        ],
    )
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"Wrote {known_clock_file}")
print(f"Wrote {summary_file}")
