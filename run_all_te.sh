#!/bin/sh
# Reproduces every number, table and figure of the Tourism Economics manuscript.
# Usage: sh run_all_te.sh data_tourism_2024.pdf      (needs pdftotext, python3 with numpy, scipy, pandas, matplotlib)
set -e
python3 build_dataset.py "${1:-data_tourism_2024.pdf}" provinces_2024.csv
python3 te_analysis.py provinces_2024.csv --outdir results_te
python3 verify_claims.py results_te/results.json
