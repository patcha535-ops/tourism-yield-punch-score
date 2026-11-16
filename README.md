# Tourism yield and spatial concentration in Thailand (code, data, manuscript sources)

Code, derived data and LaTeX sources for *"Tourism Yield and Spatial Concentration in Thailand: A Provincial Assessment Using the
Punch Score"* (77 provinces, 2024).

## Contents
| Path | What it is |
|---|---|
| `build_dataset.py` | Extracts the 77-province dataset from the compiled data tables (`data_tourism_2024.pdf`) and validates it against the published regional and national totals |
| `te_analysis.py` | All computations: Punch Score, Gini, Theil and between-region decomposition, within-region bootstrap, visitor-composition decomposition, global and local Moran, Spearman correlations, robustness, tables and figures |
| `verify_claims.py` | Compares every number quoted in the manuscript with the recomputation (171 checks) |
| `run_all_te.sh` | Reproduces everything |
| `provinces_2024.csv` | The analysis dataset (77 rows) |
| `results_te/` | Outputs the manuscript is built from (LaTeX tables, figures, `results.json`, manifest) |
| `main_te.tex`, `title_page.tex`, `supplementary_te.tex` (+ PDFs) | Manuscript, title page, supplementary material (reproduces all code inline) |

## Reproduce
```
pip install -r requirements.txt            # also needs the pdftotext command (poppler-utils)
sh run_all_te.sh data_tourism_2024.pdf     # about one minute
python make_supplement_te.py               # regenerates supplementary_te.tex
pdflatex main_te.tex && pdflatex main_te.tex    # twice
```
Default seed 20260101; run identifier and checksums are in `results_te/manifest.json`.

## Data notice
The data are compiled from public sources: Ministry of Tourism and Sports (provincial tourism statistics, revised 28 August 2025;
Domestic Tourism Statistics 2024), NESDC (gross provincial product), National Statistical Office (registered population, area) and
OCHA (province boundaries, used for centroids). Check each provider's terms of reuse before redistributing. The MIT licence covers the code only.

## Archiving / citing
Create a GitHub *release* with a tag that does not exist yet (e.g. `v1.0.0`), let Zenodo archive it, and use the DOI of that release
in the manuscript (`\RepoDOI` in `main_te.tex` and `title_page.tex`). Keep the version in `CITATION.cff` equal to the tag.
