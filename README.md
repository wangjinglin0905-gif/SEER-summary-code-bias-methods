# Bias from treatment summary codes

Version 1.0.3 — 22 September 2026

Reproducibility materials for **Bias and negative-control behavior when treatment sequences are observed only as summary codes: an ADEMP simulation study with two SEER illustrations**.

This package contains analysis code, synthetic simulation outputs, figure source tables, aggregate-only SEER results and historical verification records. The manuscript is undergoing author review. No individual-level SEER data, ethics documents or manuscript files are included. GitHub release publication and Zenodo archiving are distinct; do not infer a DOI from this package.

## Changes in 1.0.3

- Corrected Figure 6 panel a from "No prior therapy" to "RT/CT codes" in the R source and all four exports. The filter retains radiation coded None/Unknown and chemotherapy coded No/Unknown; these categories do not establish absence of treatment or reconstruct treatment timing.
- No cohort counts, estimates, analysis inputs or statistical code changed. This is a scientific-label correction, not a statistical rerun. The other five figures are unchanged.
- The unchanged large synthetic archives remain available with [release v1.0.2](https://github.com/wangjinglin0905-gif/SEER-summary-code-bias-methods/releases/tag/v1.0.2). Use them together with the v1.0.3 code package; their original checksums remain applicable.
- See `verification/release_delta_checks_v1.0.3.json` for the bounded source/export and unchanged-data checks. Earlier releases remain available.

## Changes in 1.0.2

- Added author-confirmed affiliations, ORCIDs and licenses.
- Integrated the approved Figure 2, 3 and 6 spacing changes into the portable R script and exports. Figure 1 retains its approved eight-node, ten-edge structure.
- Designated `analysis/strict_qba/strict_qba_summary_36_rows.csv` and its matching `source_data/figure4_qba_cell_results.csv` as the current QBA summaries. The unused legacy aliases `source_data/qba_36_scenario_prior_rows.csv` and `analysis/qba_prior_summary_verified.csv` are not distributed in this version. Historical inputs explicitly under `data/legacy/` remain only for traceable comparisons.
- Retained all scientific results and simulation source code. No new statistical simulation was run for this release. Earlier verification reports are dated historical evidence, not new full-pipeline validation.

## What is reproducible

- 72 primary/case-informed simulation scenarios and 252,000 replicate datasets;
- seven H2 negative-control boundary scenarios;
- four QBA calibration distributions evaluated in nine target cells;
- corrected five-year partial restricted-mean truths and rescoring;
- the six submission figures from frozen source tables;
- a 1,499-check independent aggregation audit;
- an independent 1,200-fit comparison with `survival::coxph(..., robust=TRUE)`.

## Repository map

- `analysis/`: frozen aggregate and source-level results used by the manuscript;
- `code/python/`: strict Cox, H2, QBA, pRMST, and audit scripts;
- `code/R/`: figure generation, Cox cross-implementation validation, and restricted-data SEER analysis;
- `source_data/`: concise figure/table source data; the obsolete pre-correction `rmst_12_scenarios.csv` is deliberately excluded;
- `figures/`: PNG, TIFF, PDF, and SVG exports for Figures 1–6;
- `verification/`: audit reports, check ledgers, R session information, and the 1,200-fit comparison;
- `data/synthetic/rmst_archived_intervals/`: small archived synthetic interval outputs required for pRMST rescoring;
- `data/seer_aggregate/`: aggregate-only duplicate of the SEER reporting outputs for convenient discovery;
- `docs/`: audit, data-access, reporting, and pre-upload documentation.

## Large synthetic archives

The following files are prepared separately for Zenodo and are not kept in Git history:

1. `synthetic_primary_h2_replicates.zip` — all 72 primary and seven H2 replicate files;
2. `qba_replicate_and_calibration_outputs.zip` — 180,000 QBA evaluation rows and independent oracle-calibration runs;
3. `cox_validation_source_1200_datasets.zip` — the large synthetic source datasets used for the 1,200 R/Python cross-implementation fits.

After downloading, extract each archive at the repository root. The internal paths are already arranged under `analysis/` or `data/large/`.

On Windows, choose a short working-directory path. Deeply nested extraction can exceed legacy path-length limits; the current clean-extraction test used a short directory. Do not interpret an extraction path error as an analytical failure.

## Python environment

Use Python 3.12.13 with NumPy 2.3.5 and pandas 3.0.1 to match the recorded Cox/QBA execution environment. The package requirements now pin those two dependencies; the previous `pandas<3` range was inconsistent with the execution provenance. This metadata correction is not a claim that the entire pipeline was rerun in a fresh environment on 8 September. The current figure-only smoke test and source syntax checks are documented separately in `verification/release_delta_checks_v1.0.1.json`.

```bash
python -m venv .venv
python -m pip install -r environment/requirements.txt
```

Full strict reruns are computationally intensive. From the repository root:

```bash
python code/python/strict_cox_rerun.py --workers 12
python code/python/strict_qba_reanalysis.py --workers 4
python code/python/recompute_rmst_standard_ato.py
```

After extracting the replicate archives, verify the frozen evidence:

```bash
python code/python/verify_frozen_evidence.py
```

## R workflows

To regenerate only the approved Figure 1 without running any statistical analysis, from the repository root:

```bash
Rscript -e 'source("code/R/make_figure1.R"); make_figure1(getwd())'
```

Install the packages listed in `environment/R_packages.txt`, then run:

```bash
Rscript code/R/make_submission_figures.R .
Rscript code/R/verify_cox_vs_r_frozen.R .
python code/python/reverify_full_cox_sandwich.py
```

The second and third commands require `cox_validation_source_1200_datasets.zip` to have been extracted.

## Restricted SEER data

No patient-level SEER record is included or redistributable. Eligible users must obtain the same SEER Research Data release under their own data-use agreement. Once the two local case-list CSV files are available, the aggregate-only analysis is run as:

```bash
Rscript code/R/seer_descriptive_sensitivity.R <crc.csv> <bladder.csv> analysis/seer_descriptive_sensitivity
```

The script writes aggregate results only. See `docs/SEER_DATA_ACCESS.md` before running it.

## Integrity

`MANIFEST.csv` and `checksums/SHA256SUMS.txt` identify every public-package file. Separate SHA-256 checksums for the upload archives are stored beside them in the Zenodo payload.

## Licensing

Code is licensed under MIT. Original synthetic data, aggregate outputs, figures and documentation are licensed under CC BY 4.0. See `LICENSES.md` for scope. Restricted SEER data remain governed by the NCI data-use agreement.

## Citation and access

See `CITATION.cff` for the confirmed authors and identifiers. The repository is https://github.com/wangjinglin0905-gif/SEER-summary-code-bias-methods. Any Zenodo DOI will be added only after a record is created and verified.
