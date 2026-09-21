# Reproducibility tiers

## Tier 1: immediate, repository-only

The frozen summary tables, figure source data, figures, verification ledgers, and corrected pRMST truths can be inspected without restricted data or large archives.

## Tier 2: synthetic replicate verification

Extract `synthetic_primary_h2_replicates.zip` and `qba_replicate_and_calibration_outputs.zip` from the Zenodo payload at repository root, then run `python code/python/verify_frozen_evidence.py`. The expected result is 1,499 passed checks and zero failures.

## Tier 3: full strict regeneration

Run the strict Cox/H2 and QBA scripts. These jobs are intentionally expensive; deterministic seed schemes and the executed scenario registry are retained. Regenerated outputs should be compared by checksum or numerical tolerance, not by file timestamp.

## Tier 4: R/Python Cox cross-implementation

Extract `cox_validation_source_1200_datasets.zip`, run the R script, then the Python comparison. Expected tolerances are documented in the frozen verification summary.

## Tier 5: restricted SEER reconstruction

This tier requires an independent SEER data-use agreement and is not possible from public files alone. Only aggregate outputs may be moved into the public package.
