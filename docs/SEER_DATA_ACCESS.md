# SEER data access and redistribution boundary

The manuscript uses restricted-access SEER Research Data, 17 Registries, November 2025 submission (2000–2023), with diagnoses restricted to 2004–2022. Access requires an independent application and execution of the applicable NCI SEER data-use agreement.

## Included here

- cohort-flow counts;
- aggregate baseline summaries;
- endpoint estimates and confidence intervals;
- effective sample sizes and propensity-score diagnostics;
- unknown/missing-category counts;
- aggregate cumulative incidence, event counts, surgery-code composition, and piecewise estimates.

## Not included

- case-list exports or patient-level records;
- row-level identifiers or dates;
- any file from which a case list can be reconstructed;
- the local raw-data paths used by the authors.

## Recalculation

Authorized users should provide their own CRC and bladder case-list CSV paths to `code/R/seer_descriptive_sensitivity.R`. The script writes only aggregate files to the requested output directory. Review every output against the current SEER data-use agreement before sharing.
