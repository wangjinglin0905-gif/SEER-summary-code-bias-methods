# Frozen-evidence verification report

- Checks executed: 1,499
- Checks passed: 1,499
- Checks failed: 0
- Primary retained datasets: 252,000
- H2 retained datasets: 14,000
- QBA evaluation rows: 180,000 across four priors
- pRMST scenario rows: 12

## Key rebuilt results

- Maximum primary score-residual sum: 3.379e-09
- Core pRMST coverage range: 0.750–0.956
- QBA pooled means:

| Prior | Bias (log HR) | RMSE (log HR) | Coverage | Width ratio |
|---|---:|---:|---:|---:|
| global_mechanism_prior | 0.023 | 0.077 | 0.737 | 1.202 |
| oracle_cell_specific | -0.000 | 0.041 | 0.952 | 1.003 |
| oracle_mean_x1.25 | 0.009 | 0.045 | 0.926 | 1.003 |
| oracle_mean_x1.50 | 0.018 | 0.054 | 0.852 | 1.003 |

## Failed checks

None.

## Scope

This audit independently re-aggregated frozen synthetic replicate-level outputs and checked frozen SEER aggregate tables. It did not re-read patient-level SEER files because the external G:/E: data volumes were not mounted during this audit.
