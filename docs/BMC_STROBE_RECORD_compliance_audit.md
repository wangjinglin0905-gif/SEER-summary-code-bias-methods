# BMC MRM, STROBE, and RECORD compliance audit

Audit date: 21 August 2026  
Target: BMC Medical Research Methodology, Research article  
Primary reporting framework: ADEMP for the simulation study, with STROBE and RECORD applied to the two routinely collected SEER illustrations.

## Journal-level requirements

| Requirement | Evidence checked | Implementation in revised manuscript | Status |
|---|---|---|---|
| In-scope methodological research | BMC MRM states that it publishes original methodological research and empirical studies linking methodological choice to study outcomes | The title and abstract identify an ADEMP simulation study; SEER analyses are explicitly illustrations | Met |
| Structured abstract, maximum 350 words | BMC research-article guidance | Background, Methods, Results, Conclusions; no citations; target 300-330 words | Met in draft |
| Three to ten keywords | BMC research-article guidance | Six keywords | Met in draft |
| Required article structure | BMC research-article guidance | Background, Methods, Results, Discussion, Conclusions, abbreviations, Declarations | Met in draft |
| Complete Declarations headings | BMC research-article guidance | All required headings retained, including Acknowledgements | Met in draft |
| Ethics statement for human data, including committee/exemption detail | BMC guidance requires the committee name and reference when applicable, including exemptions | A clearly flagged author-confirmation field is retained; no exemption is invented | **Author action required** |
| Data availability with access conditions and persistent location | BMC research-article guidance | SEER access restrictions are stated; code/source-data repository remains flagged until a persistent URL/DOI exists | **Author action required** |
| Transparent AI use and human accountability | BMC editorial policy treats interpretive AI use as permissible only with human oversight, verification, confidentiality protection, and disclosure | Methods and Declarations state the uses, verification boundary, and that no patient-level SEER record was sent to an external AI service; exact tool/model wording remains for author confirmation | **Author action required** |
| Potential duplicate/overlapping publication | BMC requires transparent declaration and citation of potentially overlapping work | Companion-manuscript status and degree of overlap are flagged for confirmation and editor disclosure | **Author action required; submission blocker** |

Official sources:

- BMC Medical Research Methodology research-article instructions: https://link.springer.com/journal/12874/submission-guidelines/research-article
- Journal aims and scope: https://link.springer.com/journal/12874/submission-guidelines
- BMC editorial policies, including ethics, AI, citations, and duplicate publication: https://link.springer.com/brands/bmc/editorial-policies

## STROBE cohort-study items

| Item | Required reporting | Revised location | Assessment |
|---|---|---|---|
| 1 | Design in title/abstract and balanced abstract | Title; Abstract | Met; SEER component described as a retrospective registry illustration |
| 2-3 | Rationale and objectives | Background, final two paragraphs | Met |
| 4-5 | Design, setting, dates | Methods: Study architecture; SEER data source and setting | Met; November 2025 SEER release and diagnosis years 2004-2022 stated |
| 6 | Eligibility, selection, and follow-up | Methods: SEER cohorts; Figure 6a flow | Met; sequential counts are reported |
| 7-8 | Variables and measurement | Methods: Exposure, outcomes, covariates | Met; observed surgery code is distinguished from latent strategy; period-specific stage variables described |
| 9 | Bias | Methods and Discussion | Met; exposure misclassification, residual confounding, cause-of-death error, missing-category limitations, and case-informed calibration stated |
| 10 | Study size | Methods | Met; census of eligible records in the requested SEER extracts, not a power-selected cohort |
| 11 | Quantitative variables | Methods; Supplementary table | Met; spline terms and prespecified categories stated |
| 12a | Statistical methods/confounding control | Methods: SEER statistical analysis | Met |
| 12b | Subgroups/interactions | Methods and Results | Not applicable; no confirmatory subgroup effect estimates |
| 12c | Missing data | Methods; Results; Table 4 | Met; unknown-as-category primary analysis and complete-size-and-grade restriction sensitivity reported; explicitly not called imputation |
| 12d | Loss to follow-up | Methods | Met; five-year censoring and reverse Kaplan-Meier follow-up stated; no unsupported informative-censoring claim |
| 12e | Sensitivity analyses | Methods; Results | Met; PH and known-size/grade analyses reported symmetrically for both cohorts |
| 13 | Participant flow | Results; Figure 6a | Met |
| 14 | Characteristics, missingness, follow-up | Results; Table 4; Supplementary baseline table | Met |
| 15 | Outcome events | Results; Table 4 | Met; cancer, other-cause, and all-death totals reported |
| 16 | Crude/adjusted estimates, precision, absolute risk | Results; Figure 6b-c; Table 5 | Met; overlap-weighted HRs, cross-outcome ratios, and weighted five-year CIFs reported; unweighted estimates are not labelled as causal comparators |
| 17 | Other analyses | Results | Met; piecewise PH sensitivity identified as post hoc |
| 18-21 | Summary, limitations, interpretation, generalisability | Discussion; Conclusions | Met; no transport of simulated bias magnitude to SEER |
| 22 | Funding | Declarations | Present but grant-number year conflict requires author confirmation |

## RECORD add-on items

| RECORD domain | Revised implementation | Assessment |
|---|---|---|
| Routinely collected data named in title/abstract | SEER is named in title and abstract | Met |
| Database name, release, coverage, and access | SEER Research Data, 17 Registries, November 2025 submission (2000-2023); analysis diagnosis years 2004-2022; DUA restriction stated | Met |
| Population-selection algorithms and codes | Sequential cohort criteria, ICD-O-3 histologies, site definitions, stage harmonisation, and surgery-code groups are reported; executable script and aggregate flow table are archived | Met |
| Validation of codes/algorithms | Aggregate cohort reconstruction and raw-file checksums reported; no external clinical validation of the treatment-code construct is claimed | Met with limitation |
| Data cleaning and linkage | No record linkage; deterministic recoding and range checks described; patient-level exports are not redistributed | Met |
| Flow and attrition | Counts at every selection step reported for both cancers | Met |
| Data access and reproducibility | SEER access conditions described; aggregate tables/code prepared | Persistent public repository still required |

Reporting-guideline sources:

- STROBE: https://www.equator-network.org/reporting-guidelines/strobe/
- RECORD: https://www.equator-network.org/reporting-guidelines/record/

## Submission decision from reporting audit

The revised scientific and reporting structure is suitable for journal-level review. It is not yet a clean submission package until four author-owned items are resolved: ethics determination wording, persistent repository link/DOI, exact funding numbers, and the current status/citation/editor disclosure for the overlapping Chinese-language manuscript.
