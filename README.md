# Predicted nodal risk and surgical extent in nonmetastatic appendiceal adenocarcinoma

Analysis code for:

> Li J, Wu J, Kuo ZC, Cai N, Liu Z. Predicted Nodal Risk and Cancer-Specific Mortality After Limited
> Resection or Oncologic Colectomy for Nonmetastatic Appendiceal Adenocarcinoma. *Submitted*, 2026.

The study links an individualized nodal-metastasis risk model to absolute 5-year cancer-specific
mortality differences between limited resection and oncologic colectomy, using risk-stratum-specific
propensity-score overlap weighting, competing-risk estimation, and full-pipeline bootstrap inference
against a prespecified 3-percentage-point decision bound.

## Data

This repository contains **code only**. It does not contain patient-level data.

The analysis uses the Surveillance, Epidemiology, and End Results (SEER) Program, SEER 17 Registries
Research Data, November 2025 submission. SEER research data are released by the US National Cancer
Institute to investigators who sign the SEER data-use agreement and cannot be redistributed by us.
Eligible investigators can request the same extract from <https://seer.cancer.gov/data/>.

To reproduce the analysis, export the appendiceal cohort from SEER*Stat with the variables listed in
`DATA_DICTIONARY.md`, place the file in `data/`, and run the scripts in the order given below.

## Reproducing the analysis

```r
Sys.setenv(A6_ROOT = "/path/to/your/project")   # scripts read this; default is the working directory
source("R/01B_risk_stratum_overlap_weighting_PRIMARY.R")
```

Directory layout the scripts expect:

```
<A6_ROOT>/
  02_data/      analytic_cohort_A.csv, Aline_frozen_pipeline_v2.RData
  03_results/   written by the scripts
  04_figures/   written by the figure scripts
```

### Order

| Script | Produces |
|---|---|
| `R/01A_no_summary_stage_primary_balance.R` | corrected propensity-score specification |
| `R/01B_risk_stratum_overlap_weighting_PRIMARY.R` | **primary analysis**: Table 2, Figure 2 |
| `R/02A_multiple_imputation_m20.R`, `R/02B_..._bootstrap_rubin_B1000.R` | Table S4 |
| `R/03_followup_support_calendar_eligible.R` | Table S5 |
| `R/04_continuous_risk_difference_B1000.R` | Figure S3 |
| `R/05_adequate_nodal_assessment_sensitivity.R` | Table S8 |
| `R/06_alternative_surgery_definitions_B1000.R` | Table S7 |
| `R/07_landmark_and_histology_B1000.R` | Tables S9, S10, S12-S14 |
| `R/08_global_consistency_audit.R` | cross-checks every reported number |
| `R/T1_*`, `R/T2_*`, `R/S*_*`, `R/F*_*` | the corresponding tables and figures |

Every inferential script re-estimates the whole pipeline inside each bootstrap replicate: the
nodal-risk model, individual predictions, percentile cut points, stratum-specific propensity scores,
overlap weights, and the outcome contrast. Runtime for the B = 1000 scripts is on the order of hours.

`R/01B_...` and `R/07_...` begin with a fingerprint gate that reproduces the frozen primary estimates
and stops if they do not match to 1e-9, so a broken environment fails loudly rather than silently.

## Reference results

`results_reference/` holds the summary outputs (CSV and text reports) from the frozen run reported in
the manuscript, so that a re-run can be compared line by line. No patient-level records are included.

## Environment

R 4.5.2. Packages: survival, splines, mice, survey, ggplot2 (plus officer and flextable for table
export). Exact versions are in `sessionInfo.txt`.

## License

Code: MIT (see `LICENSE`). The SEER data are governed by the SEER data-use agreement, not by this license.
