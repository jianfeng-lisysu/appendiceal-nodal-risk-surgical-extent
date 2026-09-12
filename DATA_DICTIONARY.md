# Data dictionary

The analysis expects `02_data/analytic_cohort_A.csv`, built from a SEER*Stat case listing of the
SEER 17 Registries Research Data (November 2025 submission) by `R/00_build_cohort.R` (see note below).
One row per patient; 25 columns.

## Cohort construction (frozen)

Starting from 24,632 malignant appendiceal tumours:

| Step | Restriction | Remaining |
|---|---|---|
| E1 | diagnosis 2004-2023 | 23,021 |
| E2 | prespecified epithelial histologies | 13,972 |
| E3 | microscopically confirmed | 13,877 |
| E4 | not death-certificate-only or autopsy-only | 13,868 |
| E5 | localized or regional Summary Stage | 8,223 |
| E6 | harmonized T1-T4 available | 7,689 |
| E7 | surgery code 30-32 or 40-41 | 6,951 |
| E8 | survival information available | **6,951** |

Final cohort: limited resection 2,877; oncologic colectomy 4,074; 6,906 evaluable for the
competing-risk outcome.

## Variables

| Column | Definition |
|---|---|
| `patient_id` | SEER patient identifier |
| `year`, `era` | year of diagnosis; era 2004-2009 / 2010-2017 / 2018-2023 |
| `age`, `female` | age at diagnosis; 1 = female |
| `race_eth`, `income`, `rural_urban`, `marital` | SEER recodes used in the propensity-score model |
| `hist_group` | nonmucinous (8140, 8263, 8261, 8262, 8255, 8010, 8144, 8210, 8211); mucinous (8480, 8481, 8470); GCA (8243, 8244, 8245, 8249); SRCC (8490) |
| `grade` | 1-4; Grade Recode through 2017 and Derived Summary Grade 2018+ reconciled; missing kept as its own category |
| `T` | harmonized T1-T4: EOD 2018 T recode, then Derived SEER Combined T (2016-2017), then AJCC 7th, then AJCC 6th |
| `size_mm` | CS tumor size (2004-2015) then Tumor Size Summary (2016+), valid 1-988 mm |
| `stage_sum` | localized or regional |
| `surg_code`, `exposure` | surgery of primary site; limited resection = 30-32, oncologic colectomy = 40-41. 2023 alphanumeric codes (A300 etc.) decoded to their numeric equivalents |
| `ln_examined`, `n_positive`, `adequate_ln` | regional nodes examined (<90); 1 if any positive node, 0 if none examined positive, missing if not assessed; adequate = at least 12 examined |
| `cea`, `chemo` | CEA interpretation (2010+); chemotherapy recode |
| `os_months`, `os_event` | survival months; 1 = dead from any cause |
| `css_event` | 1 = death attributable to this cancer; 0 = alive or dead of other cause; missing when cause of death is unknown (excluded from cancer-specific analyses, retained for overall survival) |
| `other_death` | competing event: died of another cause |

## Note on the cohort-building script

`R/00_build_cohort.R` is not included because it is written against a specific SEER*Stat export layout
and contains the raw column headers of that export. The rules above are complete: any equivalent export
containing the listed SEER variables reproduces the cohort. The frozen cohort is identified by
N = 6,951 with 1,212 cancer-specific and 2,111 all-cause deaths; `R/01B_...` and `R/07_...` will stop if
these are not reproduced.
