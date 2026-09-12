## =============================================================================
## A6 - SUPPLEMENTARY TABLE S8 SOURCE AUDIT ONLY
## Diagnosis-era x predicted nodal-risk-half sensitivity analysis
##
## Root:
##   <project root>
##
## INPUT:
##   03_results/06_negative_control_histology_era_FINAL_B1000.csv
##
## OUTPUT:
##   03_results/S8_era_source_audit.txt
##
## PURPOSE:
##   Audit the existing final B=1000 era x risk-half results before
##   constructing Supplementary Table S8.
##
## IMPORTANT:
##   This script does NOT re-run bootstrap and does NOT alter estimates.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
RESULT_DIR <- file.path(ROOT, "03_results")

dir.create(
    RESULT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

INPUT_FILE <- file.path(
    RESULT_DIR,
    "06_negative_control_histology_era_FINAL_B1000.csv"
)

if (!file.exists(INPUT_FILE)) {
    INPUT_FILE <- file.path(
        ROOT,
        "06_negative_control_histology_era_FINAL_B1000.csv"
    )
}

if (!file.exists(INPUT_FILE)) {
    stop(
        "找不到 06_negative_control_histology_era_FINAL_B1000.csv"
    )
}

d <- read.csv(
    INPUT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

needed <- c(
    "analysis",
    "estimate",
    "lo",
    "hi",
    "effective_B"
)

missing_cols <- setdiff(
    needed,
    names(d)
)

if (length(missing_cols) > 0) {
    stop(
        "CSV缺少字段：",
        paste(missing_cols, collapse = ", ")
    )
}

## Expected era rows from final Step 06:
expected_era_ids <- c(
    "2004-2009_low",
    "2004-2009_high",
    "2010-2017_low",
    "2010-2017_high",
    "2018-2023_low",
    "2018-2023_high"
)

era_rows <- d[
    d$analysis %in% expected_era_ids,
    ,
    drop = FALSE
]

if (nrow(era_rows) != 6) {
    stop(
        paste0(
            "Expected 6 era rows, found ",
            nrow(era_rows),
            ". Available analysis labels:\n",
            paste(d$analysis, collapse = " | ")
        )
    )
}

era_rows <- era_rows[
    match(
        expected_era_ids,
        era_rows$analysis
    ),
    ,
    drop = FALSE
]

## Basic validity checks.
if (any(!is.finite(era_rows$estimate))) {
    stop("At least one era point estimate is non-finite.")
}

if (any(!is.finite(era_rows$lo)) || any(!is.finite(era_rows$hi))) {
    stop("At least one era confidence interval is non-finite.")
}

if (any(era_rows$effective_B < 900, na.rm = TRUE)) {
    warning(
        "At least one era estimate has effective bootstrap B <900."
    )
}

## Split labels.
era_rows$era <- sub(
    "_(low|high)$",
    "",
    era_rows$analysis
)

era_rows$risk_half <- ifelse(
    grepl("_low$", era_rows$analysis),
    "Below median risk",
    "At/above median risk"
)

## Convert to percentage points for easier audit.
era_rows$estimate_pp <- 100 * era_rows$estimate
era_rows$lo_pp <- 100 * era_rows$lo
era_rows$hi_pp <- 100 * era_rows$hi

## Within-era high-minus-low contrast in point estimates only.
contrast_rows <- data.frame(
    era = c(
        "2004-2009",
        "2010-2017",
        "2018-2023"
    ),
    low_pp = NA_real_,
    high_pp = NA_real_,
    high_minus_low_pp = NA_real_,
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(contrast_rows))) {

    ee <- contrast_rows$era[i]

    low_val <- era_rows$estimate_pp[
        era_rows$era == ee &
        era_rows$risk_half == "Below median risk"
    ]

    high_val <- era_rows$estimate_pp[
        era_rows$era == ee &
        era_rows$risk_half == "At/above median risk"
    ]

    contrast_rows$low_pp[i] <- low_val
    contrast_rows$high_pp[i] <- high_val
    contrast_rows$high_minus_low_pp[i] <- high_val - low_val
}

OUT_FILE <- file.path(
    RESULT_DIR,
    "S8_era_source_audit.txt"
)

report <- c(
    "== S8 diagnosis-era x predicted nodal-risk-half source audit ==",
    paste0("Source: ", INPUT_FILE),
    "",
    "Interpretation:",
    "RD = oncologic colectomy - limited resection in 5-year cancer-specific mortality.",
    "Negative RD indicates lower cancer-specific mortality associated with oncologic colectomy.",
    "These existing era analyses used the cohort-wide overlap weights from the primary model.",
    "The overall median predicted nodal risk defined low vs high risk within each bootstrap replicate.",
    "",
    "== Final era x risk-half B=1000 results =="
)

for (i in seq_len(nrow(era_rows))) {

    report <- c(
        report,
        sprintf(
            "%-18s | %-20s | RD %+.1f pp (95%% CI %+.1f to %+.1f) | effective B=%d",
            era_rows$era[i],
            era_rows$risk_half[i],
            era_rows$estimate_pp[i],
            era_rows$lo_pp[i],
            era_rows$hi_pp[i],
            era_rows$effective_B[i]
        )
    )
}

report <- c(
    report,
    "",
    "== Point-estimate high-minus-low contrasts (descriptive only) =="
)

for (i in seq_len(nrow(contrast_rows))) {

    report <- c(
        report,
        sprintf(
            "%-10s | low %+.1f pp | high %+.1f pp | high-low %+.1f pp",
            contrast_rows$era[i],
            contrast_rows$low_pp[i],
            contrast_rows$high_pp[i],
            contrast_rows$high_minus_low_pp[i]
        )
    )
}

report <- c(
    report,
    "",
    "Decision rule before building S8:",
    "1. Confirm all six era rows are present and effective B is adequate.",
    "2. Confirm the direction/pattern is clinically interpretable.",
    "3. Do not claim a formal era interaction unless a valid interaction test was prespecified/performed.",
    "4. Treat this as a robustness/temporal-consistency sensitivity analysis."
)

writeLines(
    report,
    OUT_FILE,
    useBytes = TRUE
)

cat(
    paste(report, collapse = "\n"),
    "\n"
)

cat("\n============================================================\n")
cat("S8 ERA SOURCE AUDIT COMPLETE\n")
cat("============================================================\n")
cat("Saved: ", OUT_FILE, "\n", sep = "")
cat("No analysis was re-run.\n")
cat("============================================================\n")
