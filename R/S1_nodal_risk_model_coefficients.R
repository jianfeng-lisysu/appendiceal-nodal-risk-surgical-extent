## =============================================================================
## A6 - SUPPLEMENTARY TABLE S1 ONLY
## Frozen nodal-metastasis risk model coefficients
## JOURNAL-READY THREE-LINE TABLE — FINAL LOCK VERSION
##
## Root:
##   <project root>
##
## INPUTS
##   02_data/analytic_cohort_A.csv
##   02_data/analytic_cohort_A_riskscored.csv
##   02_data/Aline_frozen_pipeline_v2.RData
##
## OUTPUTS
##   05_tables/Supplementary_Table_S1_Nodal_Risk_Model_FINAL.docx
##   05_tables/Supplementary_Table_S1_Nodal_Risk_Model_FINAL.csv
##
## ANALYSIS
##   Development cohort:
##     oncologic colectomy + >=12 examined lymph nodes + known nodal status
##   Frozen expected n = 3,252
##   Frozen expected node-positive n = 738
##
## IMPORTANT
##   This script reproduces the already-frozen risk model.
##   It does NOT change the formula, cohort, coding, or estimand.
## =============================================================================


## =============================================================================
## 0. Setup
## =============================================================================

rm(list = ls())

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory

DATA_DIR <- file.path(
    ROOT,
    "02_data"
)

TABLE_DIR <- file.path(
    ROOT,
    "05_tables"
)

dir.create(
    TABLE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


## =============================================================================
## 1. Packages
## =============================================================================

need_pkg <- c(
    "officer",
    "flextable"
)

miss_pkg <- need_pkg[
    !vapply(
        need_pkg,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (
    length(
        miss_pkg
    ) >
    0
) {

    stop(
        "请先安装：",
        paste(
            miss_pkg,
            collapse = ", "
        )
    )
}

library(
    officer
)

library(
    flextable
)


## =============================================================================
## 2. File locator
## =============================================================================

locate_file <- function(
    filename
) {

    p1 <- file.path(
        DATA_DIR,
        filename
    )

    p2 <- file.path(
        ROOT,
        filename
    )

    if (
        file.exists(
            p1
        )
    ) {
        return(
            p1
        )
    }

    if (
        file.exists(
            p2
        )
    ) {
        return(
            p2
        )
    }

    stop(
        "找不到文件：",
        filename
    )
}


RAW_FILE <- locate_file(
    "analytic_cohort_A.csv"
)

RISK_FILE <- locate_file(
    "analytic_cohort_A_riskscored.csv"
)

PIPE_FILE <- locate_file(
    "Aline_frozen_pipeline_v2.RData"
)


## =============================================================================
## 3. Load frozen objects
## =============================================================================

load(
    PIPE_FILE
)

required_objects <- c(
    "prep",
    "RISK_FORMULA"
)

missing_objects <- required_objects[
    !vapply(
        required_objects,
        exists,
        logical(1),
        inherits = TRUE
    )
]

if (
    length(
        missing_objects
    ) >
    0
) {

    stop(
        "Frozen RData missing: ",
        paste(
            missing_objects,
            collapse = ", "
        )
    )
}


RAW <- read.csv(
    RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

RS <- read.csv(
    RISK_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


## =============================================================================
## 4. Frozen cohort fingerprints
## =============================================================================

if (
    nrow(
        RAW
    ) !=
    6951
) {

    stop(
        "Frozen analytic cohort N should be 6951; current N=",
        nrow(
            RAW
        )
    )
}

if (
    nrow(
        RS
    ) !=
    6951
) {

    stop(
        "Frozen risk-scored cohort N should be 6951; current N=",
        nrow(
            RS
        )
    )
}


## =============================================================================
## 5. Recreate frozen prepared data
## =============================================================================

d <- prep(
    RAW
)


## =============================================================================
## 6. Frozen development cohort
## =============================================================================

dev_mask <- (
    d$exposure ==
    "colectomy" &
    !is.na(
        d$ln_examined
    ) &
    d$ln_examined >=
    12 &
    !is.na(
        d$n_positive
    )
)

dev <- d[
    dev_mask,
    ,
    drop = FALSE
]


if (
    nrow(
        dev
    ) !=
    3252
) {

    stop(
        "Development cohort fingerprint failed: expected n=3252; current n=",
        nrow(
            dev
        )
    )
}


## Determine node-positive outcome directly from nodal count.
NODE_POSITIVE_N <- sum(
    dev$n_positive >
    0,
    na.rm = TRUE
)

if (
    NODE_POSITIVE_N !=
    738
) {

    stop(
        "Node-positive fingerprint failed: expected 738; current n=",
        NODE_POSITIVE_N
    )
}


## =============================================================================
## 7. Refit the frozen model
## =============================================================================

fit <- glm(
    RISK_FORMULA,
    data = dev,
    family = binomial()
)


## =============================================================================
## 8. Full-cohort risk-score fingerprint
## =============================================================================

pred_full <- predict(
    fit,
    newdata = d,
    type = "response"
)

if (
    any(
        !is.finite(
            pred_full
        )
    )
) {

    stop(
        "Non-finite full-cohort predicted risks."
    )
}


if (
    "patient_id" %in%
    names(
        RS
    ) &&
    "patient_id" %in%
    names(
        d
    ) &&
    identical(
        as.character(
            RS$patient_id
        ),
        as.character(
            d$patient_id
        )
    )
) {

    risk_fp <- max(
        abs(
            pred_full -
            RS$nodal_risk
        ),
        na.rm = TRUE
    )

} else {

    ## Cohort ordering was previously frozen; if patient ID serialization differs,
    ## use the row-order risk vector fingerprint.
    if (
        length(
            pred_full
        ) !=
        length(
            RS$nodal_risk
        )
    ) {

        stop(
            "Risk fingerprint length mismatch."
        )
    }

    risk_fp <- max(
        abs(
            pred_full -
            RS$nodal_risk
        ),
        na.rm = TRUE
    )
}


if (
    risk_fp >
    1e-8
) {

    stop(
        sprintf(
            "Frozen risk-score fingerprint failed: max error = %.12g",
            risk_fp
        )
    )
}


## =============================================================================
## 9. Extract model coefficients
## =============================================================================

coef_matrix <- summary(
    fit
)$coefficients

coef_df <- data.frame(
    term = rownames(
        coef_matrix
    ),
    beta = coef_matrix[
        ,
        1
    ],
    se = coef_matrix[
        ,
        2
    ],
    z = coef_matrix[
        ,
        3
    ],
    p = coef_matrix[
        ,
        4
    ],
    stringsAsFactors = FALSE,
    row.names = NULL
)

coef_df$ci_beta_lo <- coef_df$beta -
    qnorm(
        0.975
    ) *
    coef_df$se

coef_df$ci_beta_hi <- coef_df$beta +
    qnorm(
        0.975
    ) *
    coef_df$se

coef_df$or <- exp(
    coef_df$beta
)

coef_df$or_lo <- exp(
    coef_df$ci_beta_lo
)

coef_df$or_hi <- exp(
    coef_df$ci_beta_hi
)


## =============================================================================
## 10. Publication labels
## =============================================================================

label_map <- c(
    "(Intercept)" =
        "Intercept",
    "T_fT2" =
        "T2 vs T1",
    "T_fT3" =
        "T3 vs T1",
    "T_fT4" =
        "T4 vs T1",
    "grade_f1" =
        "Grade 1 vs unknown",
    "grade_f2" =
        "Grade 2 vs unknown",
    "grade_f3" =
        "Grade 3 vs unknown",
    "grade_f4" =
        "Grade 4 vs unknown",
    "hist_fSRCC" =
        "Signet-ring cell carcinoma vs GCA",
    "hist_fmucinous" =
        "Mucinous adenocarcinoma vs GCA",
    "hist_fnonmucinous" =
        "Nonmucinous adenocarcinoma vs GCA",
    "logsize" =
        "Log-transformed tumor size",
    "size_miss" =
        "Tumor size missing",
    "age" =
        "Age, years",
    "I(age^2)" =
        "Age squared",
    "female" =
        "Female sex",
    "year" =
        "Year of diagnosis"
)


coef_df$Predictor <- ifelse(
    coef_df$term %in%
    names(
        label_map
    ),
    unname(
        label_map[
            coef_df$term
        ]
    ),
    coef_df$term
)


## =============================================================================
## 11. Preferred display order
## =============================================================================

preferred_order <- c(
    "age",
    "I(age^2)",
    "female",
    "year",
    "logsize",
    "size_miss",
    "T_fT2",
    "T_fT3",
    "T_fT4",
    "grade_f1",
    "grade_f2",
    "grade_f3",
    "grade_f4",
    "hist_fnonmucinous",
    "hist_fmucinous",
    "hist_fSRCC",
    "(Intercept)"
)


coef_df$order_id <- match(
    coef_df$term,
    preferred_order
)

coef_df$order_id[
    is.na(
        coef_df$order_id
    )
] <- length(
    preferred_order
) +
    seq_len(
        sum(
            is.na(
                coef_df$order_id
            )
        )
    )

coef_df <- coef_df[
    order(
        coef_df$order_id
    ),
    ,
    drop = FALSE
]


## =============================================================================
## 12. Formatting helpers
## =============================================================================

fmt_beta <- function(
    x
) {

    sprintf(
        "%.4f",
        x
    )
}


fmt_se <- function(
    x
) {

    sprintf(
        "%.4f",
        x
    )
}


fmt_or <- function(
    x
) {

    if (
        !is.finite(
            x
        )
    ) {
        return(
            "—"
        )
    }

    if (
        x >=
        1000 ||
        x <
        0.001
    ) {

        return(
            format(
                x,
                scientific = TRUE,
                digits = 3
            )
        )
    }

    sprintf(
        "%.3f",
        x
    )
}


fmt_or_ci <- function(
    lo,
    hi
) {

    if (
        !is.finite(
            lo
        ) ||
        !is.finite(
            hi
        )
    ) {
        return(
            "—"
        )
    }

    paste0(
        fmt_or(
            lo
        ),
        " to ",
        fmt_or(
            hi
        )
    )
}


fmt_p <- function(
    x
) {

    if (
        is.na(
            x
        )
    ) {
        return(
            "—"
        )
    }

    if (
        x <
        0.001
    ) {
        return(
            "<0.001"
        )
    }

    sprintf(
        "%.3f",
        x
    )
}


## =============================================================================
## 13. Final display table
## =============================================================================

table_s1 <- data.frame(
    Predictor =
        coef_df$Predictor,
    `β coefficient` =
        vapply(
            coef_df$beta,
            fmt_beta,
            character(
                1
            )
        ),
    SE =
        vapply(
            coef_df$se,
            fmt_se,
            character(
                1
            )
        ),
    `Odds ratio` =
        vapply(
            coef_df$or,
            fmt_or,
            character(
                1
            )
        ),
    `95% CI for odds ratio` =
        mapply(
            fmt_or_ci,
            coef_df$or_lo,
            coef_df$or_hi,
            USE.NAMES = FALSE
        ),
    `P value` =
        vapply(
            coef_df$p,
            fmt_p,
            character(
                1
            )
        ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## Intercept OR is not clinically interpretable; retain beta for reproducibility.
intercept_row <- which(
    coef_df$term ==
    "(Intercept)"
)

if (
    length(
        intercept_row
    ) ==
    1
) {

    table_s1[
        intercept_row,
        "Odds ratio"
    ] <- "—"

    table_s1[
        intercept_row,
        "95% CI for odds ratio"
    ] <- "—"
}


## =============================================================================
## 14. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S1_Nodal_Risk_Model_FINAL.csv"
)

write.csv(
    table_s1,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 15. Build journal-ready three-line table
## =============================================================================

ft <- flextable(
    table_s1
)

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.3,
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.6,
    part = "header"
)

ft <- bold(
    ft,
    part = "header"
)

ft <- align(
    ft,
    j = 1,
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 2:6,
    align = "center",
    part = "all"
)

ft <- valign(
    ft,
    valign = "center",
    part = "all"
)


## -----------------------------------------------------------------------------
## Three-line table: remove all borders, then restore top/header/bottom rules.
## -----------------------------------------------------------------------------

ft <- border_remove(
    ft
)

top_border <- fp_border(
    color = "black",
    width = 1.2
)

mid_border <- fp_border(
    color = "black",
    width = 0.8
)

bottom_border <- fp_border(
    color = "black",
    width = 1.2
)

ft <- hline_top(
    ft,
    part = "header",
    border = top_border
)

ft <- hline(
    ft,
    i = 1,
    part = "header",
    border = mid_border
)

ft <- hline_bottom(
    ft,
    part = "body",
    border = bottom_border
)


## -----------------------------------------------------------------------------
## Compact spacing and widths.
## -----------------------------------------------------------------------------

ft <- padding(
    ft,
    padding.top = 1.2,
    padding.bottom = 1.2,
    padding.left = 1.7,
    padding.right = 1.7,
    part = "all"
)

ft <- width(
    ft,
    j = 1,
    width = 2.45
)

ft <- width(
    ft,
    j = 2,
    width = 1.00
)

ft <- width(
    ft,
    j = 3,
    width = 0.80
)

ft <- width(
    ft,
    j = 4,
    width = 1.00
)

ft <- width(
    ft,
    j = 5,
    width = 1.70
)

ft <- width(
    ft,
    j = 6,
    width = 0.85
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 16. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S1_Nodal_Risk_Model_FINAL.docx"
)

doc <- read_docx()


## Title
doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S1. Coefficients of the final nodal-metastasis risk model",
            fp_text(
                font.family = "Times New Roman",
                font.size = 9.5,
                bold = TRUE
            )
        )
    )
)


## Table
doc <- body_add_flextable(
    doc,
    value = ft,
    align = "center"
)


## Footnote
note_text <- paste0(
    "The logistic model was developed among patients undergoing oncologic colectomy with ",
    "at least 12 examined lymph nodes and known nodal status (n = 3,252; node-positive n = 738). ",
    "T1, unknown grade, goblet cell adenocarcinoma (GCA), and male sex are the reference categories ",
    "for the corresponding categorical predictors. ",
    "Age and age squared were entered jointly and should not be interpreted as separate marginal effects. ",
    "The intercept odds ratio is omitted because it is not clinically interpretable; ",
    "the beta coefficient is retained for model reproducibility. Odds ratios and P values are descriptive ",
    "model coefficients and are not intended as causal effect estimates."
)

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            note_text,
            fp_text(
                font.family = "Times New Roman",
                font.size = 7.6
            )
        )
    )
)


## Landscape section for stable single-page presentation.
doc <- body_end_section_landscape(
    doc
)


print(
    doc,
    target = DOCX_FILE
)


## =============================================================================
## 17. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S1 FINAL LOCK VERSION COMPLETE\n")
cat("============================================================\n")

cat(
    "Development cohort n = ",
    nrow(
        dev
    ),
    "\n",
    sep = ""
)

cat(
    "Node-positive n = ",
    NODE_POSITIVE_N,
    "\n",
    sep = ""
)

cat(
    sprintf(
        "Frozen risk-score max error = %.12g\n",
        risk_fp
    )
)

cat(
    "\nFrozen formula:\n"
)

print(
    RISK_FORMULA
)

cat(
    "\nModel coefficients:\n"
)

print(
    coef_df[
        ,
        c(
            "term",
            "beta",
            "se",
            "p"
        )
    ],
    row.names = FALSE
)


cat("
WORDING AUDIT:
")
cat("Title: Supplementary Table S1. Coefficients of the final nodal-metastasis risk model
")
cat("Age note: Age and age squared were entered jointly and should not be interpreted as separate marginal effects.
")

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Standard three-line table; no vertical rules.\n")
cat("No frozen model component was changed.\n")
cat("============================================================\n")
