## =============================================================================
## A6 - SUPPLEMENTARY TABLE S7 ONLY
## Sensitivity analyses using alternative surgery definitions
## JOURNAL-READY THREE-LINE TABLE
##
## Root:
##   <project root>
##
## INPUT:
##   03_results/06_S2_alternative_surgery_FINAL_B1000.csv
##
## OUTPUT:
##   05_tables/Supplementary_Table_S7_Alternative_Surgery_Definitions_FINAL.docx
##   05_tables/Supplementary_Table_S7_Alternative_Surgery_Definitions_FINAL.csv
##
## FINAL ANALYSIS SOURCE
##   Each alternative surgery-definition analysis used full-pipeline
##   bootstrap B=1000:
##     resample alternative-surgery cohort
##     -> prep
##     -> refit nodal-risk model
##     -> recompute nodal risk
##     -> refit PS/OW
##     -> re-estimate risk quintiles
##     -> weighted 5-year cancer-specific mortality RD
##
## DEFINITIONS
##   A. Limited codes 20-39 vs colectomy codes 40-79
##   B. Exclude code 32: limited codes 30-31 vs colectomy codes 40-41
##   C. Limited codes 30-32 vs colectomy codes 40-79
##
## IMPORTANT
##   This script formats the already-completed final B=1000 results.
##   It does NOT reconstruct PRE-E7 data or re-run bootstrap.
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
RESULT_DIR <- file.path(ROOT, "03_results")
TABLE_DIR <- file.path(ROOT, "05_tables")

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

if (length(miss_pkg) > 0) {
    stop(
        "请先安装：",
        paste(
            miss_pkg,
            collapse = ", "
        )
    )
}

library(officer)
library(flextable)


## =============================================================================
## 2. Locate final B=1000 result file
## =============================================================================

INPUT_FILE <- file.path(
    RESULT_DIR,
    "06_S2_alternative_surgery_FINAL_B1000.csv"
)

if (!file.exists(INPUT_FILE)) {
    INPUT_FILE <- file.path(
        ROOT,
        "06_S2_alternative_surgery_FINAL_B1000.csv"
    )
}

if (!file.exists(INPUT_FILE)) {
    stop(
        paste0(
            "找不到最终替代手术定义结果：",
            "06_S2_alternative_surgery_FINAL_B1000.csv"
        )
    )
}


d <- read.csv(
    INPUT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


## =============================================================================
## 3. Required fields
## =============================================================================

needed <- c(
    "definition",
    "description",
    "n",
    "limited_n",
    "colectomy_n",
    "quintile",
    "rd",
    "rd_pp",
    "ci_lo",
    "ci_hi",
    "ci_lo_pp",
    "ci_hi_pp",
    "effective_B"
)

missing_cols <- setdiff(
    needed,
    names(d)
)

if (length(missing_cols) > 0) {
    stop(
        "最终 S2 CSV 缺少字段：",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


## =============================================================================
## 4. Freeze definitions and order
## =============================================================================

definition_order <- c(
    "S2_A_20_39_vs_40_79",
    "S2_B_exclude_code32",
    "S2_C_colectomy_40_79"
)

quintile_order <- paste0(
    "Q",
    1:5
)

if (!all(definition_order %in% d$definition)) {
    stop(
        "最终 CSV 缺少一个或多个替代手术定义。"
    )
}

if (nrow(d) != 15) {
    stop(
        "Expected 15 rows (3 definitions × 5 quintiles); current n=",
        nrow(d)
    )
}


d$definition_factor <- factor(
    d$definition,
    levels = definition_order
)

d$quintile_factor <- factor(
    d$quintile,
    levels = quintile_order
)

d <- d[
    order(
        d$definition_factor,
        d$quintile_factor
    ),
    ,
    drop = FALSE
]


## =============================================================================
## 5. Hard population fingerprints
## =============================================================================

pop_check <- unique(
    d[
        ,
        c(
            "definition",
            "n",
            "limited_n",
            "colectomy_n"
        )
    ]
)

expected_pop <- data.frame(
    definition = definition_order,
    n = c(
        7450,
        6676,
        7115
    ),
    limited_n = c(
        3212,
        2602,
        2877
    ),
    colectomy_n = c(
        4238,
        4074,
        4238
    ),
    stringsAsFactors = FALSE
)

pop_check <- pop_check[
    match(
        definition_order,
        pop_check$definition
    ),
    ,
    drop = FALSE
]

if (
    !identical(
        as.integer(pop_check$n),
        as.integer(expected_pop$n)
    ) ||
    !identical(
        as.integer(pop_check$limited_n),
        as.integer(expected_pop$limited_n)
    ) ||
    !identical(
        as.integer(pop_check$colectomy_n),
        as.integer(expected_pop$colectomy_n)
    )
) {
    stop(
        "Alternative-surgery population-size fingerprint failed."
    )
}


## =============================================================================
## 6. B=1000 fingerprint
## =============================================================================

if (
    any(
        d$effective_B != 1000,
        na.rm = TRUE
    )
) {
    stop(
        "At least one alternative-surgery estimate does not have effective B=1000."
    )
}


## =============================================================================
## 7. Final point-estimate fingerprints
##
## Use the actual final B=1000 results, not the older legacy point estimates.
## =============================================================================

expected_final_rd_pp <- c(
    ## A: limited 20-39 vs colectomy 40-79
    -0.4405943,
    -2.3862292,
    -0.9833170,
    -3.8303532,
    -6.8441614,

    ## B: limited 30-31 vs colectomy 40-41
    +0.2239056,
    -3.4173858,
    -0.8417034,
    -4.2104291,
    -5.9163382,

    ## C: limited 30-32 vs colectomy 40-79
    +0.1829689,
    -2.3039216,
    -1.2361682,
    -3.5295610,
    -7.0899623
)

if (
    max(
        abs(
            d$rd_pp -
            expected_final_rd_pp
        ),
        na.rm = TRUE
    ) >
    0.01
) {
    stop(
        "Final alternative-surgery RD fingerprint failed."
    )
}


## =============================================================================
## 8. Publication labels
## =============================================================================

definition_label <- c(
    S2_A_20_39_vs_40_79 =
        "Limited codes 20–39 vs colectomy codes 40–79",
    S2_B_exclude_code32 =
        "Exclude code 32: limited codes 30–31 vs colectomy codes 40–41",
    S2_C_colectomy_40_79 =
        "Limited codes 30–32 vs colectomy codes 40–79"
)

quintile_label <- c(
    Q1 =
        "Q1 (lowest risk)",
    Q2 =
        "Q2",
    Q3 =
        "Q3",
    Q4 =
        "Q4",
    Q5 =
        "Q5 (highest risk)"
)


## =============================================================================
## 9. Formatting helpers
## =============================================================================

fmt_n <- function(x) {
    format(
        as.integer(x),
        big.mark = ",",
        scientific = FALSE
    )
}

fmt_pp <- function(x) {
    sprintf(
        "%+.1f",
        x
    )
}

fmt_ci <- function(lo, hi) {
    sprintf(
        "%+.1f to %+.1f",
        lo,
        hi
    )
}

fmt_B <- function(x) {
    format(
        as.integer(x),
        big.mark = ",",
        scientific = FALSE
    )
}


## =============================================================================
## 10. Build final display table
## =============================================================================

table_s7 <- data.frame(
    `Alternative surgery definition` =
        unname(
            definition_label[
                d$definition
            ]
        ),
    `N (limited / colectomy)` =
        paste0(
            vapply(
                d$n,
                fmt_n,
                character(1)
            ),
            " (",
            vapply(
                d$limited_n,
                fmt_n,
                character(1)
            ),
            " / ",
            vapply(
                d$colectomy_n,
                fmt_n,
                character(1)
            ),
            ")"
        ),
    `Predicted nodal-risk quintile` =
        unname(
            quintile_label[
                d$quintile
            ]
        ),
    `RD, pp` =
        vapply(
            d$rd_pp,
            fmt_pp,
            character(1)
        ),
    `95% CI, pp` =
        mapply(
            fmt_ci,
            d$ci_lo_pp,
            d$ci_hi_pp,
            USE.NAMES = FALSE
        ),
    `Effective bootstrap B` =
        vapply(
            d$effective_B,
            fmt_B,
            character(1)
        ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 11. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S7_Alternative_Surgery_Definitions_FINAL.csv"
)

write.csv(
    table_s7,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 12. Build journal-ready three-line table
## =============================================================================

ft <- flextable(
    table_s7
)

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.0,
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.2,
    part = "header"
)

ft <- bold(
    ft,
    part = "header"
)

ft <- align(
    ft,
    j = c(
        1,
        3
    ),
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = c(
        2,
        4,
        5,
        6
    ),
    align = "center",
    part = "all"
)

ft <- valign(
    ft,
    valign = "center",
    part = "all"
)

ft <- border_remove(
    ft
)


## Three-line borders
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


## Compact padding
ft <- padding(
    ft,
    padding.top = 0.9,
    padding.bottom = 0.9,
    padding.left = 1.3,
    padding.right = 1.3,
    part = "all"
)


## Slight visual separation before definitions B and C.
ft <- padding(
    ft,
    i = c(
        6,
        11
    ),
    padding.top = 3.2,
    part = "body"
)


## Merge repeated definition and N values vertically.
ft <- merge_v(
    ft,
    j = c(
        1,
        2
    ),
    part = "body"
)

ft <- valign(
    ft,
    j = c(
        1,
        2
    ),
    valign = "center",
    part = "body"
)


## Widths
ft <- width(
    ft,
    j = 1,
    width = 3.15
)

ft <- width(
    ft,
    j = 2,
    width = 1.55
)

ft <- width(
    ft,
    j = 3,
    width = 1.55
)

ft <- width(
    ft,
    j = 4,
    width = 0.75
)

ft <- width(
    ft,
    j = 5,
    width = 1.30
)

ft <- width(
    ft,
    j = 6,
    width = 1.10
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 13. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S7_Alternative_Surgery_Definitions_FINAL.docx"
)

doc <- read_docx()


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S7. Sensitivity analyses using alternative surgery definitions",
            fp_text(
                font.family = "Times New Roman",
                font.size = 9.4,
                bold = TRUE
            )
        )
    )
)


doc <- body_add_flextable(
    doc,
    value = ft,
    align = "center"
)


note_text <- paste0(
    "RD is the absolute 5-year cancer-specific mortality risk difference, defined as ",
    "oncologic colectomy minus limited resection; negative values indicate lower cancer-specific ",
    "mortality associated with oncologic colectomy. Alternative surgery definitions were evaluated ",
    "in the corresponding pre-primary-cohort populations because some alternative surgery codes were ",
    "excluded by the primary surgery-definition step. Each sensitivity analysis used a full-pipeline ",
    "1,000-replicate bootstrap: the alternative-surgery cohort was resampled, the nodal-risk model was ",
    "refitted, predicted nodal risk was recomputed, the propensity-score and overlap-weight models were ",
    "refitted, risk quintiles were re-estimated, and weighted 5-year cancer-specific mortality differences ",
    "were recalculated. These analyses assess robustness to surgery-code definitions and do not alter the ",
    "prespecified primary surgery contrast of codes 30–32 versus 40–41."
)


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            note_text,
            fp_text(
                font.family = "Times New Roman",
                font.size = 7.4
            )
        )
    )
)


doc <- body_end_section_landscape(
    doc
)


print(
    doc,
    target = DOCX_FILE
)


## =============================================================================
## 14. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S7 COMPLETE\n")
cat("============================================================\n")

cat(
    "Source: ",
    INPUT_FILE,
    "\n",
    sep = ""
)

for (def in definition_order) {

    zz <- d[
        d$definition == def,
        ,
        drop = FALSE
    ]

    cat(
        "\n",
        definition_label[def],
        "\n",
        sep = ""
    )

    cat(
        sprintf(
            "N=%s | limited=%s | colectomy=%s\n",
            fmt_n(zz$n[1]),
            fmt_n(zz$limited_n[1]),
            fmt_n(zz$colectomy_n[1])
        )
    )

    for (i in seq_len(nrow(zz))) {

        cat(
            sprintf(
                "  %-17s | RD %+.1f pp | 95%% CI %+.1f to %+.1f | B=%d\n",
                zz$quintile[i],
                zz$rd_pp[i],
                zz$ci_lo_pp[i],
                zz$ci_hi_pp[i],
                zz$effective_B[i]
            )
        )
    }
}

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Final B=1000 results were formatted without re-estimation.\n")
cat("Standard three-line table; no vertical rules.\n")
cat("============================================================\n")
