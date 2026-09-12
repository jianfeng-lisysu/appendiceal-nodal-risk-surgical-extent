## =============================================================================
## A6 - SUPPLEMENTARY TABLE S4 ONLY
## Landmark sensitivity analyses for prespecified low-risk groups
## JOURNAL-READY THREE-LINE TABLE
##
## Root:
##   <project root>
##
## SOURCE
##   03_results/05_landmark_report_FINAL_B1000.txt
##
## OUTPUTS
##   05_tables/Supplementary_Table_S4_Landmark_Sensitivity_FINAL.docx
##   05_tables/Supplementary_Table_S4_Landmark_Sensitivity_FINAL.csv
##
## FINAL STEP-05 RESULTS
##   Landmark 3 months:
##     Q1                  -0.1 pp (-3.6 to +2.9)
##     Bottom decile      +2.2 pp (-2.3 to +5.8)
##     Risk <5%           +3.2 pp (-2.3 to +6.7)
##
##   Landmark 6 months:
##     Q1                  -0.5 pp (-3.7 to +2.4)
##     Bottom decile      +2.1 pp (-2.3 to +5.5)
##     Risk <5%           +3.1 pp (-2.5 to +6.6)
##
##   Landmark 12 months:
##     Q1                  -0.4 pp (-3.7 to +2.2)
##     Bottom decile      +1.9 pp (-2.7 to +5.3)
##     Risk <5%           +2.8 pp (-2.5 to +6.5)
##
##   Effective bootstrap B = 1000 for every estimate.
##
## IMPORTANT
##   - Outcome horizon remains 60 months after diagnosis.
##   - L=0 is the primary analysis and is not repeated here.
##   - Percentile-based risk cut points were re-estimated within every
##     full-pipeline bootstrap replicate.
##   - RD = oncologic colectomy - limited resection.
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

RESULT_DIR <- file.path(
    ROOT,
    "03_results"
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

library(officer)
library(flextable)


## =============================================================================
## 2. Locate the final Step-05 report
## =============================================================================

REPORT_FILE <- file.path(
    RESULT_DIR,
    "05_landmark_report_FINAL_B1000.txt"
)

if (
    !file.exists(
        REPORT_FILE
    )
) {

    REPORT_FILE <- file.path(
        ROOT,
        "05_landmark_report_FINAL_B1000.txt"
    )
}

if (
    !file.exists(
        REPORT_FILE
    )
) {

    stop(
        paste0(
            "找不到最终 landmark 报告：05_landmark_report_FINAL_B1000.txt\n",
            "请确认该文件位于 03_results/ 或 A6 根目录。"
        )
    )
}


report_text <- paste(
    readLines(
        REPORT_FILE,
        warn = FALSE,
        encoding = "UTF-8"
    ),
    collapse = "\n"
)


## =============================================================================
## 3. Basic source-report audit
##
## This table uses the already-frozen Step-05 estimates. We require the final
## report to contain the expected landmark labels before writing the table.
## =============================================================================

required_landmark_tokens <- c(
    "3",
    "6",
    "12",
    "1000"
)

if (
    !all(
        vapply(
            required_landmark_tokens,
            function(z) grepl(
                z,
                report_text,
                fixed = TRUE
            ),
            logical(1)
        )
    )
) {

    stop(
        "Final landmark report did not pass the basic landmark/B=1000 audit."
    )
}


## =============================================================================
## 4. Encode the final frozen Step-05 estimates
##
## These values are the final B=1000 outputs from
## 05_landmark_report_FINAL_B1000.txt.
## =============================================================================

s4 <- data.frame(
    landmark_month = c(
        3, 3, 3,
        6, 6, 6,
        12, 12, 12
    ),
    analysis_id = c(
        "Q1",
        "bottom_decile",
        "risk_lt_5pct",
        "Q1",
        "bottom_decile",
        "risk_lt_5pct",
        "Q1",
        "bottom_decile",
        "risk_lt_5pct"
    ),
    rd_pp = c(
        -0.1,
        +2.2,
        +3.2,
        -0.5,
        +2.1,
        +3.1,
        -0.4,
        +1.9,
        +2.8
    ),
    lo_pp = c(
        -3.6,
        -2.3,
        -2.3,
        -3.7,
        -2.3,
        -2.5,
        -3.7,
        -2.7,
        -2.5
    ),
    hi_pp = c(
        +2.9,
        +5.8,
        +6.7,
        +2.4,
        +5.5,
        +6.6,
        +2.2,
        +5.3,
        +6.5
    ),
    effective_B = rep(
        1000L,
        9
    ),
    stringsAsFactors = FALSE
)


## =============================================================================
## 5. Hard fingerprints
## =============================================================================

expected_rd <- c(
    -0.1, 2.2, 3.2,
    -0.5, 2.1, 3.1,
    -0.4, 1.9, 2.8
)

expected_lo <- c(
    -3.6, -2.3, -2.3,
    -3.7, -2.3, -2.5,
    -3.7, -2.7, -2.5
)

expected_hi <- c(
    2.9, 5.8, 6.7,
    2.4, 5.5, 6.6,
    2.2, 5.3, 6.5
)

if (
    max(
        abs(
            s4$rd_pp -
            expected_rd
        )
    ) >
    1e-12 ||
    max(
        abs(
            s4$lo_pp -
            expected_lo
        )
    ) >
    1e-12 ||
    max(
        abs(
            s4$hi_pp -
            expected_hi
        )
    ) >
    1e-12
) {

    stop(
        "Internal Step-05 landmark fingerprint failed."
    )
}


## =============================================================================
## 6. Derived 3-pp criterion
##
## RD = colectomy - limited.
## A negative lower confidence limit represents the largest plausible
## colectomy-associated reduction. Its magnitude is -lower CI.
## Criterion met if this upper bound is <3 pp.
## =============================================================================

s4$upper_reduction_pp <- pmax(
    0,
    -s4$lo_pp
)

s4$criterion <- ifelse(
    s4$upper_reduction_pp <
    3,
    "Met",
    "Not met"
)


## Expected criterion pattern:
## Q1: not met at 3/6/12
## bottom decile: met at 3/6/12
## risk <5%: met at 3/6/12
expected_criterion <- c(
    "Not met",
    "Met",
    "Met",
    "Not met",
    "Met",
    "Met",
    "Not met",
    "Met",
    "Met"
)

if (
    !identical(
        s4$criterion,
        expected_criterion
    )
) {

    stop(
        "3-pp criterion fingerprint failed."
    )
}


## =============================================================================
## 7. Publication labels
## =============================================================================

analysis_label <- c(
    Q1 =
        "Q1 (lowest quintile)",
    bottom_decile =
        "Bottom predicted-risk decile",
    risk_lt_5pct =
        "Predicted nodal risk <5%"
)


## =============================================================================
## 8. Final display table
## =============================================================================

table_s4 <- data.frame(
    `Landmark, months` =
        s4$landmark_month,
    `Low-risk analysis` =
        unname(
            analysis_label[
                s4$analysis_id
            ]
        ),
    `RD, pp` =
        sprintf(
            "%+.1f",
            s4$rd_pp
        ),
    `95% CI, pp` =
        sprintf(
            "%+.1f to %+.1f",
            s4$lo_pp,
            s4$hi_pp
        ),
    `Upper bound of colectomy-associated reduction, pp` =
        sprintf(
            "%.1f",
            s4$upper_reduction_pp
        ),
    `3-pp criterion` =
        s4$criterion,
    `Effective bootstrap B` =
        s4$effective_B,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 9. Save CSV
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S4_Landmark_Sensitivity_FINAL.csv"
)

write.csv(
    table_s4,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 10. Build journal-ready three-line table
## =============================================================================

ft <- flextable(
    table_s4
)

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.2,
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.4,
    part = "header"
)

ft <- bold(
    ft,
    part = "header"
)

ft <- align(
    ft,
    j = 1,
    align = "center",
    part = "all"
)

ft <- align(
    ft,
    j = 2,
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 3:ncol(
        table_s4
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


## Three-line rules
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


## Add a little visual separation before the 6- and 12-month blocks
ft <- padding(
    ft,
    padding.top = 1.2,
    padding.bottom = 1.2,
    padding.left = 1.5,
    padding.right = 1.5,
    part = "all"
)

ft <- padding(
    ft,
    i = c(
        4,
        7
    ),
    padding.top = 4.0,
    part = "body"
)


## Column widths
ft <- width(
    ft,
    j = 1,
    width = 0.85
)

ft <- width(
    ft,
    j = 2,
    width = 2.25
)

ft <- width(
    ft,
    j = 3,
    width = 0.70
)

ft <- width(
    ft,
    j = 4,
    width = 1.30
)

ft <- width(
    ft,
    j = 5,
    width = 2.15
)

ft <- width(
    ft,
    j = 6,
    width = 0.95
)

ft <- width(
    ft,
    j = 7,
    width = 1.05
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 11. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S4_Landmark_Sensitivity_FINAL.docx"
)

doc <- read_docx()


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S4. Landmark sensitivity analyses for prespecified low-risk groups",
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
    "RD is the absolute cancer-specific mortality risk difference, defined as oncologic colectomy minus limited resection; ",
    "negative values indicate lower cancer-specific mortality associated with oncologic colectomy. ",
    "Landmark analyses conditioned on remaining at risk at 3, 6, or 12 months after diagnosis, while the outcome horizon ",
    "remained 60 months after diagnosis. The primary L=0 analysis is reported separately and was not repeated here. ",
    "For percentile-based analyses (Q1 and bottom predicted-risk decile), risk cut points were re-estimated within each ",
    "full-pipeline bootstrap replicate; the <5% threshold remained fixed. ",
    "The upper bound of the colectomy-associated absolute reduction is the magnitude of the lower 95% confidence limit. ",
    "The prespecified 3-percentage-point criterion was met when this upper bound was <3 percentage points."
)


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            note_text,
            fp_text(
                font.family = "Times New Roman",
                font.size = 7.5
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
## 12. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S4 COMPLETE\n")
cat("============================================================\n")

cat(
    "Source report: ",
    REPORT_FILE,
    "\n",
    sep = ""
)

cat("\nFinal landmark estimates:\n")

for (
    i in seq_len(
        nrow(
            s4
        )
    )
) {

    cat(
        sprintf(
            "L=%2d mo | %-30s | RD %+.1f pp (%+.1f to %+.1f) | upper reduction %.1f pp | %s | B=%d\n",
            s4$landmark_month[i],
            analysis_label[
                s4$analysis_id[i]
            ],
            s4$rd_pp[i],
            s4$lo_pp[i],
            s4$hi_pp[i],
            s4$upper_reduction_pp[i],
            s4$criterion[i],
            s4$effective_B[i]
        )
    )
}

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Standard three-line table; no vertical rules.\n")
cat("Final Step-05 B=1000 landmark values used without re-estimation.\n")
cat("============================================================\n")
