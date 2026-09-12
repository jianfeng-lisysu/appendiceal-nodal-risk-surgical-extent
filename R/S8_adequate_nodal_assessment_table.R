## =============================================================================
## A6 - SUPPLEMENTARY TABLE S9 ONLY
## Sensitivity analysis restricting the oncologic colectomy arm to
## adequate nodal assessment (>=12 examined lymph nodes)
##
## Root:
##   <project root>
##
## INPUT:
##   03_results/06_sensitivity_core_report_FINAL.txt
##
## OUTPUT:
##   05_tables/Supplementary_Table_S9_Adequate_Nodal_Assessment_FINAL.docx
##   05_tables/Supplementary_Table_S9_Adequate_Nodal_Assessment_FINAL.csv
##
## FINAL SOURCE:
##   STEP 06 S1
##   Oncologic colectomy arm restricted to LN >=12.
##   Limited-resection arm retained without an LN>=12 restriction.
##   Input n = 6,130.
##   Effective bootstrap B = 500.
##
## FINAL RESULTS:
##   Q1 -0.6 pp (-4.0 to +2.7)
##   Q2 -5.0 pp (-9.3 to -0.2)
##   Q3 -1.5 pp (-7.9 to +3.2)
##   Q4 -6.6 pp (-10.7 to +1.2)
##   Q5 -7.0 pp (-13.7 to -0.1)
##
## IMPORTANT:
##   This script formats the already-completed sensitivity analysis.
##   It does NOT re-run bootstrap or change any estimate.
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
## 2. Locate final sensitivity report
## =============================================================================

REPORT_FILE <- file.path(
    RESULT_DIR,
    "06_sensitivity_core_report_FINAL.txt"
)

if (!file.exists(REPORT_FILE)) {
    REPORT_FILE <- file.path(
        ROOT,
        "06_sensitivity_core_report_FINAL.txt"
    )
}

if (!file.exists(REPORT_FILE)) {
    stop(
        "找不到 06_sensitivity_core_report_FINAL.txt"
    )
}


report_lines <- readLines(
    REPORT_FILE,
    warn = FALSE,
    encoding = "UTF-8"
)

report_text <- paste(
    report_lines,
    collapse = "\n"
)


## =============================================================================
## 3. Source audit
## =============================================================================

required_tokens <- c(
    "colectomy arm restricted to LN>=12",
    "input n=6130",
    "effective B=500",
    "Q1: -0.6 pp",
    "Q2: -5.0 pp",
    "Q3: -1.5 pp",
    "Q4: -6.6 pp",
    "Q5: -7.0 pp"
)

source_pass <- vapply(
    required_tokens,
    function(z) {
        grepl(
            z,
            report_text,
            fixed = TRUE
        )
    },
    logical(1)
)

if (!all(source_pass)) {

    stop(
        paste0(
            "最终 S9 来源审计未通过。未找到：",
            paste(
                required_tokens[!source_pass],
                collapse = " | "
            )
        )
    )
}


## =============================================================================
## 4. Final frozen values
## =============================================================================

s9 <- data.frame(
    quintile = paste0(
        "Q",
        1:5
    ),
    rd_pp = c(
        -0.6,
        -5.0,
        -1.5,
        -6.6,
        -7.0
    ),
    lo_pp = c(
        -4.0,
        -9.3,
        -7.9,
        -10.7,
        -13.7
    ),
    hi_pp = c(
        +2.7,
        -0.2,
        +3.2,
        +1.2,
        -0.1
    ),
    effective_B = rep(
        500L,
        5
    ),
    stringsAsFactors = FALSE
)


## =============================================================================
## 5. Hard fingerprints
## =============================================================================

if (nrow(s9) != 5) {
    stop(
        "Internal S9 row-count fingerprint failed."
    )
}

if (any(s9$effective_B != 500L)) {
    stop(
        "Internal S9 bootstrap fingerprint failed."
    )
}


## =============================================================================
## 6. Publication labels
## =============================================================================

quintile_label <- c(
    Q1 = "Q1 (lowest risk)",
    Q2 = "Q2",
    Q3 = "Q3",
    Q4 = "Q4",
    Q5 = "Q5 (highest risk)"
)


## =============================================================================
## 7. Final display table
## =============================================================================

table_s9 <- data.frame(
    `Predicted nodal-risk quintile` =
        unname(
            quintile_label[
                s9$quintile
            ]
        ),
    `RD, pp` =
        sprintf(
            "%+.1f",
            s9$rd_pp
        ),
    `95% CI, pp` =
        sprintf(
            "%+.1f to %+.1f",
            s9$lo_pp,
            s9$hi_pp
        ),
    `Effective bootstrap B` =
        format(
            s9$effective_B,
            big.mark = ",",
            scientific = FALSE
        ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 8. Save CSV
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S9_Adequate_Nodal_Assessment_FINAL.csv"
)

write.csv(
    table_s9,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 9. Journal-ready three-line table
## =============================================================================

ft <- flextable(
    table_s9
)

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.6,
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.8,
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
    j = 2:ncol(
        table_s9
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


ft <- padding(
    ft,
    padding.top = 1.4,
    padding.bottom = 1.4,
    padding.left = 1.7,
    padding.right = 1.7,
    part = "all"
)


ft <- width(
    ft,
    j = 1,
    width = 2.65
)

ft <- width(
    ft,
    j = 2,
    width = 0.90
)

ft <- width(
    ft,
    j = 3,
    width = 1.55
)

ft <- width(
    ft,
    j = 4,
    width = 1.35
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 10. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S9_Adequate_Nodal_Assessment_FINAL.docx"
)

doc <- read_docx()


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            paste0(
                "Supplementary Table S9. Sensitivity analysis restricting the ",
                "oncologic colectomy arm to adequate nodal assessment"
            ),
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
    "RD is the absolute 5-year cancer-specific mortality risk difference, defined as oncologic colectomy minus limited resection; ",
    "negative values indicate lower cancer-specific mortality associated with oncologic colectomy. ",
    "In this sensitivity analysis, only the oncologic colectomy arm was restricted to patients with at least 12 examined lymph nodes, ",
    "whereas the limited-resection arm was retained without an analogous lymph-node-count restriction. ",
    "The resulting analysis included 6,130 patients. The nodal-risk model, propensity-score model, overlap weights, and predicted-risk ",
    "quintile cut points were re-estimated within each bootstrap replicate. ",
    "Confidence intervals were obtained from 500 full-pipeline bootstrap replicates. ",
    "This analysis was designed to assess sensitivity to the adequacy of nodal assessment in the oncologic colectomy arm and should not ",
    "be described as a cohort in which all patients had at least 12 examined lymph nodes."
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


doc <- body_end_section_landscape(
    doc
)


print(
    doc,
    target = DOCX_FILE
)


## =============================================================================
## 11. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S9 COMPLETE\n")
cat("============================================================\n")

cat(
    "Source report: ",
    REPORT_FILE,
    "\n",
    sep = ""
)

cat(
    "Sensitivity input N = 6,130\n"
)

cat(
    "Definition: colectomy arm restricted to LN>=12; limited arm retained.\n\n"
)

for (i in seq_len(nrow(s9))) {

    cat(
        sprintf(
            "%-17s | RD %+.1f pp | 95%% CI %+.1f to %+.1f | B=%d\n",
            quintile_label[
                s9$quintile[i]
            ],
            s9$rd_pp[i],
            s9$lo_pp[i],
            s9$hi_pp[i],
            s9$effective_B[i]
        )
    )
}

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Final B=500 sensitivity estimates were formatted without re-estimation.\n")
cat("Standard three-line table; no vertical rules.\n")
cat("============================================================\n")
