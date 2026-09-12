## =============================================================================
## A6 - TABLE 2 ONLY
## Primary risk-stratified 5-year cancer-specific mortality differences
## JOURNAL-READY THREE-LINE TABLE — FINAL LOCK VERSION
##
## Root:
##   <project root>
##
## INPUTS
##   03_results/04_curve_data_FINAL_B1000.csv
##   02_data/analytic_cohort_A_weighted.csv
##
## OUTPUT
##   05_tables/Table2_Primary_5yr_CSM_RD_FINAL.docx
##   05_tables/Table2_Primary_5yr_CSM_RD_FINAL.csv
##
## STANDARD
##   - Standard three-line table
##   - No vertical rules
##   - Times New Roman
##   - Single title only
##   - One-page portrait target
##   - Frozen B=1000 results only
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

needed_packages <- c(
    "officer",
    "flextable"
)

missing_packages <- needed_packages[
    !vapply(
        needed_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (
    length(
        missing_packages
    ) >
    0
) {

    stop(
        "请先安装：",
        paste(
            missing_packages,
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
## 2. Locate files
## =============================================================================

MAIN_FILE <- file.path(
    RESULT_DIR,
    "04_curve_data_FINAL_B1000.csv"
)

if (
    !file.exists(
        MAIN_FILE
    )
) {

    MAIN_FILE <- file.path(
        ROOT,
        "04_curve_data_FINAL_B1000.csv"
    )
}

if (
    !file.exists(
        MAIN_FILE
    )
) {

    stop(
        "找不到 04_curve_data_FINAL_B1000.csv"
    )
}


WT_FILE <- file.path(
    DATA_DIR,
    "analytic_cohort_A_weighted.csv"
)

if (
    !file.exists(
        WT_FILE
    )
) {

    WT_FILE <- file.path(
        ROOT,
        "analytic_cohort_A_weighted.csv"
    )
}

if (
    !file.exists(
        WT_FILE
    )
) {

    stop(
        "找不到 analytic_cohort_A_weighted.csv"
    )
}


## =============================================================================
## 3. Read frozen data
## =============================================================================

MAIN <- read.csv(
    MAIN_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

WT <- read.csv(
    WT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


needed_main <- c(
    "stratum",
    "rd",
    "lo",
    "hi"
)

if (
    !all(
        needed_main %in%
        names(
            MAIN
        )
    )
) {

    stop(
        "04_curve_data_FINAL_B1000.csv 缺少字段：",
        paste(
            setdiff(
                needed_main,
                names(
                    MAIN
                )
            ),
            collapse = ", "
        )
    )
}

if (
    nrow(
        WT
    ) !=
    6951
) {

    stop(
        "Frozen cohort fingerprint failed: expected N=6951."
    )
}


## =============================================================================
## 4. Freeze order + hard fingerprint
## =============================================================================

expected_order <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

if (
    !all(
        expected_order %in%
        MAIN$stratum
    )
) {

    stop(
        "7个主分析strata不完整。"
    )
}

MAIN <- MAIN[
    match(
        expected_order,
        MAIN$stratum
    ),
    ,
    drop = FALSE
]


## Final frozen point estimates and CIs from STEP 04 B=1000
EXPECTED_RD_PP <- c(
    0.3,
    -3.6,
    -0.7,
    -5.1,
    -6.0,
    3.1,
    4.2
)

EXPECTED_LO_PP <- c(
    -3.4,
    -7.8,
    -7.4,
    -9.6,
    -12.7,
    -2.2,
    -2.0
)

EXPECTED_HI_PP <- c(
    3.0,
    1.6,
    3.8,
    2.1,
    0.5,
    7.0,
    7.8
)


OBS_RD_PP <- MAIN$rd *
    100

OBS_LO_PP <- MAIN$lo *
    100

OBS_HI_PP <- MAIN$hi *
    100


if (
    max(
        abs(
            OBS_RD_PP -
            EXPECTED_RD_PP
        )
    ) >
    0.15
) {

    stop(
        "RD point-estimate fingerprint failed."
    )
}


if (
    max(
        abs(
            OBS_LO_PP -
            EXPECTED_LO_PP
        )
    ) >
    0.15
) {

    stop(
        "Lower CI fingerprint failed."
    )
}


if (
    max(
        abs(
            OBS_HI_PP -
            EXPECTED_HI_PP
        )
    ) >
    0.15
) {

    stop(
        "Upper CI fingerprint failed."
    )
}


## =============================================================================
## 5. Original-sample risk boundaries for display only
##
## IMPORTANT:
## Percentile cut points were re-estimated inside every bootstrap replicate.
## These boundaries are only descriptive values from the frozen original sample.
## =============================================================================

risk_q <- as.numeric(
    quantile(
        WT$nodal_risk,
        probs = c(
            0.10,
            0.20,
            0.40,
            0.60,
            0.80
        ),
        na.rm = TRUE,
        names = FALSE
    )
) *
    100


P10 <- risk_q[
    1
]

P20 <- risk_q[
    2
]

P40 <- risk_q[
    3
]

P60 <- risk_q[
    4
]

P80 <- risk_q[
    5
]


## =============================================================================
## 6. Build display table
## =============================================================================

risk_label <- c(
    Q1 =
        "Q1 (lowest risk)",
    Q2 =
        "Q2",
    Q3 =
        "Q3",
    Q4 =
        "Q4",
    Q5 =
        "Q5 (highest risk)",
    bottom_decile =
        "Bottom predicted-risk decile",
    risk_lt_5pct =
        "Predicted nodal risk <5%"
)


risk_definition <- c(
    Q1 =
        sprintf(
            "≤%.1f%%",
            P20
        ),
    Q2 =
        sprintf(
            ">%.1f%% to ≤%.1f%%",
            P20,
            P40
        ),
    Q3 =
        sprintf(
            ">%.1f%% to ≤%.1f%%",
            P40,
            P60
        ),
    Q4 =
        sprintf(
            ">%.1f%% to ≤%.1f%%",
            P60,
            P80
        ),
    Q5 =
        sprintf(
            ">%.1f%%",
            P80
        ),
    bottom_decile =
        sprintf(
            "≤P10 (%.1f%%)",
            P10
        ),
    risk_lt_5pct =
        "<5.0%"
)


## Colectomy-associated advantage upper bound
## RD = colectomy - limited
## Negative RD favors colectomy
## Thus the largest plausible absolute colectomy advantage = - lower CI bound
upper_advantage_pp <- (
    -MAIN$lo *
    100
)


criterion <- rep(
    "—",
    nrow(
        MAIN
    )
)

criterion_rows <- MAIN$stratum %in%
    c(
        "Q1",
        "bottom_decile",
        "risk_lt_5pct"
    )

criterion[
    criterion_rows
] <- ifelse(
    upper_advantage_pp[
        criterion_rows
    ] <
    3,
    "Met",
    "Not met"
)


table2 <- data.frame(
    `Risk stratum` =
        unname(
            risk_label[
                MAIN$stratum
            ]
        ),
    `Predicted nodal risk` =
        unname(
            risk_definition[
                MAIN$stratum
            ]
        ),
    `RD, pp` =
        sprintf(
            "%+.1f",
            MAIN$rd *
            100
        ),
    `95% CI, pp` =
        sprintf(
            "%+.1f to %+.1f",
            MAIN$lo *
            100,
            MAIN$hi *
            100
        ),
    `Upper bound of colectomy-associated reduction, pp` =
        sprintf(
            "%.1f",
            upper_advantage_pp
        ),
    `3-pp criterion` =
        criterion,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 7. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Table2_Primary_5yr_CSM_RD_FINAL.csv"
)

write.csv(
    table2,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 8. Build journal-style three-line table
## =============================================================================

ft <- flextable(
    table2
)


## -----------------------------------------------------------------------------
## Font
## -----------------------------------------------------------------------------

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.8,
    part = "all"
)

ft <- fontsize(
    ft,
    size = 9.0,
    part = "header"
)

ft <- bold(
    ft,
    part = "header"
)


## -----------------------------------------------------------------------------
## Alignment
## -----------------------------------------------------------------------------

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
## Three-line borders only
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
## Slight separation before prespecified extreme-low-risk anchors
##
## We keep a pure three-line table: no internal horizontal rule.
## Instead use a small top padding increase for rows 6 and 7.
## -----------------------------------------------------------------------------

ft <- padding(
    ft,
    padding.top = 2.1,
    padding.bottom = 2.1,
    padding.left = 2.0,
    padding.right = 2.0,
    part = "all"
)

ft <- padding(
    ft,
    i = 6,
    padding.top = 5.5,
    part = "body"
)


## -----------------------------------------------------------------------------
## Column widths
## -----------------------------------------------------------------------------

ft <- width(
    ft,
    j = 1,
    width = 1.75
)

ft <- width(
    ft,
    j = 2,
    width = 1.55
)

ft <- width(
    ft,
    j = 3,
    width = 0.70
)

ft <- width(
    ft,
    j = 4,
    width = 1.25
)

ft <- width(
    ft,
    j = 5,
    width = 1.75
)

ft <- width(
    ft,
    j = 6,
    width = 1.00
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 9. Build Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Table2_Primary_5yr_CSM_RD_FINAL.docx"
)

doc <- read_docx()


## -----------------------------------------------------------------------------
## Single title only
## -----------------------------------------------------------------------------

title_text <- fpar(
    ftext(
        "Table 2. Primary risk-stratified 5-year cancer-specific mortality differences",
        fp_text(
            font.family = "Times New Roman",
            font.size = 10,
            bold = TRUE
        )
    )
)

doc <- body_add_fpar(
    doc,
    title_text
)


## -----------------------------------------------------------------------------
## Table
## -----------------------------------------------------------------------------

doc <- body_add_flextable(
    doc,
    value = ft,
    align = "center"
)


## -----------------------------------------------------------------------------
## Footnote
## -----------------------------------------------------------------------------

note_text <- paste0(
    "RD is the absolute 5-year cancer-specific mortality risk difference, defined as ",
    "oncologic colectomy minus limited resection; negative values indicate lower cancer-specific mortality ",
    "associated with oncologic colectomy. The upper bound of the colectomy-associated absolute reduction is the magnitude of the lower ",
    "95% confidence limit. The prespecified 3-percentage-point criterion was evaluated for Q1, the bottom ",
    "predicted-risk decile, and predicted nodal risk <5%. Displayed percentile boundaries are descriptive cut ",
    "points from the original analytic sample; percentile cut points were re-estimated within every full-pipeline ",
    "bootstrap replicate."
)

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            note_text,
            fp_text(
                font.family = "Times New Roman",
                font.size = 8.2
            )
        )
    )
)


## -----------------------------------------------------------------------------
## Portrait section with compact margins
## -----------------------------------------------------------------------------

section_portrait <- prop_section(
    page_size = page_size(
        orient = "portrait"
    ),
    page_margins = page_mar(
        top = 0.55,
        bottom = 0.55,
        left = 0.55,
        right = 0.55
    ),
    type = "continuous"
)

doc <- body_end_block_section(
    doc,
    value = block_section(
        section_portrait
    )
)


print(
    doc,
    target = DOCX_FILE
)



## =============================================================================
## 10A. Wording audit
## =============================================================================

wording_audit <- paste0(
    "Header: Upper bound of colectomy-associated reduction, pp | ",
    "Footnote: The upper bound of the colectomy-associated absolute reduction ",
    "is the magnitude of the lower 95% confidence limit."
)

cat("\nWORDING AUDIT\n")
cat(wording_audit, "\n")

## =============================================================================
## 10. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("TABLE 2 FINAL LOCK VERSION COMPLETE\n")
cat("============================================================\n")

for (
    i in seq_len(
        nrow(
            table2
        )
    )
) {

    cat(
        sprintf(
            "%-30s | %-18s | RD %s pp | 95%% CI %s | upper %s | %s\n",
            table2[
                i,
                "Risk stratum"
            ],
            table2[
                i,
                "Predicted nodal risk"
            ],
            table2[
                i,
                "RD, pp"
            ],
            table2[
                i,
                "95% CI, pp"
            ],
            table2[
                i,
                "Upper bound of colectomy-associated reduction, pp"
            ],
            table2[
                i,
                "3-pp criterion"
            ]
        )
    )
}

cat("\nOriginal-sample descriptive risk cut points:\n")

cat(
    sprintf(
        "P10=%.3f%% | P20=%.3f%% | P40=%.3f%% | P60=%.3f%% | P80=%.3f%%\n",
        P10,
        P20,
        P40,
        P60,
        P80
    )
)

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Standard: journal-ready three-line table, no vertical rules.\n")
cat("Frozen B=1000 RD/CI fingerprint: PASS.\n")
cat("============================================================\n")
