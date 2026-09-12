## =============================================================================
## A6 - TABLE 1 ONLY
## Baseline characteristics of the final analytic cohort
## JOURNAL-READY THREE-LINE TABLE — FINAL LOCK VERSION
##
## Root:
##   <project root>
##
## INPUT
##   02_data/analytic_cohort_A_weighted.csv
##
## OUTPUT
##   05_tables/Table1_Baseline_Characteristics_FINAL.docx
##   05_tables/Table1_Baseline_Characteristics_FINAL.csv
##
## STANDARD
##   - Standard three-line table
##   - No vertical rules
##   - Times New Roman
##   - Landscape page
##   - Compact one-page target
##   - No duplicate title
##   - Clean journal-style headings
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
TABLE_DIR <- file.path(ROOT, "05_tables")

dir.create(
    TABLE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

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

if (length(missing_packages) > 0) {
    stop(
        "请先安装：",
        paste(
            missing_packages,
            collapse = ", "
        )
    )
}

library(officer)
library(flextable)


## =============================================================================
## 1. Read frozen weighted cohort
## =============================================================================

WT_FILE <- file.path(
    DATA_DIR,
    "analytic_cohort_A_weighted.csv"
)

if (!file.exists(WT_FILE)) {
    WT_FILE <- file.path(
        ROOT,
        "analytic_cohort_A_weighted.csv"
    )
}

if (!file.exists(WT_FILE)) {
    stop(
        "找不到 analytic_cohort_A_weighted.csv"
    )
}

WT <- read.csv(
    WT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(WT) != 6951) {
    stop(
        "Frozen cohort fingerprint failed: expected N=6951."
    )
}

if (
    sum(WT$exposure == "limited") != 2877 ||
    sum(WT$exposure == "colectomy") != 4074
) {
    stop(
        "Frozen treatment-arm fingerprint failed."
    )
}

WT$treat <- as.integer(
    WT$exposure ==
    "colectomy"
)

limited_mask <- (
    WT$treat ==
    0
)

colectomy_mask <- (
    WT$treat ==
    1
)


## =============================================================================
## 2. Statistical helpers
## =============================================================================

weighted_mean_local <- function(
    x,
    w
) {

    ok <- (
        !is.na(x) &
        !is.na(w) &
        is.finite(w)
    )

    if (sum(ok) == 0) {
        return(
            NA_real_
        )
    }

    sum(
        x[ok] *
        w[ok]
    ) /
    sum(
        w[ok]
    )
}


weighted_var_local <- function(
    x,
    w
) {

    ok <- (
        !is.na(x) &
        !is.na(w) &
        is.finite(w)
    )

    x <- x[ok]
    w <- w[ok]

    if (
        length(x) <
        2 ||
        sum(w) <=
        0
    ) {
        return(
            NA_real_
        )
    }

    m <- weighted_mean_local(
        x,
        w
    )

    sum(
        w *
        (
            x -
            m
        )^2
    ) /
    sum(w)
}


smd_local <- function(
    x,
    treat,
    w
) {

    ok <- (
        !is.na(x) &
        !is.na(treat) &
        !is.na(w) &
        is.finite(w)
    )

    x <- x[ok]
    treat <- treat[ok]
    w <- w[ok]

    if (
        sum(treat == 1) ==
        0 ||
        sum(treat == 0) ==
        0
    ) {
        return(
            NA_real_
        )
    }

    m1 <- weighted_mean_local(
        x[
            treat ==
            1
        ],
        w[
            treat ==
            1
        ]
    )

    m0 <- weighted_mean_local(
        x[
            treat ==
            0
        ],
        w[
            treat ==
            0
        ]
    )

    v1 <- weighted_var_local(
        x[
            treat ==
            1
        ],
        w[
            treat ==
            1
        ]
    )

    v0 <- weighted_var_local(
        x[
            treat ==
            0
        ],
        w[
            treat ==
            0
        ]
    )

    den <- sqrt(
        (
            v1 +
            v0
        ) /
        2
    )

    if (
        !is.finite(den) ||
        den ==
        0
    ) {
        return(
            0
        )
    }

    (
        m1 -
        m0
    ) /
    den
}


fmt_cont <- function(
    x
) {

    x <- x[
        is.finite(x)
    ]

    if (length(x) == 0) {
        return("")
    }

    sprintf(
        "%.1f (%.1f)",
        mean(x),
        sd(x)
    )
}


fmt_pct <- function(
    x
) {

    x <- as.numeric(x)

    denom <- sum(
        !is.na(x)
    )

    num <- sum(
        x ==
        1,
        na.rm = TRUE
    )

    if (denom == 0) {
        return("")
    }

    sprintf(
        "%d (%.1f%%)",
        num,
        100 *
        num /
        denom
    )
}


fmt_smd <- function(
    x
) {

    z <- round(
        x,
        3
    )

    if (
        is.finite(z) &&
        abs(z) <
        0.0005
    ) {
        z <- 0
    }

    if (!is.finite(z)) {
        return("")
    }

    sprintf(
        "%.3f",
        z
    )
}


## =============================================================================
## 3. Build Table 1 rows
## =============================================================================

rows <- data.frame(
    Characteristic = character(0),
    Limited = character(0),
    Colectomy = character(0),
    Unweighted_SMD = character(0),
    Weighted_SMD = character(0),
    row_type = character(0),
    stringsAsFactors = FALSE
)


add_group <- function(
    label
) {

    rows <<- rbind(
        rows,
        data.frame(
            Characteristic = label,
            Limited = "",
            Colectomy = "",
            Unweighted_SMD = "",
            Weighted_SMD = "",
            row_type = "group",
            stringsAsFactors = FALSE
        )
    )
}


add_continuous <- function(
    label,
    x
) {

    rows <<- rbind(
        rows,
        data.frame(
            Characteristic = label,
            Limited =
                fmt_cont(
                    x[
                        limited_mask
                    ]
                ),
            Colectomy =
                fmt_cont(
                    x[
                        colectomy_mask
                    ]
                ),
            Unweighted_SMD =
                fmt_smd(
                    smd_local(
                        x,
                        WT$treat,
                        rep(
                            1,
                            nrow(WT)
                        )
                    )
                ),
            Weighted_SMD =
                fmt_smd(
                    smd_local(
                        x,
                        WT$treat,
                        WT$ow
                    )
                ),
            row_type = "data",
            stringsAsFactors = FALSE
        )
    )
}


add_binary <- function(
    label,
    indicator
) {

    indicator <- as.numeric(
        indicator
    )

    rows <<- rbind(
        rows,
        data.frame(
            Characteristic = label,
            Limited =
                fmt_pct(
                    indicator[
                        limited_mask
                    ]
                ),
            Colectomy =
                fmt_pct(
                    indicator[
                        colectomy_mask
                    ]
                ),
            Unweighted_SMD =
                fmt_smd(
                    smd_local(
                        indicator,
                        WT$treat,
                        rep(
                            1,
                            nrow(WT)
                        )
                    )
                ),
            Weighted_SMD =
                fmt_smd(
                    smd_local(
                        indicator,
                        WT$treat,
                        WT$ow
                    )
                ),
            row_type = "data",
            stringsAsFactors = FALSE
        )
    )
}


## -----------------------------------------------------------------------------
## Core characteristics
## -----------------------------------------------------------------------------

add_continuous(
    "Age, years",
    WT$age
)

add_binary(
    "Female sex",
    WT$female
)

add_continuous(
    "Year of diagnosis",
    WT$year
)

add_continuous(
    "Tumor size, mm (imputed)",
    WT$size_i
)

add_binary(
    "Tumor size missing",
    WT$size_miss
)

add_continuous(
    "Predicted nodal risk, %",
    WT$nodal_risk *
    100
)


## -----------------------------------------------------------------------------
## T stage
## -----------------------------------------------------------------------------

add_group(
    "T stage"
)

for (
    lev in c(
        "T1",
        "T2",
        "T3",
        "T4"
    )
) {

    add_binary(
        paste0(
            "    ",
            lev
        ),
        WT$T ==
        lev
    )
}


## -----------------------------------------------------------------------------
## Grade
## -----------------------------------------------------------------------------

add_group(
    "Tumor grade"
)

grade_chr <- as.character(
    WT$grade_f
)

grade_labels <- c(
    "0" = "    Unknown",
    "1" = "    Grade 1",
    "2" = "    Grade 2",
    "3" = "    Grade 3",
    "4" = "    Grade 4"
)

for (
    lev in names(
        grade_labels
    )
) {

    add_binary(
        grade_labels[
            lev
        ],
        grade_chr ==
        lev
    )
}


## -----------------------------------------------------------------------------
## Histology
## -----------------------------------------------------------------------------

add_group(
    "Histology"
)

hist_labels <- c(
    nonmucinous =
        "    Nonmucinous adenocarcinoma",
    mucinous =
        "    Mucinous adenocarcinoma",
    GCA =
        "    Goblet cell adenocarcinoma",
    SRCC =
        "    Signet-ring cell carcinoma"
)

for (
    lev in names(
        hist_labels
    )
) {

    add_binary(
        hist_labels[
            lev
        ],
        WT$hist_group ==
        lev
    )
}


## -----------------------------------------------------------------------------
## Stage
## -----------------------------------------------------------------------------

add_binary(
    "Regional stage",
    WT$stage_sum ==
    "Regional"
)


## -----------------------------------------------------------------------------
## Race / ethnicity
## -----------------------------------------------------------------------------

add_group(
    "Race/ethnicity"
)

race_values <- as.character(
    WT$race4
)

race_levels <- c(
    "Hispanic (All Races)",
    "Non-Hispanic Black",
    "Non-Hispanic White",
    "Other/Unknown"
)

race_labels <- c(
    "Hispanic (All Races)" =
        "    Hispanic",
    "Non-Hispanic Black" =
        "    Non-Hispanic Black",
    "Non-Hispanic White" =
        "    Non-Hispanic White",
    "Other/Unknown" =
        "    Other/unknown"
)

for (
    lev in race_levels
) {

    if (
        lev %in%
        race_values
    ) {

        add_binary(
            race_labels[
                lev
            ],
            race_values ==
            lev
        )
    }
}


## -----------------------------------------------------------------------------
## Socioeconomic / marital variables
## -----------------------------------------------------------------------------

add_binary(
    "Higher-income category",
    WT$income_hi
)

add_binary(
    "Metropolitan residence",
    WT$metro
)

add_binary(
    "Married",
    WT$married
)


## =============================================================================
## 4. Final clean table object
## =============================================================================

table1 <- rows[
    ,
    c(
        "Characteristic",
        "Limited",
        "Colectomy",
        "Unweighted_SMD",
        "Weighted_SMD"
    )
]

names(
    table1
) <- c(
    "Characteristic",
    "Limited resection\n(n = 2,877)",
    "Oncologic colectomy\n(n = 4,074)",
    "Unweighted\nSMD",
    "Overlap-weighted\nSMD"
)


## =============================================================================
## 5. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Table1_Baseline_Characteristics_FINAL.csv"
)

write.csv(
    table1,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 6. Build journal-style flextable
## =============================================================================

ft <- flextable(
    table1
)

## Font
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
    size = 8.3,
    part = "header"
)

## Alignment
ft <- align(
    ft,
    j = 1,
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 2:5,
    align = "center",
    part = "all"
)

ft <- valign(
    ft,
    valign = "center",
    part = "all"
)

ft <- bold(
    ft,
    part = "header"
)

## Remove every border first
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

## Three-line table
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

## Group headings
group_rows <- which(
    rows$row_type ==
    "group"
)

ft <- bold(
    ft,
    i = group_rows,
    j = 1,
    bold = TRUE,
    part = "body"
)

## No horizontal rule inside the body except the three principal rules
## Group rows are distinguished only by bold text and spacing.

## Padding / compactness
ft <- padding(
    ft,
    padding.top = 1.1,
    padding.bottom = 1.1,
    padding.left = 1.8,
    padding.right = 1.8,
    part = "all"
)

## Column widths
ft <- width(
    ft,
    j = 1,
    width = 2.75
)

ft <- width(
    ft,
    j = 2,
    width = 1.55
)

ft <- width(
    ft,
    j = 3,
    width = 1.70
)

ft <- width(
    ft,
    j = 4,
    width = 1.05
)

ft <- width(
    ft,
    j = 5,
    width = 1.20
)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)



## =============================================================================
## 7. Build Word document — single title only
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Table1_Baseline_Characteristics_FINAL.docx"
)

doc <- read_docx()


## Title
title_text <- fpar(
    ftext(
        "Table 1. Baseline characteristics of the final analytic cohort",
        fp_text(
            font.family = "Times New Roman",
            font.size = 9.5,
            bold = TRUE
        )
    )
)

doc <- body_add_fpar(
    doc,
    title_text
)

## Table
doc <- body_add_flextable(
    doc,
    value = ft,
    align = "center"
)

## Footnote
note_text <- paste0(
    "Values are mean (SD) or n (%). ",
    "GCA indicates goblet cell adenocarcinoma; SMD, standardized mean difference; ",
    "SRCC, signet-ring cell carcinoma. ",
    "Overlap-weighted SMDs are shown after weighting using the frozen primary propensity-score model."
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

## Apply a landscape section using officer's broadly compatible helper.
## The section break is appended at the end and applies to the preceding content.
doc <- body_end_section_landscape(
    doc
)

print(
    doc,
    target = DOCX_FILE
)


## =============================================================================
## 8. Audit
## =============================================================================

weighted_smd_values <- suppressWarnings(
    as.numeric(
        table1[
            ,
            "Overlap-weighted\nSMD"
        ]
    )
)

MAX_WEIGHTED_SMD <- max(
    abs(
        weighted_smd_values
    ),
    na.rm = TRUE
)

N_WEIGHTED_GT_010 <- sum(
    abs(
        weighted_smd_values
    ) >
    0.10,
    na.rm = TRUE
)


cat("\n")
cat("============================================================\n")
cat("TABLE 1 COMPLETE\n")
cat("============================================================\n")

cat(
    "N = ",
    nrow(WT),
    " | limited = ",
    sum(limited_mask),
    " | colectomy = ",
    sum(colectomy_mask),
    "\n",
    sep = ""
)

cat(
    sprintf(
        "Maximum weighted |SMD| = %.3f\n",
        MAX_WEIGHTED_SMD
    )
)

cat(
    "Weighted |SMD| > 0.10: ",
    N_WEIGHTED_GT_010,
    "\n",
    sep = ""
)

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Standard: journal-ready three-line table, no vertical rules.\n")
cat("Target layout: FINAL LOCK one-page landscape table with footnote on page 1.\n")
cat("============================================================\n")
