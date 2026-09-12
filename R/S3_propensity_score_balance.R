## =============================================================================
## A6 - SUPPLEMENTARY TABLE S3 ONLY
## Propensity-score overlap and covariate balance after overlap weighting
## JOURNAL-READY TWO-PANEL THREE-LINE TABLE — COMPATIBILITY FIX
##
## Root:
##   <project root>
##
## INPUT
##   02_data/analytic_cohort_A_weighted.csv
##
## OUTPUTS
##   05_tables/Supplementary_Table_S3_PS_Balance_FINAL.docx
##   05_tables/Supplementary_Table_S3_PS_Diagnostics_FINAL.csv
##   05_tables/Supplementary_Table_S3_Covariate_Balance_FINAL.csv
##
## FROZEN PRIMARY COHORT
##   N = 6,951
##   Limited resection = 2,877
##   Oncologic colectomy = 4,074
##
## FROZEN WEIGHTING TARGETS
##   ESS colectomy ≈ 3,780
##   ESS limited ≈ 2,735
##   Maximum overlap weight ≈ 0.876
##
## IMPORTANT
##   This script does not refit or alter the propensity-score model.
##   It audits the already-frozen PS and overlap weights.
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
## 2. Locate and read weighted cohort
## =============================================================================

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


WT <- read.csv(
    WT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


## =============================================================================
## 3. Frozen cohort fingerprints
## =============================================================================

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


WT$treat <- as.integer(
    WT$exposure ==
    "colectomy"
)


N_LIMITED <- sum(
    WT$treat ==
    0
)

N_COLECTOMY <- sum(
    WT$treat ==
    1
)


if (
    N_LIMITED !=
    2877 ||
    N_COLECTOMY !=
    4074
) {

    stop(
        "Frozen treatment-arm fingerprint failed."
    )
}


required_cols <- c(
    "ps",
    "ow",
    "age",
    "female",
    "year",
    "size_i",
    "size_miss",
    "logsize",
    "T",
    "grade_f",
    "hist_group",
    "stage_sum",
    "race4",
    "income_hi",
    "metro",
    "married",
    "nodal_risk"
)

missing_cols <- setdiff(
    required_cols,
    names(
        WT
    )
)

if (
    length(
        missing_cols
    ) >
    0
) {

    stop(
        "Weighted cohort missing columns: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


## =============================================================================
## 4. Validate overlap-weight identity
##
## PS = probability of oncologic colectomy.
## Overlap weight:
##   colectomy: 1 - PS
##   limited:    PS
## =============================================================================

ow_expected <- ifelse(
    WT$treat ==
    1,
    1 -
    WT$ps,
    WT$ps
)


OW_IDENTITY_ERROR <- max(
    abs(
        WT$ow -
        ow_expected
    ),
    na.rm = TRUE
)


if (
    OW_IDENTITY_ERROR >
    1e-10
) {

    stop(
        sprintf(
            "Overlap-weight identity failed: max error %.12g",
            OW_IDENTITY_ERROR
        )
    )
}


## =============================================================================
## 5. Helpers
## =============================================================================

ess_fun <- function(
    w
) {

    w <- w[
        !is.na(
            w
        ) &
        is.finite(
            w
        )
    ]

    if (
        length(
            w
        ) ==
        0 ||
        sum(
            w^2
        ) ==
        0
    ) {

        return(
            NA_real_
        )
    }

    sum(
        w
    )^2 /
    sum(
        w^2
    )
}


weighted_mean_local <- function(
    x,
    w
) {

    ok <- (
        !is.na(
            x
        ) &
        !is.na(
            w
        ) &
        is.finite(
            w
        )
    )

    if (
        sum(
            ok
        ) ==
        0
    ) {

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
        !is.na(
            x
        ) &
        !is.na(
            w
        ) &
        is.finite(
            w
        )
    )

    x <- x[ok]
    w <- w[ok]

    if (
        length(
            x
        ) <
        2 ||
        sum(
            w
        ) <=
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
    sum(
        w
    )
}


smd_fun <- function(
    x,
    treat,
    w
) {

    ok <- (
        !is.na(
            x
        ) &
        !is.na(
            treat
        ) &
        !is.na(
            w
        ) &
        is.finite(
            w
        )
    )

    x <- x[ok]
    treat <- treat[ok]
    w <- w[ok]

    if (
        sum(
            treat ==
            1
        ) ==
        0 ||
        sum(
            treat ==
            0
        ) ==
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
        !is.finite(
            den
        ) ||
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


fmt3 <- function(
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
        abs(
            x
        ) <
        0.0005
    ) {
        x <- 0
    }

    sprintf(
        "%.3f",
        x
    )
}


fmt_range <- function(
    lo,
    hi
) {

    sprintf(
        "%.3f–%.3f",
        lo,
        hi
    )
}


## =============================================================================
## 6. Panel A — PS / weight diagnostics
## =============================================================================

ps_lim <- WT$ps[
    WT$treat ==
    0
]

ps_col <- WT$ps[
    WT$treat ==
    1
]

ow_lim <- WT$ow[
    WT$treat ==
    0
]

ow_col <- WT$ow[
    WT$treat ==
    1
]


ESS_LIMITED <- ess_fun(
    ow_lim
)

ESS_COLECTOMY <- ess_fun(
    ow_col
)

MAX_WEIGHT <- max(
    WT$ow,
    na.rm = TRUE
)


## Hard fingerprints with sensible tolerance.
if (
    abs(
        ESS_LIMITED -
        2735
    ) >
    5
) {

    stop(
        sprintf(
            "Limited ESS fingerprint failed: %.2f",
            ESS_LIMITED
        )
    )
}


if (
    abs(
        ESS_COLECTOMY -
        3780
    ) >
    5
) {

    stop(
        sprintf(
            "Colectomy ESS fingerprint failed: %.2f",
            ESS_COLECTOMY
        )
    )
}


if (
    abs(
        MAX_WEIGHT -
        0.876
    ) >
    0.01
) {

    stop(
        sprintf(
            "Max-weight fingerprint failed: %.4f",
            MAX_WEIGHT
        )
    )
}


panel_a <- data.frame(
    `Diagnostic` = c(
        "Observed sample size",
        "Propensity score, P1–P99",
        "Propensity score, median",
        "Overlap weight, mean",
        "Overlap weight, maximum",
        "Effective sample size"
    ),
    `Limited resection` = c(
        format(
            N_LIMITED,
            big.mark = ","
        ),
        fmt_range(
            as.numeric(
                quantile(
                    ps_lim,
                    0.01,
                    na.rm = TRUE
                )
            ),
            as.numeric(
                quantile(
                    ps_lim,
                    0.99,
                    na.rm = TRUE
                )
            )
        ),
        fmt3(
            median(
                ps_lim,
                na.rm = TRUE
            )
        ),
        fmt3(
            mean(
                ow_lim,
                na.rm = TRUE
            )
        ),
        fmt3(
            max(
                ow_lim,
                na.rm = TRUE
            )
        ),
        sprintf(
            "%.0f",
            ESS_LIMITED
        )
    ),
    `Oncologic colectomy` = c(
        format(
            N_COLECTOMY,
            big.mark = ","
        ),
        fmt_range(
            as.numeric(
                quantile(
                    ps_col,
                    0.01,
                    na.rm = TRUE
                )
            ),
            as.numeric(
                quantile(
                    ps_col,
                    0.99,
                    na.rm = TRUE
                )
            )
        ),
        fmt3(
            median(
                ps_col,
                na.rm = TRUE
            )
        ),
        fmt3(
            mean(
                ow_col,
                na.rm = TRUE
            )
        ),
        fmt3(
            max(
                ow_col,
                na.rm = TRUE
            )
        ),
        sprintf(
            "%.0f",
            ESS_COLECTOMY
        )
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 7. Panel B — balance components
##
## Include all core clinical balance components plus social/race variables
## that were included in the final PS specification, and nodal risk as an
## additional prognostic-balance diagnostic.
## =============================================================================

grade_chr <- as.character(
    WT$grade_f
)

race_chr <- as.character(
    WT$race4
)


balance_matrix <- data.frame(
    `Age, years` =
        WT$age,
    `Age squared` =
        WT$age^2,
    `Female sex` =
        WT$female,
    `Year of diagnosis` =
        WT$year,
    `Tumor size, mm (imputed)` =
        WT$size_i,
    `Log-transformed tumor size` =
        WT$logsize,
    `Tumor size missing` =
        WT$size_miss,
    `T2` =
        as.numeric(
            WT$T ==
            "T2"
        ),
    `T3` =
        as.numeric(
            WT$T ==
            "T3"
        ),
    `T4` =
        as.numeric(
            WT$T ==
            "T4"
        ),
    `Grade 1` =
        as.numeric(
            grade_chr ==
            "1"
        ),
    `Grade 2` =
        as.numeric(
            grade_chr ==
            "2"
        ),
    `Grade 3` =
        as.numeric(
            grade_chr ==
            "3"
        ),
    `Grade 4` =
        as.numeric(
            grade_chr ==
            "4"
        ),
    `Grade unknown` =
        as.numeric(
            grade_chr ==
            "0"
        ),
    `Nonmucinous adenocarcinoma` =
        as.numeric(
            WT$hist_group ==
            "nonmucinous"
        ),
    `Mucinous adenocarcinoma` =
        as.numeric(
            WT$hist_group ==
            "mucinous"
        ),
    `Goblet cell adenocarcinoma` =
        as.numeric(
            WT$hist_group ==
            "GCA"
        ),
    `Signet-ring cell carcinoma` =
        as.numeric(
            WT$hist_group ==
            "SRCC"
        ),
    `Regional stage` =
        as.numeric(
            WT$stage_sum ==
            "Regional"
        ),
    `Hispanic` =
        as.numeric(
            race_chr ==
            "Hispanic (All Races)"
        ),
    `Non-Hispanic Black` =
        as.numeric(
            race_chr ==
            "Non-Hispanic Black"
        ),
    `Non-Hispanic White` =
        as.numeric(
            race_chr ==
            "Non-Hispanic White"
        ),
    `Other/unknown race-ethnicity` =
        as.numeric(
            race_chr ==
            "Other/Unknown"
        ),
    `Higher-income category` =
        WT$income_hi,
    `Metropolitan residence` =
        WT$metro,
    `Married` =
        WT$married,
    `Predicted nodal risk` =
        WT$nodal_risk,
    check.names = FALSE
)


balance_rows <- vector(
    "list",
    ncol(
        balance_matrix
    )
)


for (
    j in seq_len(
        ncol(
            balance_matrix
        )
    )
) {

    x <- balance_matrix[
        ,
        j
    ]

    balance_rows[[j]] <- data.frame(
        Covariate =
            names(
                balance_matrix
            )[
                j
            ],
        Unweighted_SMD =
            smd_fun(
                x,
                WT$treat,
                rep(
                    1,
                    nrow(
                        WT
                    )
                )
            ),
        Weighted_SMD =
            smd_fun(
                x,
                WT$treat,
                WT$ow
            ),
        stringsAsFactors = FALSE
    )
}


balance_raw <- do.call(
    rbind,
    balance_rows
)


MAX_ABS_WEIGHTED_SMD <- max(
    abs(
        balance_raw$Weighted_SMD
    ),
    na.rm = TRUE
)

N_WEIGHTED_GT_010 <- sum(
    abs(
        balance_raw$Weighted_SMD
    ) >
    0.10,
    na.rm = TRUE
)


if (
    MAX_ABS_WEIGHTED_SMD >
    0.02
) {

    stop(
        sprintf(
            "Unexpected weighted imbalance: max |SMD| = %.4f",
            MAX_ABS_WEIGHTED_SMD
        )
    )
}


panel_b <- data.frame(
    Covariate =
        balance_raw$Covariate,
    `Unweighted SMD` =
        vapply(
            balance_raw$Unweighted_SMD,
            fmt3,
            character(
                1
            )
        ),
    `Overlap-weighted SMD` =
        vapply(
            balance_raw$Weighted_SMD,
            fmt3,
            character(
                1
            )
        ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 8. Save CSV audit copies
## =============================================================================

PANEL_A_CSV <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_PS_Diagnostics_FINAL.csv"
)

PANEL_B_CSV <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_Covariate_Balance_FINAL.csv"
)


write.csv(
    panel_a,
    PANEL_A_CSV,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)

write.csv(
    panel_b,
    PANEL_B_CSV,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 9. Three-line table helper
## =============================================================================

make_three_line <- function(
    dat,
    font_size = 8.2
) {

    ft <- flextable(
        dat
    )

    ft <- font(
        ft,
        fontname = "Times New Roman",
        part = "all"
    )

    ft <- fontsize(
        ft,
        size = font_size,
        part = "all"
    )

    ft <- fontsize(
        ft,
        size = font_size +
        0.2,
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

    if (
        ncol(
            dat
        ) >
        1
    ) {

        ft <- align(
            ft,
            j = 2:ncol(
                dat
            ),
            align = "center",
            part = "all"
        )
    }

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
        padding.top = 0.9,
        padding.bottom = 0.9,
        padding.left = 1.5,
        padding.right = 1.5,
        part = "all"
    )

    ft
}


## =============================================================================
## 10. Format Panel A
## =============================================================================

ft_a <- make_three_line(
    panel_a,
    font_size = 8.2
)

ft_a <- width(
    ft_a,
    j = 1,
    width = 2.60
)

ft_a <- width(
    ft_a,
    j = 2,
    width = 1.60
)

ft_a <- width(
    ft_a,
    j = 3,
    width = 1.70
)

ft_a <- set_table_properties(
    ft_a,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 11. Format Panel B
## =============================================================================

ft_b <- make_three_line(
    panel_b,
    font_size = 7.9
)

ft_b <- width(
    ft_b,
    j = 1,
    width = 3.30
)

ft_b <- width(
    ft_b,
    j = 2,
    width = 1.35
)

ft_b <- width(
    ft_b,
    j = 3,
    width = 1.55
)

ft_b <- set_table_properties(
    ft_b,
    layout = "fixed",
    width = 1
)


## =============================================================================
## 12. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_PS_Balance_FINAL.docx"
)


doc <- read_docx()


## Main title
doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S3. Propensity-score overlap and covariate balance after overlap weighting",
            fp_text(
                font.family = "Times New Roman",
                font.size = 9.5,
                bold = TRUE
            )
        )
    )
)


## Panel A title
doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Panel A. Propensity-score and weighting diagnostics",
            fp_text(
                font.family = "Times New Roman",
                font.size = 8.6,
                bold = TRUE
            )
        )
    )
)


doc <- body_add_flextable(
    doc,
    value = ft_a,
    align = "center"
)


## Panel B title
doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Panel B. Standardized mean differences before and after overlap weighting",
            fp_text(
                font.family = "Times New Roman",
                font.size = 8.6,
                bold = TRUE
            )
        )
    )
)


doc <- body_add_flextable(
    doc,
    value = ft_b,
    align = "center"
)


## Footnote
note_text <- paste0(
    "PS indicates propensity score; SMD, standardized mean difference. ",
    "The propensity score represents the modeled probability of oncologic colectomy. ",
    "Overlap weights were defined as 1−PS for oncologic colectomy and PS for limited resection. ",
    "Effective sample size was calculated as (sum of weights)^2 / sum of squared weights. ",
    "Absolute SMD <0.10 was prespecified as acceptable balance. ",
    "Panel B includes the modeled baseline covariates and predicted nodal risk as an additional prognostic-balance diagnostic."
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
## 13. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S3 COMPLETE\n")
cat("Compatibility note: T_f was not required; raw T was used for balance indicators.\n")
cat("============================================================\n")

cat(
    "N = ",
    nrow(
        WT
    ),
    " | limited = ",
    N_LIMITED,
    " | colectomy = ",
    N_COLECTOMY,
    "\n",
    sep = ""
)

cat(
    sprintf(
        "OW identity max error = %.12g\n",
        OW_IDENTITY_ERROR
    )
)

cat(
    sprintf(
        "ESS limited = %.1f | ESS colectomy = %.1f\n",
        ESS_LIMITED,
        ESS_COLECTOMY
    )
)

cat(
    sprintf(
        "Maximum overlap weight = %.4f\n",
        MAX_WEIGHT
    )
)

cat(
    sprintf(
        "Maximum weighted |SMD| = %.4f\n",
        MAX_ABS_WEIGHTED_SMD
    )
)

cat(
    "Weighted |SMD| >0.10 = ",
    N_WEIGHTED_GT_010,
    "\n",
    sep = ""
)

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("Panel A CSV: ", PANEL_A_CSV, "\n", sep = "")
cat("Panel B CSV: ", PANEL_B_CSV, "\n", sep = "")

cat("\n")
cat("Standard three-line tables; no vertical rules.\n")
cat("Frozen PS and overlap weights were audited, not refitted.\n")
cat("============================================================\n")
