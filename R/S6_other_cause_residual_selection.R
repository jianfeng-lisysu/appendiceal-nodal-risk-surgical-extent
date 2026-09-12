## =============================================================================
## A6 - SUPPLEMENTARY TABLE S5 ONLY
## Negative-control analysis of 5-year other-cause mortality
## FULL-PIPELINE BOOTSTRAP B=1000 — SYNTAX-CHECKED VERSION
##
## Root:
##   <project root>
##
## INPUTS
##   02_data/analytic_cohort_A.csv
##   02_data/analytic_cohort_A_weighted.csv
##   02_data/Aline_frozen_pipeline_v2.RData
##
## ANALYSIS
##   Outcome:
##     5-year cumulative incidence of other-cause death
##     with cancer-specific death as the competing event
##
##   Exposure contrast:
##     oncologic colectomy - limited resection
##
##   Strata:
##     Q1-Q5 of predicted nodal risk
##     bottom predicted-risk decile
##     predicted nodal risk <5%
##
##   Bootstrap:
##     B = 1000
##     complete frozen pipeline re-run in every bootstrap replicate
##     percentile-based risk cut points re-estimated in every replicate
##     <5% threshold fixed
##
## OUTPUTS
##   03_results/06_negative_control_othercause_FINAL_B1000.csv
##   03_results/06_negative_control_othercause_bootstrap_FINAL_B1000.csv
##   03_results/06_negative_control_othercause_FINAL_B1000.txt
##   05_tables/Supplementary_Table_S5_Negative_Control_FINAL.docx
##
## INTERPRETATION
##   Negative RD = lower other-cause mortality associated with oncologic colectomy.
##   This is used qualitatively as a residual-selection diagnostic only.
##   It is NOT used for numerical bias subtraction or causal correction.
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

RESULT_DIR <- file.path(
    ROOT,
    "03_results"
)

TABLE_DIR <- file.path(
    ROOT,
    "05_tables"
)

dir.create(
    RESULT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
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

WT_FILE <- locate_file(
    "analytic_cohort_A_weighted.csv"
)

PIPE_FILE <- locate_file(
    "Aline_frozen_pipeline_v2.RData"
)


## =============================================================================
## 3. Load frozen pipeline and data
## =============================================================================

load(
    PIPE_FILE
)

required_objects <- c(
    "prep",
    "pipeline",
    "RISK_FORMULA",
    "PS_FORMULA"
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

WT <- read.csv(
    WT_FILE,
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
        "Frozen RAW cohort N should be 6951."
    )
}

if (
    nrow(
        WT
    ) !=
    6951
) {

    stop(
        "Frozen weighted cohort N should be 6951."
    )
}


if (
    sum(
        RAW$exposure ==
        "limited"
    ) !=
    2877 ||
    sum(
        RAW$exposure ==
        "colectomy"
    ) !=
    4074
) {

    stop(
        "Frozen treatment-arm fingerprint failed."
    )
}


## =============================================================================
## 5. Re-run original sample through frozen pipeline
## =============================================================================

d0 <- pipeline(
    RAW
)


required_cols <- c(
    "treat",
    "ps",
    "ow",
    "nodal_risk",
    "os_months",
    "css_event",
    "other_death"
)

missing_cols <- setdiff(
    required_cols,
    names(
        d0
    )
)

if (
    length(
        missing_cols
    ) >
    0
) {

    stop(
        "pipeline() output missing: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}


## =============================================================================
## 6. Exact pipeline fingerprint against frozen weighted file
## =============================================================================

if (
    "patient_id" %in%
    names(
        d0
    ) &&
    "patient_id" %in%
    names(
        WT
    )
) {

    if (
        !identical(
            as.character(
                d0$patient_id
            ),
            as.character(
                WT$patient_id
            )
        )
    ) {

        stop(
            "patient_id ordering differs from frozen weighted file."
        )
    }
}


RISK_FP <- max(
    abs(
        d0$nodal_risk -
        WT$nodal_risk
    ),
    na.rm = TRUE
)

PS_FP <- max(
    abs(
        d0$ps -
        WT$ps
    ),
    na.rm = TRUE
)

OW_FP <- max(
    abs(
        d0$ow -
        WT$ow
    ),
    na.rm = TRUE
)


if (
    RISK_FP >
    1e-8 ||
    PS_FP >
    1e-8 ||
    OW_FP >
    1e-8
) {

    stop(
        sprintf(
            paste0(
                "Frozen pipeline fingerprint failed: ",
                "risk=%.12g, PS=%.12g, OW=%.12g"
            ),
            RISK_FP,
            PS_FP,
            OW_FP
        )
    )
}


## =============================================================================
## 7. Weighted Aalen-Johansen estimator
##
## event1 = event of interest
## event2 = competing event
## start = -1 includes month-0 records, consistent with final primary analysis
## tau = 60 months after diagnosis
## =============================================================================

aj_w <- function(
    time,
    event1,
    event2,
    wt,
    start = -1,
    tau = 60
) {

    ok <- (
        !is.na(
            time
        ) &
        !is.na(
            event1
        ) &
        !is.na(
            event2
        ) &
        !is.na(
            wt
        ) &
        is.finite(
            wt
        ) &
        wt >=
        0
    )

    time <- time[ok]
    event1 <- event1[ok]
    event2 <- event2[ok]
    wt <- wt[ok]

    keep <- time >
        start

    time <- time[keep]
    event1 <- event1[keep]
    event2 <- event2[keep]
    wt <- wt[keep]

    if (
        length(
            time
        ) ==
        0 ||
        sum(
            wt
        ) <=
        0
    ) {

        return(
            NA_real_
        )
    }


    ut <- sort(
        unique(
            time
        )
    )

    S <- 1
    cif <- 0
    at_risk <- sum(
        wt
    )


    for (
        tt in ut
    ) {

        idx <- time ==
            tt

        d1 <- sum(
            wt[idx] *
            event1[idx]
        )

        d2 <- sum(
            wt[idx] *
            event2[idx]
        )

        dw <- sum(
            wt[idx]
        )


        if (
            tt <=
            tau &&
            at_risk >
            0
        ) {

            cif <- cif +
                S *
                d1 /
                at_risk

            S <- S *
                (
                    1 -
                    (
                        d1 +
                        d2
                    ) /
                    at_risk
                )
        }


        at_risk <- at_risk -
            dw

        if (
            tt >
            tau
        ) {
            break
        }
    }


    as.numeric(
        cif
    )
}


## =============================================================================
## 8. Other-cause mortality estimator within one subset
## =============================================================================

estimate_othercause <- function(
    g
) {

    if (
        nrow(
            g
        ) <
        20
    ) {

        return(
            c(
                limited = NA_real_,
                colectomy = NA_real_,
                rd = NA_real_
            )
        )
    }


    if (
        sum(
            g$treat ==
            0
        ) <
        10 ||
        sum(
            g$treat ==
            1
        ) <
        10
    ) {

        return(
            c(
                limited = NA_real_,
                colectomy = NA_real_,
                rd = NA_real_
            )
        )
    }


    ## Other-cause death is the event of interest.
    ## Cancer-specific death is the competing event.
    event_other <- as.integer(
        g$other_death ==
        1
    )

    event_cancer <- as.integer(
        g$css_event ==
        1
    )


    lim <- g$treat ==
        0

    col <- g$treat ==
        1


    cif_limited <- aj_w(
        time = g$os_months[lim],
        event1 = event_other[lim],
        event2 = event_cancer[lim],
        wt = g$ow[lim],
        start = -1,
        tau = 60
    )


    cif_colectomy <- aj_w(
        time = g$os_months[col],
        event1 = event_other[col],
        event2 = event_cancer[col],
        wt = g$ow[col],
        start = -1,
        tau = 60
    )


    c(
        limited = cif_limited,
        colectomy = cif_colectomy,
        rd = cif_colectomy -
            cif_limited
    )
}


## =============================================================================
## 9. Build the seven prespecified risk strata
## =============================================================================

make_strata <- function(
    d
) {

    q20 <- as.numeric(
        quantile(
            d$nodal_risk,
            probs = c(
                0.20,
                0.40,
                0.60,
                0.80
            ),
            na.rm = TRUE,
            names = FALSE
        )
    )


    p10 <- as.numeric(
        quantile(
            d$nodal_risk,
            probs = 0.10,
            na.rm = TRUE,
            names = FALSE
        )
    )


    if (
        any(
            !is.finite(
                q20
            )
        ) ||
        !is.finite(
            p10
        )
    ) {

        stop(
            "Non-finite risk cut point."
        )
    }


    if (
        any(
            diff(
                q20
            ) <=
            0
        )
    ) {

        stop(
            "Non-unique risk-quintile cut points."
        )
    }


    quintile <- cut(
        d$nodal_risk,
        breaks = c(
            -Inf,
            q20,
            Inf
        ),
        labels = paste0(
            "Q",
            1:5
        ),
        include.lowest = TRUE,
        right = TRUE
    )


    list(
        Q1 = quintile ==
            "Q1",
        Q2 = quintile ==
            "Q2",
        Q3 = quintile ==
            "Q3",
        Q4 = quintile ==
            "Q4",
        Q5 = quintile ==
            "Q5",
        bottom_decile = d$nodal_risk <=
            p10,
        risk_lt_5pct = d$nodal_risk <
            0.05,
        q20 = q20,
        p10 = p10
    )
}


## =============================================================================
## 10. Estimate all seven strata
## =============================================================================

estimate_all7 <- function(
    d
) {

    st <- make_strata(
        d
    )


    ids <- c(
        "Q1",
        "Q2",
        "Q3",
        "Q4",
        "Q5",
        "bottom_decile",
        "risk_lt_5pct"
    )


    out <- matrix(
        NA_real_,
        nrow = length(
            ids
        ),
        ncol = 3,
        dimnames = list(
            ids,
            c(
                "limited",
                "colectomy",
                "rd"
            )
        )
    )


    for (
        i in seq_along(
            ids
        )
    ) {

        idx <- st[[ids[i]]]


        out[
            i,
            ] <- estimate_othercause(
                d[
                    idx,
                    ,
                    drop = FALSE
                ]
            )
    }


    list(
        est = out,
        q20 = st$q20,
        p10 = st$p10
    )
}


## =============================================================================
## 11. Original-sample point estimates
## =============================================================================

point_obj <- estimate_all7(
    d0
)

point_est <- point_obj$est

Q20_ORIGINAL <- point_obj$q20

P10_ORIGINAL <- point_obj$p10


## =============================================================================
## 12. Full-pipeline bootstrap B=1000
## =============================================================================

B <- 1000

set.seed(
    20260830
)


ids <- rownames(
    point_est
)


boot_rd <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(
        ids
    ),
    dimnames = list(
        NULL,
        ids
    )
)


t0 <- Sys.time()


for (
    b in seq_len(
        B
    )
) {

    idx_b <- sample(
        seq_len(
            nrow(
                RAW
            )
        ),
        size = nrow(
            RAW
        ),
        replace = TRUE
    )


    db <- RAW[
        idx_b,
        ,
        drop = FALSE
    ]


    boot_vals <- tryCatch(
        {

            dx <- pipeline(
                db
            )


            z <- estimate_all7(
                dx
            )


            z$est[
                ,
                "rd"
            ]
        },
        error = function(
            e
        ) {

            rep(
                NA_real_,
                length(
                    ids
                )
            )
        }
    )


    boot_rd[
        b,
        ] <- boot_vals


    if (
        b %%
        100 ==
        0
    ) {

        elapsed <- as.numeric(
            difftime(
                Sys.time(),
                t0,
                units = "mins"
            )
        )


        cat(
            "Bootstrap ",
            b,
            "/",
            B,
            " | elapsed ",
            sprintf(
                "%.1f",
                elapsed
            ),
            " min\n",
            sep = ""
        )
    }
}


## =============================================================================
## 13. Bootstrap percentile intervals
## =============================================================================

ci_lo <- rep(
    NA_real_,
    length(
        ids
    )
)

ci_hi <- rep(
    NA_real_,
    length(
        ids
    )
)

effective_B <- rep(
    0L,
    length(
        ids
    )
)


for (
    i in seq_along(
        ids
    )
) {

    v <- boot_rd[
        ,
        i
    ]

    v <- v[
        is.finite(
            v
        )
    ]


    effective_B[i] <- length(
        v
    )


    if (
        length(
            v
        ) >=
        100
    ) {

        ci_lo[i] <- as.numeric(
            quantile(
                v,
                0.025,
                na.rm = TRUE,
                names = FALSE
            )
        )

        ci_hi[i] <- as.numeric(
            quantile(
                v,
                0.975,
                na.rm = TRUE,
                names = FALSE
            )
        )
    }
}


if (
    any(
        effective_B <
        900
    )
) {

    warning(
        paste0(
            "Some negative-control strata have effective bootstrap B <900: ",
            paste(
                ids[
                    effective_B <
                    900
                ],
                collapse = ", "
            )
        )
    )
}


## =============================================================================
## 14. Final analysis result object
## =============================================================================

result <- data.frame(
    stratum = ids,
    limited_cif = point_est[
        ,
        "limited"
    ],
    colectomy_cif = point_est[
        ,
        "colectomy"
    ],
    rd = point_est[
        ,
        "rd"
    ],
    lo = ci_lo,
    hi = ci_hi,
    effective_B = effective_B,
    stringsAsFactors = FALSE
)


## =============================================================================
## 15. Publication labels
## =============================================================================

label_map <- c(
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


## =============================================================================
## 16. Save analysis CSV + bootstrap CSV
## =============================================================================

RESULT_CSV <- file.path(
    RESULT_DIR,
    "06_negative_control_othercause_FINAL_B1000.csv"
)

BOOT_CSV <- file.path(
    RESULT_DIR,
    "06_negative_control_othercause_bootstrap_FINAL_B1000.csv"
)


write.csv(
    result,
    RESULT_CSV,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


boot_out <- data.frame(
    replicate = seq_len(
        B
    ),
    boot_rd,
    check.names = FALSE
)

write.csv(
    boot_out,
    BOOT_CSV,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 17. Build Supplementary Table S5
## =============================================================================

table_s5 <- data.frame(
    `Predicted nodal-risk stratum` =
        unname(
            label_map[
                result$stratum
            ]
        ),
    `Limited resection, %` =
        sprintf(
            "%.1f",
            100 *
            result$limited_cif
        ),
    `Oncologic colectomy, %` =
        sprintf(
            "%.1f",
            100 *
            result$colectomy_cif
        ),
    `RD, pp` =
        sprintf(
            "%+.1f",
            100 *
            result$rd
        ),
    `95% CI, pp` =
        sprintf(
            "%+.1f to %+.1f",
            100 *
            result$lo,
            100 *
            result$hi
        ),
    `Effective bootstrap B` =
        result$effective_B,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


## =============================================================================
## 18. Journal-ready three-line table
## =============================================================================

ft <- flextable(
    table_s5
)

ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
)

ft <- fontsize(
    ft,
    size = 8.4,
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
    j = 2:ncol(
        table_s5
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
    padding.top = 1.3,
    padding.bottom = 1.3,
    padding.left = 1.6,
    padding.right = 1.6,
    part = "all"
)


ft <- padding(
    ft,
    i = 6,
    padding.top = 4.0,
    part = "body"
)


ft <- width(
    ft,
    j = 1,
    width = 2.55
)

ft <- width(
    ft,
    j = 2,
    width = 1.25
)

ft <- width(
    ft,
    j = 3,
    width = 1.40
)

ft <- width(
    ft,
    j = 4,
    width = 0.75
)

ft <- width(
    ft,
    j = 5,
    width = 1.35
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
## 19. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S5_Negative_Control_FINAL.docx"
)

doc <- read_docx()


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S5. Other-cause mortality as a residual-selection diagnostic by predicted nodal risk",
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
    "Values in the two treatment columns are overlap-weighted 5-year cumulative incidences of other-cause death, ",
    "with cancer-specific death treated as the competing event. RD is oncologic colectomy minus limited resection; ",
    "negative values indicate lower other-cause mortality associated with oncologic colectomy. ",
    "The prespecified risk-model and propensity-score pipeline was re-estimated within each of 1,000 bootstrap replicates. ",
    "Percentile-based risk cut points were re-estimated in each bootstrap replicate, whereas the <5% threshold remained fixed. ",
    "This analysis was used qualitatively as a residual-selection diagnostic and was not used for numerical bias subtraction ",
    "or causal correction. Because surgery can plausibly influence noncancer mortality, the analysis should not be interpreted ",
    "as a strict negative control satisfying a no-effect assumption."
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
## 20. Human-readable report
## =============================================================================

REPORT_FILE <- file.path(
    RESULT_DIR,
    "06_negative_control_othercause_FINAL_B1000.txt"
)


report <- c(
    "== FINAL other-cause mortality residual-selection diagnostic ==",
    paste0(
        "RAW N=",
        nrow(
            RAW
        ),
        "; bootstrap B=",
        B,
        "."
    ),
    paste0(
        "RD = oncologic colectomy - limited resection in 5-year other-cause mortality."
    ),
    paste0(
        "Negative RD indicates lower other-cause mortality associated with oncologic colectomy."
    ),
    paste0(
        "Interpretation: qualitative residual-selection diagnostic only; not a strict no-effect negative control and no numerical bias subtraction."
    ),
    ""
)


for (
    i in seq_len(
        nrow(
            result
        )
    )
) {

    report <- c(
        report,
        sprintf(
            "%-30s: limited %.1f%% | colectomy %.1f%% | RD %+.1f pp (95%% CI %+.1f to %+.1f; effective B=%d)",
            label_map[
                result$stratum[i]
            ],
            100 *
            result$limited_cif[i],
            100 *
            result$colectomy_cif[i],
            100 *
            result$rd[i],
            100 *
            result$lo[i],
            100 *
            result$hi[i],
            result$effective_B[i]
        )
    )
}


report <- c(
    report,
    "",
    "Original-sample risk cut points:",
    paste0(
        "P10 = ",
        sprintf(
            "%.2f%%",
            100 *
            P10_ORIGINAL
        )
    ),
    paste0(
        "Quintile boundaries = ",
        paste(
            sprintf(
                "%.2f%%",
                100 *
                Q20_ORIGINAL
            ),
            collapse = ", "
        )
    ),
    "",
    paste0(
        "Frozen pipeline fingerprints: risk=",
        format(
            RISK_FP,
            scientific = TRUE,
            digits = 3
        ),
        "; PS=",
        format(
            PS_FP,
            scientific = TRUE,
            digits = 3
        ),
        "; OW=",
        format(
            OW_FP,
            scientific = TRUE,
            digits = 3
        ),
        "."
    )
)


writeLines(
    report,
    REPORT_FILE,
    useBytes = TRUE
)


## =============================================================================
## 21. Console
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S5 FINAL LOCK VERSION COMPLETE\n")
cat("============================================================\n")

cat(
    sprintf(
        "Pipeline fingerprint: risk %.12g | PS %.12g | OW %.12g\n",
        RISK_FP,
        PS_FP,
        OW_FP
    )
)

cat("\nFinal negative-control results:\n")

for (
    i in seq_len(
        nrow(
            result
        )
    )
) {

    cat(
        sprintf(
            "%-30s | limited %.1f%% | colectomy %.1f%% | RD %+.1f pp | 95%% CI %+.1f to %+.1f | B=%d\n",
            label_map[
                result$stratum[i]
            ],
            100 *
            result$limited_cif[i],
            100 *
            result$colectomy_cif[i],
            100 *
            result$rd[i],
            100 *
            result$lo[i],
            100 *
            result$hi[i],
            result$effective_B[i]
        )
    )
}

cat("\nDOCX  : ", DOCX_FILE, "\n", sep = "")
cat("Result: ", RESULT_CSV, "\n", sep = "")
cat("Boot  : ", BOOT_CSV, "\n", sep = "")
cat("Report: ", REPORT_FILE, "\n", sep = "")

cat("\n")
cat("Full-pipeline bootstrap B=1000 completed.\n")
cat("No numerical bias subtraction was performed.\n")
cat("============================================================\n")
