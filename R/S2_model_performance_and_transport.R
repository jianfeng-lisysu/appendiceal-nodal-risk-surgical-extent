## =============================================================================
## A6 - SUPPLEMENTARY TABLE S2 ONLY
## Performance and transport validation of the final nodal-metastasis risk model
## JOURNAL-READY THREE-LINE TABLE — FINAL LOCK CANDIDATE
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
##   05_tables/Supplementary_Table_S2_Model_Performance_FINAL.docx
##   05_tables/Supplementary_Table_S2_Model_Performance_FINAL.csv
##
## FROZEN PERFORMANCE TARGETS
##   Development cohort:
##     n = 3,252
##     node-positive n = 738
##     full-model apparent AUC = 0.738
##     full-model apparent Brier = 0.1506
##     bootstrap-corrected AUC = 0.732
##     bootstrap-corrected calibration slope = 0.963
##     slope 95% CI = 0.865 to 1.063
##
##   T-stage-only benchmark:
##     AUC = 0.656
##     Brier = 0.1645
##
##   Limited-resection transport with >=1 examined LN:
##     n = 1,338
##     node-positive n = 263
##     AUC = 0.757
##     calibration slope = 1.067
##
##   Adequate-node transport with >=12 examined LN:
##     n = 879
##     AUC = 0.750
##     mean predicted risk = 0.225
##     observed node-positive proportion = 0.216
##     calibration slope = 1.014
##
## IMPORTANT
##   This script reproduces/validates the already-frozen model performance.
##   It does not change the model or cohort definitions.
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
## 3. Load frozen objects and data
## =============================================================================

load(
    PIPE_FILE
)

needed_objects <- c(
    "prep",
    "RISK_FORMULA"
)

missing_objects <- needed_objects[
    !vapply(
        needed_objects,
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


if (
    nrow(
        RAW
    ) !=
    6951
) {

    stop(
        "Frozen analytic cohort N should be 6951."
    )
}

if (
    nrow(
        RS
    ) !=
    6951
) {

    stop(
        "Frozen risk-scored cohort N should be 6951."
    )
}


## =============================================================================
## 4. Prepare data and refit final risk model
## =============================================================================

d <- prep(
    RAW
)


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


dev$node_positive <- as.integer(
    dev$n_positive >
    0
)

if (
    sum(
        dev$node_positive
    ) !=
    738
) {

    stop(
        "Development node-positive fingerprint failed."
    )
}


fit_full <- glm(
    RISK_FORMULA,
    data = dev,
    family = binomial()
)


pred_dev <- predict(
    fit_full,
    newdata = dev,
    type = "response"
)


pred_full <- predict(
    fit_full,
    newdata = d,
    type = "response"
)


## =============================================================================
## 5. Frozen nodal-risk fingerprint
## =============================================================================

if (
    length(
        pred_full
    ) !=
    nrow(
        RS
    )
) {

    stop(
        "Frozen risk-score fingerprint length mismatch."
    )
}


risk_fp <- max(
    abs(
        pred_full -
        RS$nodal_risk
    ),
    na.rm = TRUE
)


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
## 6. Metric helpers
## =============================================================================

auc_rank <- function(
    y,
    p
) {

    ok <- (
        !is.na(y) &
        !is.na(p)
    )

    y <- y[ok]
    p <- p[ok]

    n1 <- sum(
        y ==
        1
    )

    n0 <- sum(
        y ==
        0
    )

    if (
        n1 ==
        0 ||
        n0 ==
        0
    ) {

        return(
            NA_real_
        )
    }

    r <- rank(
        p,
        ties.method = "average"
    )

    (
        sum(
            r[
                y ==
                1
            ]
        ) -
        n1 *
        (
            n1 +
            1
        ) /
        2
    ) /
    (
        n1 *
        n0
    )
}


brier_score <- function(
    y,
    p
) {

    ok <- (
        !is.na(y) &
        !is.na(p)
    )

    mean(
        (
            y[ok] -
            p[ok]
        )^2
    )
}


calibration_slope <- function(
    y,
    p
) {

    ok <- (
        !is.na(y) &
        !is.na(p) &
        p >
        0 &
        p <
        1
    )

    fit <- glm(
        y[ok] ~
        qlogis(
            p[ok]
        ),
        family = binomial()
    )

    unname(
        coef(
            fit
        )[2]
    )
}


calibration_intercept <- function(
    y,
    p
) {

    ok <- (
        !is.na(y) &
        !is.na(p) &
        p >
        0 &
        p <
        1
    )

    fit <- glm(
        y[ok] ~
        offset(
            qlogis(
                p[ok]
            )
        ),
        family = binomial()
    )

    unname(
        coef(
            fit
        )[1]
    )
}


## =============================================================================
## 7. Development apparent performance
## =============================================================================

AUC_DEV_APP <- auc_rank(
    dev$node_positive,
    pred_dev
)

BRIER_DEV_APP <- brier_score(
    dev$node_positive,
    pred_dev
)

SLOPE_DEV_APP <- calibration_slope(
    dev$node_positive,
    pred_dev
)

INTERCEPT_DEV_APP <- calibration_intercept(
    dev$node_positive,
    pred_dev
)


## =============================================================================
## 8. T-stage-only benchmark
##
## Use the same development cohort and T-stage factor from the frozen prepared data.
## =============================================================================

fit_t_only <- glm(
    node_positive ~
    T_f,
    data = dev,
    family = binomial()
)

pred_t_only <- predict(
    fit_t_only,
    newdata = dev,
    type = "response"
)

AUC_T_ONLY <- auc_rank(
    dev$node_positive,
    pred_t_only
)

BRIER_T_ONLY <- brier_score(
    dev$node_positive,
    pred_t_only
)


## =============================================================================
## 9. Limited-resection transport cohorts
## =============================================================================

transport_any_mask <- (
    d$exposure ==
    "limited" &
    !is.na(
        d$ln_examined
    ) &
    d$ln_examined >=
    1 &
    !is.na(
        d$n_positive
    )
)

transport_any <- d[
    transport_any_mask,
    ,
    drop = FALSE
]

transport_any$node_positive <- as.integer(
    transport_any$n_positive >
    0
)

transport_any$pred <- predict(
    fit_full,
    newdata = transport_any,
    type = "response"
)


if (
    nrow(
        transport_any
    ) !=
    1338
) {

    stop(
        "Transport >=1 LN fingerprint failed: expected n=1338; current n=",
        nrow(
            transport_any
        )
    )
}


if (
    sum(
        transport_any$node_positive
    ) !=
    263
) {

    stop(
        "Transport >=1 LN node-positive fingerprint failed."
    )
}


AUC_TRANSPORT_ANY <- auc_rank(
    transport_any$node_positive,
    transport_any$pred
)

BRIER_TRANSPORT_ANY <- brier_score(
    transport_any$node_positive,
    transport_any$pred
)

SLOPE_TRANSPORT_ANY <- calibration_slope(
    transport_any$node_positive,
    transport_any$pred
)

INTERCEPT_TRANSPORT_ANY <- calibration_intercept(
    transport_any$node_positive,
    transport_any$pred
)

PRED_MEAN_TRANSPORT_ANY <- mean(
    transport_any$pred,
    na.rm = TRUE
)

OBS_MEAN_TRANSPORT_ANY <- mean(
    transport_any$node_positive,
    na.rm = TRUE
)


## Adequate-node transport >=12 LN
transport_adequate <- transport_any[
    transport_any$ln_examined >=
    12,
    ,
    drop = FALSE
]


if (
    nrow(
        transport_adequate
    ) !=
    879
) {

    stop(
        "Adequate-node transport fingerprint failed: expected n=879; current n=",
        nrow(
            transport_adequate
        )
    )
}


AUC_TRANSPORT_ADEQ <- auc_rank(
    transport_adequate$node_positive,
    transport_adequate$pred
)

BRIER_TRANSPORT_ADEQ <- brier_score(
    transport_adequate$node_positive,
    transport_adequate$pred
)

SLOPE_TRANSPORT_ADEQ <- calibration_slope(
    transport_adequate$node_positive,
    transport_adequate$pred
)

INTERCEPT_TRANSPORT_ADEQ <- calibration_intercept(
    transport_adequate$node_positive,
    transport_adequate$pred
)

PRED_MEAN_TRANSPORT_ADEQ <- mean(
    transport_adequate$pred,
    na.rm = TRUE
)

OBS_MEAN_TRANSPORT_ADEQ <- mean(
    transport_adequate$node_positive,
    na.rm = TRUE
)

NODE_POS_ADEQ <- sum(
    transport_adequate$node_positive,
    na.rm = TRUE
)


## =============================================================================
## 10. Frozen bootstrap-corrected development metrics
##
## These values come from the already-frozen B=500 full development bootstrap.
## They are reported, not recomputed here, to avoid changing the frozen analysis.
## =============================================================================

AUC_DEV_BOOT <- 0.732
SLOPE_DEV_BOOT <- 0.963
SLOPE_DEV_BOOT_LO <- 0.865
SLOPE_DEV_BOOT_HI <- 1.063


## =============================================================================
## 11. Hard metric fingerprints
## =============================================================================

check_close <- function(
    observed,
    expected,
    tolerance,
    label
) {

    if (
        !is.finite(
            observed
        ) ||
        abs(
            observed -
            expected
        ) >
        tolerance
    ) {

        stop(
            sprintf(
                "%s fingerprint failed: observed %.6f, expected %.6f",
                label,
                observed,
                expected
            )
        )
    }
}


check_close(
    AUC_DEV_APP,
    0.738,
    0.003,
    "Development apparent AUC"
)

check_close(
    BRIER_DEV_APP,
    0.1506,
    0.002,
    "Development apparent Brier"
)

check_close(
    AUC_T_ONLY,
    0.656,
    0.005,
    "T-stage-only AUC"
)

check_close(
    BRIER_T_ONLY,
    0.1645,
    0.002,
    "T-stage-only Brier"
)

check_close(
    AUC_TRANSPORT_ANY,
    0.757,
    0.005,
    "Transport >=1 LN AUC"
)

check_close(
    SLOPE_TRANSPORT_ANY,
    1.067,
    0.03,
    "Transport >=1 LN slope"
)

check_close(
    AUC_TRANSPORT_ADEQ,
    0.750,
    0.005,
    "Adequate-node transport AUC"
)

check_close(
    PRED_MEAN_TRANSPORT_ADEQ,
    0.225,
    0.01,
    "Adequate-node transport mean predicted risk"
)

check_close(
    OBS_MEAN_TRANSPORT_ADEQ,
    0.216,
    0.01,
    "Adequate-node transport observed proportion"
)

check_close(
    SLOPE_TRANSPORT_ADEQ,
    1.014,
    0.03,
    "Adequate-node transport slope"
)


## =============================================================================
## 12. Formatting helpers
## =============================================================================

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

    sprintf(
        "%.3f",
        x
    )
}


fmt4 <- function(
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

    sprintf(
        "%.4f",
        x
    )
}


fmt_prop <- function(
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

    sprintf(
        "%.3f",
        x
    )
}


## =============================================================================
## 13. Build final S2 display table — compact journal version
## =============================================================================

clean_zero <- function(x) {
    if (!is.finite(x)) return("—")
    if (abs(x) < 0.0005) x <- 0
    sprintf("%.3f", x)
}

fmt_pred_obs <- function(pred, obs) {
    if (!is.finite(pred) || !is.finite(obs)) return("—")
    sprintf("%.3f / %.3f", pred, obs)
}

table_s2 <- data.frame(
    `Evaluation set` = c(
        "Development: T-stage-only benchmark",
        "Development: final model (apparent)",
        "Development: bootstrap-corrected",
        "Limited-resection transport: ≥1 examined LN",
        "Adequate-node transport: ≥12 examined LN"
    ),
    `n` = c(
        nrow(dev),
        nrow(dev),
        nrow(dev),
        nrow(transport_any),
        nrow(transport_adequate)
    ),
    `Node-positive, n (%)` = c(
        sprintf("%d (%.1f%%)", sum(dev$node_positive), 100 * mean(dev$node_positive)),
        sprintf("%d (%.1f%%)", sum(dev$node_positive), 100 * mean(dev$node_positive)),
        sprintf("%d (%.1f%%)", sum(dev$node_positive), 100 * mean(dev$node_positive)),
        sprintf("%d (%.1f%%)", sum(transport_any$node_positive), 100 * mean(transport_any$node_positive)),
        sprintf("%d (%.1f%%)", NODE_POS_ADEQ, 100 * OBS_MEAN_TRANSPORT_ADEQ)
    ),
    `AUC` = c(
        fmt3(AUC_T_ONLY),
        fmt3(AUC_DEV_APP),
        fmt3(AUC_DEV_BOOT),
        fmt3(AUC_TRANSPORT_ANY),
        fmt3(AUC_TRANSPORT_ADEQ)
    ),
    `Brier score` = c(
        fmt4(BRIER_T_ONLY),
        fmt4(BRIER_DEV_APP),
        "—",
        fmt4(BRIER_TRANSPORT_ANY),
        fmt4(BRIER_TRANSPORT_ADEQ)
    ),
    `Calibration slope` = c(
        "—",
        clean_zero(SLOPE_DEV_APP),
        sprintf(
            "%.3f (95%% CI %.3f–%.3f)",
            SLOPE_DEV_BOOT,
            SLOPE_DEV_BOOT_LO,
            SLOPE_DEV_BOOT_HI
        ),
        clean_zero(SLOPE_TRANSPORT_ANY),
        clean_zero(SLOPE_TRANSPORT_ADEQ)
    ),
    `Calibration intercept` = c(
        "—",
        clean_zero(INTERCEPT_DEV_APP),
        "—",
        clean_zero(INTERCEPT_TRANSPORT_ANY),
        clean_zero(INTERCEPT_TRANSPORT_ADEQ)
    ),
    `Predicted / observed risk` = c(
        "—",
        fmt_pred_obs(mean(pred_dev, na.rm = TRUE), mean(dev$node_positive)),
        "—",
        fmt_pred_obs(PRED_MEAN_TRANSPORT_ANY, OBS_MEAN_TRANSPORT_ANY),
        fmt_pred_obs(PRED_MEAN_TRANSPORT_ADEQ, OBS_MEAN_TRANSPORT_ADEQ)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

## =============================================================================
## 14. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S2_Model_Performance_FINAL.csv"
)

write.csv(
    table_s2,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 15. Build journal-style three-line table
## =============================================================================

ft <- flextable(table_s2)

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
    j = 1,
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 2:ncol(table_s2),
    align = "center",
    part = "all"
)

ft <- valign(
    ft,
    valign = "center",
    part = "all"
)

ft <- border_remove(ft)

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
    padding.top = 1.0,
    padding.bottom = 1.0,
    padding.left = 1.3,
    padding.right = 1.3,
    part = "all"
)

ft <- width(ft, j = 1, width = 2.85)
ft <- width(ft, j = 2, width = 0.55)
ft <- width(ft, j = 3, width = 1.15)
ft <- width(ft, j = 4, width = 0.60)
ft <- width(ft, j = 5, width = 0.80)
ft <- width(ft, j = 6, width = 1.55)
ft <- width(ft, j = 7, width = 1.05)
ft <- width(ft, j = 8, width = 1.30)

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
    "Supplementary_Table_S2_Model_Performance_FINAL.docx"
)

doc <- read_docx()


doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S2. Performance and transport validation of the final nodal-metastasis risk model",
            fp_text(
                font.family = "Times New Roman",
                font.size = 9.3,
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
    "AUC indicates area under the receiver operating characteristic curve; LN, lymph node. ",
    "The development cohort comprised patients undergoing oncologic colectomy with at least 12 examined lymph nodes ",
    "and known nodal status. The T-stage-only model was fitted in the same development cohort. ",
    "Bootstrap-corrected development performance was taken from the prespecified 500-replicate bootstrap; ",
    "the reported calibration-slope interval is the prespecified bootstrap interval. ",
    "Transport analyses were performed among patients undergoing limited resection with known nodal status and ",
    "at least 1 or at least 12 examined lymph nodes, respectively. ",
    "Predicted / observed risk reports the mean predicted nodal risk followed by the observed node-positive proportion. ",
    "Apparent development calibration slope and intercept are expected to be approximately 1 and 0 because they are ",
    "evaluated in the model-development sample; bootstrap-corrected estimates should be used to assess optimism."
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
## 17. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S2 COMPLETE\n")
cat("============================================================\n")

cat(
    sprintf(
        "Frozen risk-score max error = %.12g\n",
        risk_fp
    )
)

cat("\nDevelopment:\n")
cat(
    sprintf(
        "n=%d | N+=%d | AUC apparent=%.3f | Brier=%.4f | slope apparent=%.3f | intercept=%.3f\n",
        nrow(
            dev
        ),
        sum(
            dev$node_positive
        ),
        AUC_DEV_APP,
        BRIER_DEV_APP,
        SLOPE_DEV_APP,
        INTERCEPT_DEV_APP
    )
)

cat(
    sprintf(
        "T-only AUC=%.3f | Brier=%.4f\n",
        AUC_T_ONLY,
        BRIER_T_ONLY
    )
)

cat(
    sprintf(
        "Bootstrap-corrected AUC=%.3f | slope=%.3f (95%% CI %.3f to %.3f)\n",
        AUC_DEV_BOOT,
        SLOPE_DEV_BOOT,
        SLOPE_DEV_BOOT_LO,
        SLOPE_DEV_BOOT_HI
    )
)


cat("\nTransport >=1 LN:\n")
cat(
    sprintf(
        "n=%d | N+=%d | AUC=%.3f | Brier=%.4f | slope=%.3f | intercept=%.3f | pred=%.3f | obs=%.3f\n",
        nrow(
            transport_any
        ),
        sum(
            transport_any$node_positive
        ),
        AUC_TRANSPORT_ANY,
        BRIER_TRANSPORT_ANY,
        SLOPE_TRANSPORT_ANY,
        INTERCEPT_TRANSPORT_ANY,
        PRED_MEAN_TRANSPORT_ANY,
        OBS_MEAN_TRANSPORT_ANY
    )
)


cat("\nTransport >=12 LN:\n")
cat(
    sprintf(
        "n=%d | N+=%d | AUC=%.3f | Brier=%.4f | slope=%.3f | intercept=%.3f | pred=%.3f | obs=%.3f\n",
        nrow(
            transport_adequate
        ),
        NODE_POS_ADEQ,
        AUC_TRANSPORT_ADEQ,
        BRIER_TRANSPORT_ADEQ,
        SLOPE_TRANSPORT_ADEQ,
        INTERCEPT_TRANSPORT_ADEQ,
        PRED_MEAN_TRANSPORT_ADEQ,
        OBS_MEAN_TRANSPORT_ADEQ
    )
)

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Standard three-line table; no vertical rules.\n")
cat("No frozen model component was changed.\n")
cat("============================================================\n")
