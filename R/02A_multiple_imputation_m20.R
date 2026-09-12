## =============================================================================
## A6 RED-TEAM REVISION — STEP 2A v2 (FIXED)
## Multiple imputation gatekeeper (m = 20)
##
## PURPOSE
##   Prespecified missing-data sensitivity analysis before expensive MI uncertainty
##   estimation.
##
## CORRECTED ANALYTIC FRAMEWORK
##   - Summary Stage is NOT included in the propensity-score model.
##   - Propensity scores are refitted WITHIN each predicted-risk stratum.
##   - Bottom predicted-risk decile remains the PRIMARY low-risk anchor.
##   - risk <5% remains SUPPORTIVE.
##   - Grade and tumor size are multiply imputed (m=20).
##   - Missing-indicator terms are NOT used in the MI analysis.
##
## THIS STEP IS A GATEKEEPER
##   It estimates the m=20 pooled POINT ESTIMATES, between-imputation variation,
##   model performance, and balance. It does NOT claim MI-pooled 95% confidence
##   intervals. If the point estimates remain acceptable, STEP 2B will obtain
##   valid uncertainty using bootstrap-within-imputation + Rubin pooling.
##
## OUTPUTS
##   03_results/REDTEAM_02A_MI_m20_PerImputation_RD.csv
##   03_results/REDTEAM_02A_MI_m20_Pooled_PointEstimates.csv
##   03_results/REDTEAM_02A_MI_m20_ModelPerformance.csv
##   03_results/REDTEAM_02A_MI_m20_Balance.csv
##   03_results/REDTEAM_02A_MI_m20_Report.txt
##
## IMPORTANT
##   Does NOT overwrite any FINAL / Step 1A / Step 1B files.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
RESULT_DIR <- file.path(ROOT, "03_results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

pick_existing <- function(paths) {
    z <- paths[file.exists(paths)]
    if (length(z) == 0) return(NA_character_)
    z[1]
}

PIPE_FILE <- pick_existing(c(
    file.path(DATA_DIR, "Aline_frozen_pipeline_v2.RData"),
    file.path(ROOT, "Aline_frozen_pipeline_v2.RData")
))

RAW_FILE <- pick_existing(c(
    file.path(DATA_DIR, "analytic_cohort_A.csv"),
    file.path(ROOT, "analytic_cohort_A.csv")
))

STEP1B_FILE <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Primary_B1000.csv"
)

if (is.na(PIPE_FILE)) stop("Missing Aline_frozen_pipeline_v2.RData")
if (is.na(RAW_FILE)) stop("Missing analytic_cohort_A.csv")

if (!requireNamespace("mice", quietly = TRUE)) {
    install.packages("mice", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required.")
}

load(PIPE_FILE)
RAW <- read.csv(RAW_FILE, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(RAW) != 6951) stop("Frozen cohort fingerprint failed: N must be 6951.")
if (!exists("prep", mode = "function")) stop("prep() missing from frozen pipeline.")
if (!exists("pipeline", mode = "function")) stop("pipeline() missing from frozen pipeline.")

required_raw <- c(
    "year", "age", "female", "race_eth", "income", "rural_urban",
    "marital", "hist_group", "grade", "T", "size_mm", "exposure",
    "ln_examined", "n_positive", "adequate_ln", "os_months",
    "os_event", "css_event", "other_death"
)

miss_raw <- setdiff(required_raw, names(RAW))
if (length(miss_raw) > 0) {
    stop("RAW is missing required fields: ", paste(miss_raw, collapse = ", "))
}

## =============================================================================
## 1. Frozen derived covariates that do NOT depend on grade/size imputation
## =============================================================================

BASE_WT <- pipeline(RAW)

if (nrow(BASE_WT) != nrow(RAW)) stop("pipeline() changed row count unexpectedly.")

required_base <- c(
    "treat", "race4", "income_hi", "metro", "married",
    "size_i", "logsize"
)

miss_base <- setdiff(required_base, names(BASE_WT))
if (length(miss_base) > 0) {
    stop("Frozen pipeline output missing: ", paste(miss_base, collapse = ", "))
}

## =============================================================================
## 2. Infer the EXACT frozen log-size transformation from observed cases
## =============================================================================

obs_size <- which(
    is.finite(RAW$size_mm) &
    RAW$size_mm > 0 &
    is.finite(BASE_WT$logsize)
)

if (length(obs_size) < 100) stop("Too few observed tumor sizes to infer log transformation.")

candidate_transforms <- list(
    log = function(x) log(x),
    log1p = function(x) log1p(x),
    log10 = function(x) log10(x)
)

transform_error <- vapply(
    candidate_transforms,
    function(fun) {
        max(
            abs(
                fun(RAW$size_mm[obs_size]) -
                BASE_WT$logsize[obs_size]
            ),
            na.rm = TRUE
        )
    },
    numeric(1)
)

best_transform_name <- names(which.min(transform_error))
best_transform_error <- min(transform_error)

if (!is.finite(best_transform_error) || best_transform_error > 1e-6) {
    stop(
        "Could not reproduce frozen logsize transformation. Errors: ",
        paste(
            names(transform_error),
            signif(transform_error, 5),
            sep = "=",
            collapse = " | "
        )
    )
}

log_transform <- candidate_transforms[[best_transform_name]]

## =============================================================================
## 3. Build the m=20 imputation dataset
##
## Summary Stage, lymph-node ascertainment, chemotherapy, and adequate_ln
## are intentionally NOT used as imputation predictors.
##
## Treatment and survival outcomes are included as auxiliary variables, which is
## standard for MI of baseline covariates used in outcome/treatment models.
## =============================================================================

grade_chr <- as.character(RAW$grade)
grade_chr[!grade_chr %in% c("1", "2", "3", "4")] <- NA_character_

size_for_mi <- suppressWarnings(as.numeric(RAW$size_mm))
size_for_mi[!is.finite(size_for_mi) | size_for_mi <= 0] <- NA_real_

MI_DATA <- data.frame(
    size_mm = size_for_mi,
    grade = factor(grade_chr, levels = c("1", "2", "3", "4")),
    age = suppressWarnings(as.numeric(RAW$age)),
    female = factor(RAW$female, levels = c(0, 1)),
    year = suppressWarnings(as.numeric(RAW$year)),
    T = factor(RAW$T, levels = c("T1", "T2", "T3", "T4")),
    hist_group = factor(
        RAW$hist_group,
        levels = c("GCA", "SRCC", "mucinous", "nonmucinous")
    ),
    treat = factor(BASE_WT$treat, levels = c(0, 1)),
    race4 = factor(BASE_WT$race4),
    income_hi = factor(BASE_WT$income_hi, levels = c(0, 1)),
    metro = factor(BASE_WT$metro, levels = c(0, 1)),
    married = factor(BASE_WT$married, levels = c(0, 1)),
    os_months = suppressWarnings(as.numeric(RAW$os_months)),
    os_event = suppressWarnings(as.numeric(RAW$os_event)),
    stringsAsFactors = FALSE
)

## FIX:
## css_event / other_death contain missing values in this SEER cohort.
## In v1 they were used as predictors but were not imputed themselves, which can
## leave PMM recipients without a usable predictive mean. We therefore use only
## fully observed auxiliary predictors.
complete_aux <- c(
    "age", "female", "year", "T", "hist_group", "treat",
    "race4", "income_hi", "metro", "married", "os_months", "os_event"
)

aux_missing <- vapply(
    MI_DATA[complete_aux],
    function(x) sum(is.na(x)),
    integer(1)
)

if (any(aux_missing > 0)) {
    bad_aux <- names(aux_missing)[aux_missing > 0]
    stop(
        paste0(
            "An intended MI predictor is incomplete: ",
            paste(
                bad_aux,
                aux_missing[bad_aux],
                sep = "=",
                collapse = " | "
            ),
            ". Send this error message to ChatGPT."
        )
    )
}

missing_before <- c(
    grade = sum(is.na(MI_DATA$grade)),
    size_mm = sum(is.na(MI_DATA$size_mm))
)

ini <- mice::mice(
    MI_DATA,
    m = 1,
    maxit = 0,
    printFlag = FALSE
)

meth <- ini$method
pred <- ini$predictorMatrix

meth[] <- ""
meth["size_mm"] <- "pmm"
meth["grade"] <- "polyreg"

pred[,] <- 0

aux_vars <- complete_aux

pred["size_mm", c("grade", aux_vars)] <- 1
pred["grade", c("size_mm", aux_vars)] <- 1

pred["size_mm", "size_mm"] <- 0
pred["grade", "grade"] <- 0

set.seed(20260831)

IMP <- mice::mice(
    MI_DATA,
    m = 20,
    maxit = 20,
    method = meth,
    predictorMatrix = pred,
    printFlag = TRUE,
    seed = 20260831
)

## =============================================================================
## 3B. Verify that m=20 actually resolved grade and tumor-size missingness
## =============================================================================

completion_audit <- data.frame(
    imputation = 1:20,
    grade_missing = NA_integer_,
    size_missing = NA_integer_,
    size_nonpositive = NA_integer_,
    stringsAsFactors = FALSE
)

for (m_check in 1:20) {
    cc <- mice::complete(IMP, action = m_check)

    cc_size <- suppressWarnings(as.numeric(cc$size_mm))

    completion_audit$grade_missing[m_check] <- sum(is.na(cc$grade))
    completion_audit$size_missing[m_check] <- sum(is.na(cc_size))
    completion_audit$size_nonpositive[m_check] <- sum(
        is.finite(cc_size) &
        cc_size <= 0
    )
}

if (
    any(completion_audit$grade_missing > 0) ||
    any(completion_audit$size_missing > 0) ||
    any(completion_audit$size_nonpositive > 0)
) {
    print(completion_audit)

    if (!is.null(IMP$loggedEvents)) {
        cat("\nMICE loggedEvents:\n")
        print(IMP$loggedEvents)
    }

    stop(
        paste0(
            "MI completion audit failed. ",
            "Send the printed audit/loggedEvents to ChatGPT."
        )
    )
}

cat(
    "MI completion audit PASS: all 20 datasets have complete grade and positive tumor size.\n\n"
)

## =============================================================================
## 4. Helpers: AUC, calibration, weighted Aalen-Johansen
## =============================================================================

auc_rank <- function(y, p) {
    ok <- is.finite(y) & is.finite(p)
    y <- y[ok]
    p <- p[ok]

    n1 <- sum(y == 1)
    n0 <- sum(y == 0)

    if (n1 == 0 || n0 == 0) return(NA_real_)

    r <- rank(p, ties.method = "average")

    (
        sum(r[y == 1]) -
        n1 * (n1 + 1) / 2
    ) / (n1 * n0)
}

cal_slope <- function(y, p) {
    ok <- is.finite(y) & is.finite(p) & p > 0 & p < 1
    y <- y[ok]
    p <- p[ok]

    if (sum(y == 1) < 10 || sum(y == 0) < 10) return(NA_real_)

    lp <- qlogis(p)

    fit <- tryCatch(
        glm(y ~ lp, family = binomial()),
        error = function(e) NULL
    )

    if (is.null(fit)) return(NA_real_)

    as.numeric(coef(fit)["lp"])
}

aj_w <- function(time, ev1, ev2, wt, start = -1, tau = 60) {
    ok <- !is.na(time) & !is.na(ev1) & !is.na(ev2) &
        !is.na(wt) & is.finite(wt)

    time <- time[ok]
    ev1 <- ev1[ok]
    ev2 <- ev2[ok]
    wt <- wt[ok]

    m <- time > start
    time <- time[m]
    ev1 <- ev1[m]
    ev2 <- ev2[m]
    wt <- wt[m]

    if (length(time) == 0 || sum(wt) <= 0) return(NA_real_)

    dc <- rowsum(wt * ev1, time, reorder = FALSE)
    dr <- rowsum(wt * ev2, time, reorder = FALSE)
    dw <- rowsum(wt, time, reorder = FALSE)

    ut <- as.numeric(rownames(dc))
    o <- order(ut)

    ut <- ut[o]
    dc <- dc[o]
    dr <- dr[o]
    dw <- dw[o]

    S <- 1
    cif <- 0
    ar <- sum(wt)

    for (k in seq_along(ut)) {
        if (ut[k] > tau) break

        if (ar > 0) {
            cif <- cif + S * dc[k] / ar
            S <- S * (1 - (dc[k] + dr[k]) / ar)
        }

        ar <- ar - dw[k]
    }

    as.numeric(cif)
}

ess <- function(w) {
    w <- w[is.finite(w) & w > 0]
    if (length(w) == 0) return(NA_real_)
    sum(w)^2 / sum(w^2)
}

weighted_smd <- function(x, tr, w) {
    ok <- is.finite(x) & is.finite(tr) & is.finite(w)
    x <- x[ok]
    tr <- tr[ok]
    w <- w[ok]

    if (sum(tr == 1) == 0 || sum(tr == 0) == 0) return(NA_real_)

    m1 <- weighted.mean(x[tr == 1], w[tr == 1])
    m0 <- weighted.mean(x[tr == 0], w[tr == 0])

    v1 <- weighted.mean((x[tr == 1] - m1)^2, w[tr == 1])
    v0 <- weighted.mean((x[tr == 0] - m0)^2, w[tr == 0])

    den <- sqrt((v1 + v0) / 2)

    if (!is.finite(den) || den <= 0) return(0)

    (m1 - m0) / den
}

## =============================================================================
## 5. Fixed numeric design matrices for the MI risk model and PS model
##
## MI analysis drops:
##   - size missing indicator
##   - unknown-grade category
##   - Summary Stage
## =============================================================================

make_risk_X <- function(d) {
    data.frame(
        age = d$age,
        age2 = d$age^2,
        female = d$female,
        year = d$year,
        logsize = d$logsize,
        T2 = as.numeric(d$T == "T2"),
        T3 = as.numeric(d$T == "T3"),
        T4 = as.numeric(d$T == "T4"),
        grade2 = as.numeric(d$grade == 2),
        grade3 = as.numeric(d$grade == 3),
        grade4 = as.numeric(d$grade == 4),
        hist_SRCC = as.numeric(d$hist_group == "SRCC"),
        hist_mucinous = as.numeric(d$hist_group == "mucinous"),
        hist_nonmucinous = as.numeric(d$hist_group == "nonmucinous"),
        check.names = FALSE
    )
}

make_ps_X <- function(d) {
    rx <- make_risk_X(d)

    cbind(
        rx,
        race_black = as.numeric(d$race4 == "Non-Hispanic Black"),
        race_white = as.numeric(d$race4 == "Non-Hispanic White"),
        race_other = as.numeric(d$race4 == "Other/Unknown"),
        income_hi = d$income_hi,
        metro = d$metro,
        married = d$married
    )
}

logistic_predict <- function(X_train, y_train, X_new) {
    X_train <- as.matrix(X_train)
    X_new <- as.matrix(X_new)

    if (!identical(colnames(X_train), colnames(X_new))) {
        stop("Training/new design columns differ.")
    }

    keep <- vapply(
        seq_len(ncol(X_train)),
        function(j) {
            x <- X_train[, j]
            is.finite(sd(x, na.rm = TRUE)) && sd(x, na.rm = TRUE) > 0
        },
        logical(1)
    )

    Xtr <- cbind("(Intercept)" = 1, X_train[, keep, drop = FALSE])
    Xnw <- cbind("(Intercept)" = 1, X_new[, keep, drop = FALSE])

    fit <- glm.fit(
        x = Xtr,
        y = y_train,
        family = binomial()
    )

    cf <- fit$coefficients
    cf[!is.finite(cf)] <- 0

    lp <- as.numeric(Xnw %*% cf)
    p <- plogis(lp)

    pmin(pmax(p, 1e-6), 1 - 1e-6)
}

## =============================================================================
## 6. Build one complete MI analytic dataset
## =============================================================================

build_complete <- function(imp_number) {
    comp <- mice::complete(IMP, action = imp_number)

    d <- data.frame(
        patient_id = RAW$patient_id,
        year = as.numeric(RAW$year),
        age = as.numeric(RAW$age),
        female = as.numeric(RAW$female),
        hist_group = as.character(RAW$hist_group),
        grade = as.integer(as.character(comp$grade)),
        T = as.character(RAW$T),
        size_mm = as.numeric(comp$size_mm),
        exposure = as.character(RAW$exposure),
        ln_examined = as.numeric(RAW$ln_examined),
        n_positive = as.numeric(RAW$n_positive),
        adequate_ln = as.numeric(RAW$adequate_ln),
        os_months = as.numeric(RAW$os_months),
        css_event = as.numeric(RAW$css_event),
        other_death = as.numeric(RAW$other_death),
        treat = as.numeric(BASE_WT$treat),
        race4 = as.character(BASE_WT$race4),
        income_hi = as.numeric(BASE_WT$income_hi),
        metro = as.numeric(BASE_WT$metro),
        married = as.numeric(BASE_WT$married),
        stringsAsFactors = FALSE
    )

    if (
        any(!is.finite(d$size_mm)) ||
        any(d$size_mm <= 0)
    ) {
        stop(
            paste0(
                "Completed MI dataset ",
                imp_number,
                " still contains invalid tumor size: missing/nonfinite=",
                sum(!is.finite(d$size_mm)),
                ", nonpositive=",
                sum(is.finite(d$size_mm) & d$size_mm <= 0),
                "."
            )
        )
    }

    d$logsize <- log_transform(d$size_mm)

    if (any(!is.finite(d$logsize))) {
        stop(
            paste0(
                "Non-finite logsize after MI in dataset ",
                imp_number,
                "."
            )
        )
    }

    d$node_pos <- ifelse(
        is.na(d$n_positive),
        NA_real_,
        as.numeric(d$n_positive > 0)
    )

    ## Refit nodal-risk model in the original development structure.
    dev <- (
        d$treat == 1 &
        d$adequate_ln == 1 &
        !is.na(d$node_pos)
    )

    if (sum(dev) < 3000) stop("MI development cohort unexpectedly small.")

    Xall <- make_risk_X(d)
    Xdev <- Xall[dev, , drop = FALSE]

    d$nodal_risk <- logistic_predict(
        X_train = Xdev,
        y_train = d$node_pos[dev],
        X_new = Xall
    )

    attr(d, "dev_n") <- sum(dev)
    attr(d, "dev_pos") <- sum(d$node_pos[dev] == 1)

    d
}

## =============================================================================
## 7. Risk strata + subgroup-specific PS/OW
## =============================================================================

make_masks <- function(d) {
    cuts <- as.numeric(
        quantile(
            d$nodal_risk,
            probs = 0:5 / 5,
            na.rm = TRUE,
            names = FALSE
        )
    )

    cuts[1] <- cuts[1] - 1e-8
    cuts[6] <- cuts[6] + 1e-8

    p10 <- as.numeric(
        quantile(
            d$nodal_risk,
            0.10,
            na.rm = TRUE
        )
    )

    grp <- cut(
        d$nodal_risk,
        breaks = cuts,
        labels = FALSE,
        include.lowest = TRUE
    )

    list(
        Q1 = grp == 1 & !is.na(grp),
        Q2 = grp == 2 & !is.na(grp),
        Q3 = grp == 3 & !is.na(grp),
        Q4 = grp == 4 & !is.na(grp),
        Q5 = grp == 5 & !is.na(grp),
        bottom_decile = d$nodal_risk <= p10,
        risk_lt_5pct = d$nodal_risk < 0.05
    )
}

fit_subgroup_ow <- function(g) {
    X <- make_ps_X(g)

    ps <- logistic_predict(
        X_train = X,
        y_train = g$treat,
        X_new = X
    )

    g$ps_sub <- ps
    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    g
}

estimate_stratum <- function(g) {
    g <- fit_subgroup_ow(g)

    g$cancer <- as.integer(
        g$css_event == 1 &
        g$os_months <= 60
    )

    g$otherd <- as.integer(
        g$other_death == 1 &
        g$os_months <= 60
    )

    a0 <- g[g$treat == 0, , drop = FALSE]
    a1 <- g[g$treat == 1, , drop = FALSE]

    if (nrow(a0) < 10 || nrow(a1) < 10) {
        return(c(
            rd = NA_real_,
            cif0 = NA_real_,
            cif1 = NA_real_,
            max_smd = NA_real_,
            ess0 = NA_real_,
            ess1 = NA_real_
        ))
    }

    cif0 <- aj_w(
        a0$os_months,
        a0$cancer,
        a0$otherd,
        a0$ow_sub
    )

    cif1 <- aj_w(
        a1$os_months,
        a1$cancer,
        a1$otherd,
        a1$ow_sub
    )

    bm <- cbind(
        make_ps_X(g),
        nodal_risk = g$nodal_risk
    )

    sv <- vapply(
        seq_len(ncol(bm)),
        function(j) {
            weighted_smd(
                bm[, j],
                g$treat,
                g$ow_sub
            )
        },
        numeric(1)
    )

    c(
        rd = cif1 - cif0,
        cif0 = cif0,
        cif1 = cif1,
        max_smd = max(abs(sv), na.rm = TRUE),
        ess0 = ess(a0$ow_sub),
        ess1 = ess(a1$ow_sub)
    )
}

## =============================================================================
## 8. Model performance in one imputed dataset
## =============================================================================

model_performance <- function(d) {
    dev <- (
        d$treat == 1 &
        d$adequate_ln == 1 &
        !is.na(d$node_pos)
    )

    tr12 <- (
        d$treat == 0 &
        d$adequate_ln == 1 &
        !is.na(d$node_pos)
    )

    c(
        dev_n = sum(dev),
        dev_pos = sum(d$node_pos[dev] == 1),
        dev_auc = auc_rank(d$node_pos[dev], d$nodal_risk[dev]),
        dev_brier = mean(
            (d$node_pos[dev] - d$nodal_risk[dev])^2,
            na.rm = TRUE
        ),
        transport12_n = sum(tr12),
        transport12_pos = sum(d$node_pos[tr12] == 1),
        transport12_auc = auc_rank(
            d$node_pos[tr12],
            d$nodal_risk[tr12]
        ),
        transport12_brier = mean(
            (d$node_pos[tr12] - d$nodal_risk[tr12])^2,
            na.rm = TRUE
        ),
        transport12_slope = cal_slope(
            d$node_pos[tr12],
            d$nodal_risk[tr12]
        ),
        transport12_pred = mean(
            d$nodal_risk[tr12],
            na.rm = TRUE
        ),
        transport12_obs = mean(
            d$node_pos[tr12],
            na.rm = TRUE
        )
    )
}

## =============================================================================
## 9. Run all 20 imputations
## =============================================================================

labs <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

rd_rows <- list()
perf_rows <- list()
balance_rows <- list()

for (m in 1:20) {
    cat("MI dataset ", m, "/20\n", sep = "")

    d <- build_complete(m)
    masks <- make_masks(d)

    perf <- model_performance(d)

    perf_rows[[m]] <- data.frame(
        imputation = m,
        t(perf),
        check.names = FALSE
    )

    for (nm in labs) {
        g <- d[masks[[nm]], , drop = FALSE]

        est <- tryCatch(
            estimate_stratum(g),
            error = function(e) {
                c(
                    rd = NA_real_,
                    cif0 = NA_real_,
                    cif1 = NA_real_,
                    max_smd = NA_real_,
                    ess0 = NA_real_,
                    ess1 = NA_real_
                )
            }
        )

        rd_rows[[length(rd_rows) + 1]] <- data.frame(
            imputation = m,
            analysis = nm,
            n = nrow(g),
            limited = sum(g$treat == 0),
            colectomy = sum(g$treat == 1),
            RD_pp = est["rd"] * 100,
            CSM_limited_pct = est["cif0"] * 100,
            CSM_colectomy_pct = est["cif1"] * 100,
            max_abs_weighted_SMD = est["max_smd"],
            ESS_limited = est["ess0"],
            ESS_colectomy = est["ess1"],
            stringsAsFactors = FALSE
        )

        balance_rows[[length(balance_rows) + 1]] <- data.frame(
            imputation = m,
            analysis = nm,
            max_abs_weighted_SMD = est["max_smd"],
            balance_pass = is.finite(est["max_smd"]) && est["max_smd"] < 0.10,
            stringsAsFactors = FALSE
        )
    }
}

RD_LONG <- do.call(rbind, rd_rows)
PERF <- do.call(rbind, perf_rows)
BALANCE <- do.call(rbind, balance_rows)

## =============================================================================
## 10. Pool POINT ESTIMATES across m=20
##
## The SD/min/max below describe between-imputation variation only.
## They are NOT inferential 95% confidence intervals.
## =============================================================================

pool_rows <- list()

for (nm in labs) {
    z <- RD_LONG[
        RD_LONG$analysis == nm &
        is.finite(RD_LONG$RD_pp),
        ,
        drop = FALSE
    ]

    pool_rows[[nm]] <- data.frame(
        analysis = nm,
        m_effective = nrow(z),
        pooled_RD_pp = mean(z$RD_pp),
        between_imp_SD_pp = sd(z$RD_pp),
        min_RD_pp = min(z$RD_pp),
        max_RD_pp = max(z$RD_pp),
        mean_CSM_limited_pct = mean(z$CSM_limited_pct),
        mean_CSM_colectomy_pct = mean(z$CSM_colectomy_pct),
        mean_max_abs_SMD = mean(z$max_abs_weighted_SMD),
        worst_max_abs_SMD = max(z$max_abs_weighted_SMD),
        all_balance_pass = all(z$max_abs_weighted_SMD < 0.10),
        mean_ESS_limited = mean(z$ESS_limited),
        mean_ESS_colectomy = mean(z$ESS_colectomy),
        stringsAsFactors = FALSE
    )
}

POOLED <- do.call(rbind, pool_rows)

## =============================================================================
## 11. Compare with corrected Step 1B complete-case/missing-indicator analysis
## =============================================================================

if (file.exists(STEP1B_FILE)) {
    S1B <- read.csv(
        STEP1B_FILE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    rd_col <- if ("RD_pp" %in% names(S1B)) {
        "RD_pp"
    } else if ("rd_pp" %in% names(S1B)) {
        "rd_pp"
    } else {
        NA_character_
    }

    if (!is.na(rd_col)) {
        REF <- data.frame(
            analysis = S1B$analysis,
            step1B_RD_pp = S1B[[rd_col]],
            stringsAsFactors = FALSE
        )

        POOLED <- merge(
            POOLED,
            REF,
            by = "analysis",
            all.x = TRUE,
            sort = FALSE
        )

        POOLED <- POOLED[
            match(labs, POOLED$analysis),
            ,
            drop = FALSE
        ]

        POOLED$change_vs_step1B_pp <- (
            POOLED$pooled_RD_pp -
            POOLED$step1B_RD_pp
        )
    }
}

## =============================================================================
## 12. Performance summaries
## =============================================================================

PERF_SUMMARY <- data.frame(
    metric = setdiff(names(PERF), "imputation"),
    mean = vapply(
        PERF[setdiff(names(PERF), "imputation")],
        function(x) mean(as.numeric(x), na.rm = TRUE),
        numeric(1)
    ),
    sd = vapply(
        PERF[setdiff(names(PERF), "imputation")],
        function(x) sd(as.numeric(x), na.rm = TRUE),
        numeric(1)
    ),
    min = vapply(
        PERF[setdiff(names(PERF), "imputation")],
        function(x) min(as.numeric(x), na.rm = TRUE),
        numeric(1)
    ),
    max = vapply(
        PERF[setdiff(names(PERF), "imputation")],
        function(x) max(as.numeric(x), na.rm = TRUE),
        numeric(1)
    ),
    stringsAsFactors = FALSE
)

## =============================================================================
## 13. Gatekeeping / score-alert logic
## =============================================================================

get_pool <- function(nm, field) {
    POOLED[POOLED$analysis == nm, field][1]
}

bottom_rd <- get_pool("bottom_decile", "pooled_RD_pp")
risk5_rd <- get_pool("risk_lt_5pct", "pooled_RD_pp")
q5_rd <- get_pool("Q5", "pooled_RD_pp")

bottom_shift <- if ("change_vs_step1B_pp" %in% names(POOLED)) {
    get_pool("bottom_decile", "change_vs_step1B_pp")
} else {
    NA_real_
}

risk5_shift <- if ("change_vs_step1B_pp" %in% names(POOLED)) {
    get_pool("risk_lt_5pct", "change_vs_step1B_pp")
} else {
    NA_real_
}

q5_shift <- if ("change_vs_step1B_pp" %in% names(POOLED)) {
    get_pool("Q5", "change_vs_step1B_pp")
} else {
    NA_real_
}

score_alert <- "NO IMMEDIATE POINT-ESTIMATE DOWNGRADE"

if (
    is.finite(bottom_shift) &&
    abs(bottom_shift) > 1.0
) {
    score_alert <- paste0(
        "CAUTION: bottom-decile pooled MI point estimate moved ",
        sprintf("%+.2f", bottom_shift),
        " pp vs Step 1B. Proceed to Step 2B before preserving the current result-convincingness score."
    )
}

if (
    is.finite(q5_shift) &&
    abs(q5_shift) > 2.0
) {
    score_alert <- paste0(
        "PROTOCOL TRIGGER: Q5 moved ",
        sprintf("%+.2f", q5_shift),
        " pp (>2 pp). The prespecified plan requires considering MI as the primary missing-data analysis."
    )
}

if (
    is.finite(bottom_rd) &&
    bottom_rd < -1.0
) {
    score_alert <- paste0(
        "SCORE ALERT: MI moved the primary bottom-decile point estimate toward a colectomy-associated reduction (",
        sprintf("%+.2f", bottom_rd),
        " pp). Result-convincingness should be provisionally downgraded pending Step 2B uncertainty."
    )
}

## =============================================================================
## 14. Save
## =============================================================================

OUT_LONG <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_PerImputation_RD.csv"
)

OUT_POOL <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_Pooled_PointEstimates.csv"
)

OUT_PERF <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_ModelPerformance.csv"
)

OUT_BAL <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_Balance.csv"
)

OUT_REPORT <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_Report.txt"
)

write.csv(RD_LONG, OUT_LONG, row.names = FALSE)
write.csv(POOLED, OUT_POOL, row.names = FALSE)
write.csv(PERF_SUMMARY, OUT_PERF, row.names = FALSE)
write.csv(BALANCE, OUT_BAL, row.names = FALSE)

report <- c(
    "======================================================================",
    "A6 RED-TEAM STEP 2A — MULTIPLE IMPUTATION GATEKEEPER (m=20)",
    "======================================================================",
    paste0("Original grade missing: ", missing_before["grade"], " / ", nrow(RAW)),
    paste0("Original tumor-size missing: ", missing_before["size_mm"], " / ", nrow(RAW)),
    paste0(
        "Frozen logsize transformation reproduced as: ",
        best_transform_name,
        " | max error=",
        format(best_transform_error, scientific = TRUE, digits = 4)
    ),
    "",
    "MI specification:",
    "  grade: multinomial logistic imputation (polyreg)",
    "  tumor size: predictive mean matching (pmm)",
    "  m=20; maxit=20",
    "  Summary Stage NOT used in imputation or PS.",
    "  Nodal ascertainment / LN count NOT used as imputation predictors.",
    "  Missing-indicator terms NOT used in the MI risk/PS models.",
    "  Treatment, overall-survival time, and overall-death indicator included as complete auxiliary predictors.",
    "  Cause-specific outcome variables were excluded from MI predictors because they contain missing values.",
    "",
    "== m=20 POOLED POINT ESTIMATES ==",
    "NOTE: min/max and between-imputation SD are NOT 95% confidence intervals."
)

for (i in seq_len(nrow(POOLED))) {
    z <- POOLED[i, ]

    shift_text <- if ("change_vs_step1B_pp" %in% names(POOLED)) {
        sprintf(" | change vs Step1B %+.2f pp", z$change_vs_step1B_pp)
    } else {
        ""
    }

    report <- c(
        report,
        paste0(
            sprintf(
                "%-15s pooled RD %+.2f pp | between-MI SD %.2f | range %+.2f to %+.2f | worst max|SMD| %.3f | balance %s",
                z$analysis,
                z$pooled_RD_pp,
                z$between_imp_SD_pp,
                z$min_RD_pp,
                z$max_RD_pp,
                z$worst_max_abs_SMD,
                ifelse(z$all_balance_pass, "PASS", "FAIL")
            ),
            shift_text
        )
    )
}

report <- c(
    report,
    "",
    "== MODEL PERFORMANCE ACROSS IMPUTATIONS =="
)

for (i in seq_len(nrow(PERF_SUMMARY))) {
    z <- PERF_SUMMARY[i, ]

    report <- c(
        report,
        sprintf(
            "%-24s mean %.4f | SD %.4f | range %.4f to %.4f",
            z$metric,
            z$mean,
            z$sd,
            z$min,
            z$max
        )
    )
}

report <- c(
    report,
    "",
    "== SCORE / PROTOCOL ALERT ==",
    score_alert,
    "",
    "Interpretation rule:",
    "STEP 2A is only the point-estimate gatekeeper. The prespecified 3-pp criterion CANNOT be finally judged from this file because valid MI-pooled 95% uncertainty has not yet been calculated.",
    "If Step 2A is acceptable, run Step 2B for bootstrap-within-imputation + Rubin-pooled uncertainty.",
    "======================================================================"
)

writeLines(report, OUT_REPORT, useBytes = TRUE)

cat(paste(report, collapse = "\n"), "\n")
cat("\nREPORT: ", OUT_REPORT, "\n", sep = "")
