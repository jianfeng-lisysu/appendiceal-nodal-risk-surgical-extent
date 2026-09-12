## =============================================================================
## A6 RED-TEAM REVISION — STEP 2B
## Multiple imputation (m=20) + bootstrap-within-imputation + Rubin pooling
##
## PURPOSE
##   Final uncertainty analysis for the prespecified missing-data sensitivity.
##
## CORRECTED FRAMEWORK
##   - Summary Stage excluded from PS.
##   - Grade + tumor size multiply imputed (m=20).
##   - Missing-indicator terms removed in MI risk/PS models.
##   - Nodal-risk model refitted in every bootstrap replicate.
##   - Risk cut points recreated in every bootstrap replicate.
##   - PS refitted WITHIN each risk stratum.
##   - Bottom predicted-risk decile remains PRIMARY low-risk anchor.
##   - risk <5% remains SUPPORTIVE.
##
## UNCERTAINTY
##   For each of 20 imputed datasets:
##     point estimate q_m
##     bootstrap sampling variance U_m (default B=250)
##
##   Rubin pooling:
##     Qbar = mean(q_m)
##     Ubar = mean(U_m)
##     B    = var(q_m)
##     T    = Ubar + (1 + 1/m) * B
##     95% CI uses Rubin degrees of freedom.
##
## CHECKPOINT / RESUME
##   Each imputation's bootstrap results are saved separately.
##   If interrupted, rerun the same script and completed imputations are skipped.
##
## FINAL CONFIRMATION RUN: B=1000 PER IMPUTATION (total 20,000 bootstrap replicates).
##   This is the final Monte Carlo confirmation run for the primary 3-pp margin.
##   Do not escalate B further based on whether the result is favorable or unfavorable.
##
## OUTPUTS
##   03_results/REDTEAM_02B_MI_Rubin_Final.csv
##   03_results/REDTEAM_02B_MI_PerImputation.csv
##   03_results/REDTEAM_02B_MI_BootstrapVariance.csv
##   03_results/REDTEAM_02B_MI_Report.txt
##   03_results/REDTEAM_02B_checkpoints/*.csv
##
## Does NOT overwrite any existing FINAL / Step1 / Step2A files.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
RESULT_DIR <- file.path(ROOT, "03_results")
CHECK_DIR <- file.path(RESULT_DIR, "REDTEAM_02B_checkpoints")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)

M_IMP <- 20L
B_BOOT <- 1000L
MASTER_SEED <- 20260901L

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

STEP2A_FILE <- file.path(
    RESULT_DIR,
    "REDTEAM_02A_MI_m20_Pooled_PointEstimates.csv"
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

if (!exists("pipeline", mode = "function")) stop("pipeline() missing.")
if (!exists("prep", mode = "function")) stop("prep() missing.")

RAW <- read.csv(
    RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(RAW) != 6951) {
    stop("Frozen cohort fingerprint failed: expected N=6951.")
}

BASE_WT <- pipeline(RAW)

required_base <- c(
    "treat", "race4", "income_hi", "metro", "married", "logsize"
)

miss_base <- setdiff(required_base, names(BASE_WT))
if (length(miss_base) > 0) {
    stop("Frozen pipeline output missing: ", paste(miss_base, collapse = ", "))
}

## =============================================================================
## 1. Reproduce frozen log-size transformation
## =============================================================================

obs_size <- which(
    is.finite(RAW$size_mm) &
    RAW$size_mm > 0 &
    is.finite(BASE_WT$logsize)
)

candidate_transforms <- list(
    log = function(x) log(x),
    log1p = function(x) log1p(x),
    log10 = function(x) log10(x)
)

transform_error <- vapply(
    candidate_transforms,
    function(fun) {
        max(
            abs(fun(RAW$size_mm[obs_size]) - BASE_WT$logsize[obs_size]),
            na.rm = TRUE
        )
    },
    numeric(1)
)

best_transform_name <- names(which.min(transform_error))
best_transform_error <- min(transform_error)

if (!is.finite(best_transform_error) || best_transform_error > 1e-6) {
    stop("Could not reproduce frozen logsize transformation.")
}

log_transform <- candidate_transforms[[best_transform_name]]

## =============================================================================
## 2. Recreate the exact Step 2A m=20 imputation
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
    stop(
        "An intended MI auxiliary predictor is incomplete: ",
        paste(
            names(aux_missing)[aux_missing > 0],
            aux_missing[aux_missing > 0],
            sep = "=",
            collapse = " | "
        )
    )
}

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
pred["size_mm", c("grade", complete_aux)] <- 1
pred["grade", c("size_mm", complete_aux)] <- 1
pred["size_mm", "size_mm"] <- 0
pred["grade", "grade"] <- 0

set.seed(20260831)

IMP <- mice::mice(
    MI_DATA,
    m = M_IMP,
    maxit = 20,
    method = meth,
    predictorMatrix = pred,
    printFlag = TRUE,
    seed = 20260831
)

## Completion audit
for (m in seq_len(M_IMP)) {
    cc <- mice::complete(IMP, action = m)
    cc_size <- suppressWarnings(as.numeric(cc$size_mm))

    if (
        any(is.na(cc$grade)) ||
        any(!is.finite(cc_size)) ||
        any(cc_size <= 0)
    ) {
        stop("MI completion audit failed in imputation ", m)
    }
}

cat("MI completion audit PASS for all ", M_IMP, " imputations.\n\n", sep = "")

## =============================================================================
## 3. Core helpers
## =============================================================================

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
            s <- sd(X_train[, j], na.rm = TRUE)
            is.finite(s) && s > 0
        },
        logical(1)
    )

    Xtr <- cbind("(Intercept)" = 1, X_train[, keep, drop = FALSE])
    Xnw <- cbind("(Intercept)" = 1, X_new[, keep, drop = FALSE])

    fit <- suppressWarnings(
        glm.fit(
            x = Xtr,
            y = y_train,
            family = binomial()
        )
    )

    cf <- fit$coefficients
    cf[!is.finite(cf)] <- 0

    p <- plogis(as.numeric(Xnw %*% cf))
    pmin(pmax(p, 1e-6), 1 - 1e-6)
}

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

    if (any(!is.finite(d$size_mm)) || any(d$size_mm <= 0)) {
        stop("Invalid tumor size after MI in dataset ", imp_number)
    }

    d$logsize <- log_transform(d$size_mm)

    d$node_pos <- ifelse(
        is.na(d$n_positive),
        NA_real_,
        as.numeric(d$n_positive > 0)
    )

    d
}

refit_nodal_risk <- function(d) {
    dev <- (
        d$treat == 1 &
        d$adequate_ln == 1 &
        !is.na(d$node_pos)
    )

    if (sum(dev) < 500 || sum(d$node_pos[dev] == 1) < 50) {
        stop("Bootstrap development cohort insufficient.")
    }

    Xall <- make_risk_X(d)
    Xdev <- Xall[dev, , drop = FALSE]

    d$nodal_risk <- logistic_predict(
        X_train = Xdev,
        y_train = d$node_pos[dev],
        X_new = Xall
    )

    d
}

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
    if (sum(g$treat == 0) < 10 || sum(g$treat == 1) < 10) {
        stop("Treatment group too small in risk stratum.")
    }

    X <- make_ps_X(g)

    ps <- logistic_predict(
        X_train = X,
        y_train = g$treat,
        X_new = X
    )

    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    g
}

estimate_one_stratum <- function(g) {
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

    c0 <- aj_w(
        a0$os_months,
        a0$cancer,
        a0$otherd,
        a0$ow_sub
    )

    c1 <- aj_w(
        a1$os_months,
        a1$cancer,
        a1$otherd,
        a1$ow_sub
    )

    c1 - c0
}

labs <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

estimate_all <- function(d, refit_risk = TRUE) {
    if (refit_risk) {
        d <- refit_nodal_risk(d)
    }

    masks <- make_masks(d)

    out <- rep(NA_real_, length(labs))
    names(out) <- labs

    for (i in seq_along(labs)) {
        g <- d[masks[[labs[i]]], , drop = FALSE]

        out[i] <- tryCatch(
            estimate_one_stratum(g),
            error = function(e) NA_real_
        )
    }

    out
}

## =============================================================================
## 4. Per-imputation point estimates + bootstrap variances
## =============================================================================

QMAT <- matrix(
    NA_real_,
    nrow = M_IMP,
    ncol = length(labs),
    dimnames = list(
        paste0("MI", seq_len(M_IMP)),
        labs
    )
)

UMAT <- matrix(
    NA_real_,
    nrow = M_IMP,
    ncol = length(labs),
    dimnames = list(
        paste0("MI", seq_len(M_IMP)),
        labs
    )
)

EFFMAT <- matrix(
    0L,
    nrow = M_IMP,
    ncol = length(labs),
    dimnames = list(
        paste0("MI", seq_len(M_IMP)),
        labs
    )
)

set.seed(MASTER_SEED)

for (m in seq_len(M_IMP)) {
    cat("\n============================================================\n")
    cat("STEP 2B | imputation ", m, "/", M_IMP, "\n", sep = "")
    cat("============================================================\n")

    d <- build_complete(m)
    d_point <- refit_nodal_risk(d)

    QMAT[m, ] <- estimate_all(
        d_point,
        refit_risk = FALSE
    )

    check_file <- file.path(
        CHECK_DIR,
        sprintf(
            "MI%02d_BOOT_B%d.csv",
            m,
            B_BOOT
        )
    )

    boot_m <- NULL

    if (file.exists(check_file)) {
        old <- tryCatch(
            read.csv(
                check_file,
                stringsAsFactors = FALSE,
                check.names = FALSE
            ),
            error = function(e) NULL
        )

        if (
            !is.null(old) &&
            nrow(old) == B_BOOT &&
            all(labs %in% names(old))
        ) {
            boot_m <- as.matrix(old[, labs, drop = FALSE])
            storage.mode(boot_m) <- "numeric"
            cat("Checkpoint found — bootstrap skipped.\n")
        }
    }

    if (is.null(boot_m)) {
        boot_m <- matrix(
            NA_real_,
            nrow = B_BOOT,
            ncol = length(labs),
            dimnames = list(NULL, labs)
        )

        set.seed(MASTER_SEED + m)

        t0 <- Sys.time()

        for (b in seq_len(B_BOOT)) {
            idx <- sample(
                seq_len(nrow(d)),
                size = nrow(d),
                replace = TRUE
            )

            db <- d[idx, , drop = FALSE]

            boot_m[b, ] <- tryCatch(
                estimate_all(
                    db,
                    refit_risk = TRUE
                ),
                error = function(e) {
                    rep(NA_real_, length(labs))
                }
            )

            if (b %% 50 == 0) {
                cat(
                    "  bootstrap ",
                    b,
                    "/",
                    B_BOOT,
                    " | elapsed ",
                    round(
                        difftime(
                            Sys.time(),
                            t0,
                            units = "mins"
                        ),
                        1
                    ),
                    " min\n",
                    sep = ""
                )
            }
        }

        write.csv(
            data.frame(
                replicate = seq_len(B_BOOT),
                boot_m,
                check.names = FALSE
            ),
            check_file,
            row.names = FALSE
        )
    }

    for (j in seq_along(labs)) {
        v <- boot_m[, j]
        v <- v[is.finite(v)]

        EFFMAT[m, j] <- length(v)

        if (length(v) >= ceiling(0.90 * B_BOOT)) {
            UMAT[m, j] <- var(v)
        }
    }

    if (any(EFFMAT[m, ] < ceiling(0.90 * B_BOOT))) {
        warning(
            "Imputation ",
            m,
            " has <90% effective bootstrap replicates for: ",
            paste(
                labs[EFFMAT[m, ] < ceiling(0.90 * B_BOOT)],
                collapse = ", "
            )
        )
    }
}

## =============================================================================
## 5. Rubin pooling
## =============================================================================

pool_one <- function(q, u, effective_B, label) {
    ok <- is.finite(q) & is.finite(u)

    q <- q[ok]
    u <- u[ok]
    effective_B <- effective_B[ok]

    m_eff <- length(q)

    if (m_eff < 10) {
        return(
            data.frame(
                analysis = label,
                m_effective = m_eff,
                pooled_RD = NA_real_,
                pooled_RD_pp = NA_real_,
                Ubar = NA_real_,
                B_between = NA_real_,
                total_variance = NA_real_,
                SE_pp = NA_real_,
                df = NA_real_,
                lo_pp = NA_real_,
                hi_pp = NA_real_,
                upper_possible_reduction_pp = NA_real_,
                criterion_2pp = NA,
                criterion_3pp = NA,
                criterion_5pp = NA,
                mean_bootstrap_B = mean(effective_B),
                min_bootstrap_B = min(effective_B),
                stringsAsFactors = FALSE
            )
        )
    }

    qbar <- mean(q)
    ubar <- mean(u)
    bvar <- var(q)

    total <- ubar + (1 + 1 / m_eff) * bvar

    if (!is.finite(total) || total <= 0) {
        stop("Invalid total Rubin variance for ", label)
    }

    se <- sqrt(total)

    extra <- (1 + 1 / m_eff) * bvar

    if (!is.finite(extra) || extra <= 0) {
        df <- Inf
    } else {
        df <- (m_eff - 1) * (1 + ubar / extra)^2
    }

    crit <- if (is.finite(df)) {
        qt(0.975, df = df)
    } else {
        qnorm(0.975)
    }

    lo <- qbar - crit * se
    hi <- qbar + crit * se

    upper_reduction_pp <- max(0, -lo * 100)

    data.frame(
        analysis = label,
        m_effective = m_eff,
        pooled_RD = qbar,
        pooled_RD_pp = qbar * 100,
        Ubar = ubar,
        B_between = bvar,
        total_variance = total,
        SE_pp = se * 100,
        df = df,
        lo_pp = lo * 100,
        hi_pp = hi * 100,
        upper_possible_reduction_pp = upper_reduction_pp,
        criterion_2pp = upper_reduction_pp < 2,
        criterion_3pp = upper_reduction_pp < 3,
        criterion_5pp = upper_reduction_pp < 5,
        mean_bootstrap_B = mean(effective_B),
        min_bootstrap_B = min(effective_B),
        stringsAsFactors = FALSE
    )
}

FINAL_LIST <- lapply(
    seq_along(labs),
    function(j) {
        pool_one(
            q = QMAT[, j],
            u = UMAT[, j],
            effective_B = EFFMAT[, j],
            label = labs[j]
        )
    }
)

FINAL <- do.call(rbind, FINAL_LIST)

## =============================================================================
## 6. Comparisons with Step 1B and Step 2A
## =============================================================================

if (file.exists(STEP1B_FILE)) {
    s1 <- read.csv(
        STEP1B_FILE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    rd_col <- if ("RD_pp" %in% names(s1)) {
        "RD_pp"
    } else if ("rd_pp" %in% names(s1)) {
        "rd_pp"
    } else {
        NA_character_
    }

    if (!is.na(rd_col)) {
        ref1 <- data.frame(
            analysis = s1$analysis,
            Step1B_RD_pp = s1[[rd_col]],
            stringsAsFactors = FALSE
        )

        FINAL <- merge(
            FINAL,
            ref1,
            by = "analysis",
            all.x = TRUE,
            sort = FALSE
        )

        FINAL <- FINAL[
            match(labs, FINAL$analysis),
            ,
            drop = FALSE
        ]

        FINAL$change_vs_Step1B_pp <- (
            FINAL$pooled_RD_pp -
            FINAL$Step1B_RD_pp
        )
    }
}

if (file.exists(STEP2A_FILE)) {
    s2 <- read.csv(
        STEP2A_FILE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if ("pooled_RD_pp" %in% names(s2)) {
        ref2 <- data.frame(
            analysis = s2$analysis,
            Step2A_pooled_RD_pp = s2$pooled_RD_pp,
            stringsAsFactors = FALSE
        )

        FINAL <- merge(
            FINAL,
            ref2,
            by = "analysis",
            all.x = TRUE,
            sort = FALSE
        )

        FINAL <- FINAL[
            match(labs, FINAL$analysis),
            ,
            drop = FALSE
        ]

        FINAL$point_check_vs_Step2A_pp <- (
            FINAL$pooled_RD_pp -
            FINAL$Step2A_pooled_RD_pp
        )
    }
}

## =============================================================================
## 7. Per-imputation audit tables
## =============================================================================

PER_IMP <- data.frame(
    imputation = seq_len(M_IMP),
    QMAT,
    check.names = FALSE
)

for (nm in labs) {
    PER_IMP[[paste0(nm, "_pp")]] <- PER_IMP[[nm]] * 100
}

VAR_LONG <- do.call(
    rbind,
    lapply(
        seq_len(M_IMP),
        function(m) {
            data.frame(
                imputation = m,
                analysis = labs,
                point_RD_pp = QMAT[m, ] * 100,
                bootstrap_variance = UMAT[m, ],
                bootstrap_SE_pp = sqrt(UMAT[m, ]) * 100,
                effective_bootstrap_B = EFFMAT[m, ],
                stringsAsFactors = FALSE
            )
        }
    )
)

## =============================================================================
## 8. Score / decision alert
## =============================================================================

get_final <- function(nm, field) {
    FINAL[FINAL$analysis == nm, field][1]
}

bottom_rd <- get_final("bottom_decile", "pooled_RD_pp")
bottom_lo <- get_final("bottom_decile", "lo_pp")
bottom_hi <- get_final("bottom_decile", "hi_pp")
bottom_upper <- get_final("bottom_decile", "upper_possible_reduction_pp")
bottom_met3 <- get_final("bottom_decile", "criterion_3pp")

risk5_met3 <- get_final("risk_lt_5pct", "criterion_3pp")
risk5_upper <- get_final("risk_lt_5pct", "upper_possible_reduction_pp")

q5_shift <- if ("change_vs_Step1B_pp" %in% names(FINAL)) {
    get_final("Q5", "change_vs_Step1B_pp")
} else {
    NA_real_
}

if (isTRUE(bottom_met3)) {
    score_alert <- paste0(
        "PRIMARY RESULT SURVIVES MI: bottom decile pooled RD ",
        sprintf("%+.2f", bottom_rd),
        " pp (95% CI ",
        sprintf("%+.2f", bottom_lo),
        " to ",
        sprintf("%+.2f", bottom_hi),
        "); upper compatible colectomy-associated reduction ",
        sprintf("%.2f", bottom_upper),
        " pp < 3 pp."
    )
} else {
    score_alert <- paste0(
        "SCORE DOWNGRADE REQUIRED: bottom-decile MI-pooled 95% CI no longer excludes a 3-pp colectomy-associated reduction. ",
        "Pooled RD ",
        sprintf("%+.2f", bottom_rd),
        " pp (95% CI ",
        sprintf("%+.2f", bottom_lo),
        " to ",
        sprintf("%+.2f", bottom_hi),
        "); upper compatible reduction ",
        sprintf("%.2f", bottom_upper),
        " pp."
    )
}

protocol_alert <- if (
    is.finite(q5_shift) &&
    abs(q5_shift) > 2
) {
    paste0(
        "PROTOCOL TRIGGER: Q5 shifted ",
        sprintf("%+.2f", q5_shift),
        " pp vs Step 1B (>2 pp). MI should be promoted to the primary missing-data analysis."
    )
} else {
    "No >2-pp Q5 protocol trigger."
}

mc_alert <- paste0(
    "FINAL MONTE CARLO CONFIRMATION COMPLETED at B=1000 per imputation. ",
    "Do not rerun with a larger B based on whether the 3-pp criterion is met. ",
    "Lock the conclusion from this run."
)

## =============================================================================
## 9. Save final files
## =============================================================================

OUT_FINAL <- file.path(
    RESULT_DIR,
    "REDTEAM_02B_MI_Rubin_Final.csv"
)

OUT_IMP <- file.path(
    RESULT_DIR,
    "REDTEAM_02B_MI_PerImputation.csv"
)

OUT_VAR <- file.path(
    RESULT_DIR,
    "REDTEAM_02B_MI_BootstrapVariance.csv"
)

OUT_REPORT <- file.path(
    RESULT_DIR,
    "REDTEAM_02B_MI_Report.txt"
)

write.csv(FINAL, OUT_FINAL, row.names = FALSE)
write.csv(PER_IMP, OUT_IMP, row.names = FALSE)
write.csv(VAR_LONG, OUT_VAR, row.names = FALSE)

report <- c(
    "======================================================================",
    "A6 RED-TEAM STEP 2B — MI + BOOTSTRAP-WITHIN-IMPUTATION + RUBIN POOLING",
    "======================================================================",
    paste0("m = ", M_IMP),
    paste0("bootstrap per imputation = ", B_BOOT),
    paste0("total planned bootstrap replicates = ", M_IMP * B_BOOT),
    paste0(
        "logsize transformation = ",
        best_transform_name,
        " | reproduction error = ",
        format(best_transform_error, scientific = TRUE, digits = 4)
    ),
    "Summary Stage excluded from PS.",
    "PS refitted separately within each predicted-risk stratum.",
    "Bottom decile remains the primary low-risk anchor; risk<5% remains supportive.",
    "",
    "== RUBIN-POOLED CANCER-SPECIFIC MORTALITY RD =="
)

for (i in seq_len(nrow(FINAL))) {
    z <- FINAL[i, ]

    low_text <- ""

    if (z$analysis %in% c("Q1", "bottom_decile", "risk_lt_5pct")) {
        low_text <- sprintf(
            " | upper reduction %.2f pp | 2pp:%s 3pp:%s 5pp:%s",
            z$upper_possible_reduction_pp,
            ifelse(z$criterion_2pp, "MET", "NOT MET"),
            ifelse(z$criterion_3pp, "MET", "NOT MET"),
            ifelse(z$criterion_5pp, "MET", "NOT MET")
        )
    }

    change_text <- if ("change_vs_Step1B_pp" %in% names(FINAL)) {
        sprintf(
            " | change vs Step1B %+.2f pp",
            z$change_vs_Step1B_pp
        )
    } else {
        ""
    }

    report <- c(
        report,
        paste0(
            sprintf(
                "%-15s RD %+.2f pp (95%% CI %+.2f to %+.2f) | SE %.2f pp | df %.1f | m=%d | mean/min boot B %.1f/%d",
                z$analysis,
                z$pooled_RD_pp,
                z$lo_pp,
                z$hi_pp,
                z$SE_pp,
                z$df,
                z$m_effective,
                z$mean_bootstrap_B,
                z$min_bootstrap_B
            ),
            low_text,
            change_text
        )
    )
}

report <- c(
    report,
    "",
    "== PRIMARY DECISION ==",
    score_alert,
    "",
    "== SUPPORTIVE <5% THRESHOLD ==",
    paste0(
        "risk<5% 3-pp criterion: ",
        ifelse(isTRUE(risk5_met3), "MET", "NOT MET"),
        " | upper compatible reduction ",
        sprintf("%.2f", risk5_upper),
        " pp."
    ),
    "",
    "== PROTOCOL CHECK ==",
    protocol_alert,
    "",
    "== MONTE CARLO CHECK ==",
    mc_alert,
    "",
    "Interpretation:",
    "RD = oncologic colectomy minus limited resection in 5-year cancer-specific mortality.",
    "Negative values indicate lower cancer-specific mortality associated with colectomy.",
    "The primary 3-pp bounding claim is met only when the lower 95% confidence limit is greater than -3.0 pp.",
    "======================================================================"
)

writeLines(report, OUT_REPORT, useBytes = TRUE)

cat(paste(report, collapse = "\n"), "\n")
cat("\nREPORT: ", OUT_REPORT, "\n", sep = "")
