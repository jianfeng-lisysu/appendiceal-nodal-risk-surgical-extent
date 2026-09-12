## =============================================================================
## A6 RED-TEAM REVISION — STEP 1B
## Risk-stratum-specific propensity-score refitting after removing Summary Stage
##
## Purpose
##   1. Keep the corrected NO-SUMMARY-STAGE PS specification.
##   2. Refit the PS separately within Q1-Q5, bottom decile, and risk<5%.
##   3. Recalculate subgroup-specific overlap weights.
##   4. Re-estimate weighted 5-year cancer-specific mortality RD.
##   5. Audit within-stratum balance, ESS, PS overlap, and max weight.
##   6. Repeat the full pipeline in B=1000 bootstrap replicates:
##        raw resampling
##        -> nodal-risk model refit
##        -> nodal-risk strata re-created
##        -> PS refit INSIDE each risk stratum
##        -> subgroup-specific OW
##        -> weighted competing-risk RD
##
## IMPORTANT
##   - Does NOT overwrite any existing FINAL or REDTEAM_01A files.
##   - Primary low-risk anchor remains the bottom predicted-risk decile.
##   - risk<5% remains supportive.
##   - Q1-Q5 are retained as secondary risk-stratified analyses.
##
## Output
##   03_results/REDTEAM_01B_SubgroupOW_Primary_B1000.csv
##   03_results/REDTEAM_01B_SubgroupOW_Diagnostics.csv
##   03_results/REDTEAM_01B_SubgroupOW_Balance_Components.csv
##   03_results/REDTEAM_01B_SubgroupOW_Bootstrap_B1000.csv
##   03_results/REDTEAM_01B_SubgroupOW_Report.txt
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

STEP1A_FILE <- file.path(
    RESULT_DIR,
    "REDTEAM_01A_NoSummaryStage_Primary_B1000.csv"
)

if (is.na(PIPE_FILE)) stop("Missing Aline_frozen_pipeline_v2.RData")
if (is.na(RAW_FILE)) stop("Missing analytic_cohort_A.csv")

load(PIPE_FILE)

if (!exists("pipeline", mode = "function")) stop("pipeline() missing")
if (!exists("PS_FORMULA")) stop("PS_FORMULA missing")

RAW <- read.csv(
    RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(RAW) != 6951) stop("Frozen cohort fingerprint failed: N must be 6951")

## =============================================================================
## 1. Correct PS formula: remove Summary Stage
## =============================================================================

PS_ORIG <- PS_FORMULA
terms_orig <- attr(terms(PS_ORIG), "term.labels")
stage_terms <- terms_orig[grepl("stage", terms_orig, ignore.case = TRUE)]

if (length(stage_terms) == 0) {
    stop("No stage-related term found in PS_FORMULA")
}

PS_NOSTAGE <- reformulate(
    termlabels = setdiff(terms_orig, stage_terms),
    response = all.vars(PS_ORIG)[1],
    env = environment(PS_ORIG)
)

pipeline_env <- environment(pipeline)

had_env_formula <- exists(
    "PS_FORMULA",
    envir = pipeline_env,
    inherits = FALSE
)

old_env_formula <- if (had_env_formula) {
    get("PS_FORMULA", envir = pipeline_env, inherits = FALSE)
} else {
    NULL
}

restore_formula <- function() {
    assign("PS_FORMULA", PS_ORIG, envir = .GlobalEnv)
    if (had_env_formula) {
        assign("PS_FORMULA", old_env_formula, envir = pipeline_env)
    }
}
on.exit(restore_formula(), add = TRUE)

assign("PS_FORMULA", PS_NOSTAGE, envir = .GlobalEnv)
assign("PS_FORMULA", PS_NOSTAGE, envir = pipeline_env)

## =============================================================================
## 2. Common estimators
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

make_outcome <- function(d) {
    z <- d[
        !is.na(d$css_event) &
        !is.na(d$os_months) &
        !is.na(d$treat) &
        !is.na(d$nodal_risk),
        ,
        drop = FALSE
    ]

    z$cancer <- as.integer(
        z$css_event == 1 &
        z$os_months <= 60
    )

    z$otherd <- as.integer(
        z$other_death == 1 &
        z$os_months <= 60
    )

    z
}

make_cuts <- function(x) {
    z <- as.numeric(
        quantile(
            x,
            probs = 0:5 / 5,
            na.rm = TRUE,
            names = FALSE
        )
    )

    z[1] <- z[1] - 1e-8
    z[6] <- z[6] + 1e-8
    z
}

make_masks <- function(d) {
    cuts <- make_cuts(d$nodal_risk)

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

## =============================================================================
## 3. Subgroup-specific PS / overlap weights
## =============================================================================

fit_subgroup_ow <- function(g) {
    if (nrow(g) < 40) stop("Subgroup too small")
    if (sum(g$treat == 0) < 15) stop("Too few limited patients")
    if (sum(g$treat == 1) < 15) stop("Too few colectomy patients")

    fit <- suppressWarnings(
        glm(
            PS_NOSTAGE,
            data = g,
            family = binomial()
        )
    )

    ps <- suppressWarnings(
        predict(
            fit,
            newdata = g,
            type = "response"
        )
    )

    if (any(!is.finite(ps))) stop("Non-finite subgroup PS")

    ## Only protect against exact numerical 0/1 generated by separation.
    ps <- pmin(
        pmax(ps, .Machine$double.eps),
        1 - .Machine$double.eps
    )

    g$ps_sub <- ps

    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    if (sum(g$ow_sub, na.rm = TRUE) <= 0) {
        stop("Invalid subgroup overlap weights")
    }

    g
}

estimate_one <- function(g, event = c("cancer", "other")) {
    event <- match.arg(event)

    g <- fit_subgroup_ow(g)

    a1 <- g[g$treat == 1, , drop = FALSE]
    a0 <- g[g$treat == 0, , drop = FALSE]

    if (event == "cancer") {
        c1 <- aj_w(
            a1$os_months,
            a1$cancer,
            a1$otherd,
            a1$ow_sub
        )
        c0 <- aj_w(
            a0$os_months,
            a0$cancer,
            a0$otherd,
            a0$ow_sub
        )
    } else {
        c1 <- aj_w(
            a1$os_months,
            a1$otherd,
            a1$cancer,
            a1$ow_sub
        )
        c0 <- aj_w(
            a0$os_months,
            a0$otherd,
            a0$cancer,
            a0$ow_sub
        )
    }

    c(
        limited = c0,
        colectomy = c1,
        rd = c1 - c0
    )
}

estimate_all <- function(d) {
    z <- make_outcome(d)
    masks <- make_masks(z)

    cancer <- rep(NA_real_, 7)
    other <- rep(NA_real_, 7)

    names(cancer) <- names(masks)
    names(other) <- names(masks)

    for (i in seq_along(masks)) {
        g <- z[masks[[i]], , drop = FALSE]

        cancer[i] <- tryCatch(
            estimate_one(g, "cancer")["rd"],
            error = function(e) NA_real_
        )

        other[i] <- tryCatch(
            estimate_one(g, "other")["rd"],
            error = function(e) NA_real_
        )
    }

    c(
        cancer,
        other
    )
}

## =============================================================================
## 4. Point estimates
## =============================================================================

d0 <- pipeline(RAW)
point <- estimate_all(d0)

labs <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

point_cancer <- point[1:7]
point_other <- point[8:14]

## =============================================================================
## 5. Balance diagnostics
## =============================================================================

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

    if (sum(tr == 1) == 0 || sum(tr == 0) == 0) {
        return(NA_real_)
    }

    m1 <- weighted.mean(x[tr == 1], w[tr == 1])
    m0 <- weighted.mean(x[tr == 0], w[tr == 0])

    v1 <- weighted.mean(
        (x[tr == 1] - m1)^2,
        w[tr == 1]
    )

    v0 <- weighted.mean(
        (x[tr == 0] - m0)^2,
        w[tr == 0]
    )

    denom <- sqrt((v1 + v0) / 2)

    if (!is.finite(denom) || denom <= 0) {
        return(0)
    }

    (m1 - m0) / denom
}

balance_matrix <- function(d) {
    grade_num <- suppressWarnings(as.numeric(d$grade))

    data.frame(
        age = d$age,
        age2 = d$age^2,
        female = d$female,
        year = d$year,
        logsize = d$logsize,
        size_miss = d$size_miss,

        T2 = as.numeric(d$T == "T2"),
        T3 = as.numeric(d$T == "T3"),
        T4 = as.numeric(d$T == "T4"),

        grade1 = as.numeric(grade_num == 1),
        grade2 = as.numeric(grade_num == 2),
        grade3 = as.numeric(grade_num == 3),
        grade4 = as.numeric(grade_num == 4),

        nonmucinous = as.numeric(d$hist_group == "nonmucinous"),
        mucinous = as.numeric(d$hist_group == "mucinous"),
        SRCC = as.numeric(d$hist_group == "SRCC"),

        race_black = as.numeric(d$race4 == "Non-Hispanic Black"),
        race_white = as.numeric(d$race4 == "Non-Hispanic White"),
        race_other = as.numeric(d$race4 == "Other/Unknown"),

        income_hi = d$income_hi,
        metro = d$metro,
        married = d$married,

        nodal_risk = d$nodal_risk,
        stringsAsFactors = FALSE
    )
}

z0 <- make_outcome(d0)
masks0 <- make_masks(z0)

diag_list <- list()
balance_list <- list()

for (nm in names(masks0)) {
    g <- z0[masks0[[nm]], , drop = FALSE]
    gs <- fit_subgroup_ow(g)

    a0 <- gs[gs$treat == 0, , drop = FALSE]
    a1 <- gs[gs$treat == 1, , drop = FALSE]

    bm <- balance_matrix(gs)

    smd_vec <- vapply(
        names(bm),
        function(v) {
            weighted_smd(
                bm[[v]],
                gs$treat,
                gs$ow_sub
            )
        },
        numeric(1)
    )

    cancer_pair <- estimate_one(g, "cancer")
    other_pair <- estimate_one(g, "other")

    q0 <- quantile(
        a0$ps_sub,
        c(.01, .99),
        na.rm = TRUE
    )

    q1 <- quantile(
        a1$ps_sub,
        c(.01, .99),
        na.rm = TRUE
    )

    diag_list[[nm]] <- data.frame(
        analysis = nm,
        n = nrow(gs),
        limited = nrow(a0),
        colectomy = nrow(a1),

        cancer_events_limited = sum(a0$cancer == 1),
        cancer_events_colectomy = sum(a1$cancer == 1),

        cancer_CIF_limited_pct = cancer_pair["limited"] * 100,
        cancer_CIF_colectomy_pct = cancer_pair["colectomy"] * 100,
        cancer_RD_pp = cancer_pair["rd"] * 100,

        other_CIF_limited_pct = other_pair["limited"] * 100,
        other_CIF_colectomy_pct = other_pair["colectomy"] * 100,
        other_RD_pp = other_pair["rd"] * 100,

        ESS_limited = ess(a0$ow_sub),
        ESS_colectomy = ess(a1$ow_sub),

        PS_P1_limited = q0[1],
        PS_P99_limited = q0[2],
        PS_P1_colectomy = q1[1],
        PS_P99_colectomy = q1[2],

        max_weight = max(gs$ow_sub, na.rm = TRUE),

        max_abs_weighted_SMD = max(
            abs(smd_vec),
            na.rm = TRUE
        ),

        worst_balance_variable = names(
            which.max(
                abs(smd_vec)
            )
        ),

        balance_pass = max(
            abs(smd_vec),
            na.rm = TRUE
        ) < 0.10,

        stringsAsFactors = FALSE
    )

    balance_list[[nm]] <- data.frame(
        analysis = nm,
        variable = names(smd_vec),
        weighted_SMD = as.numeric(smd_vec),
        stringsAsFactors = FALSE
    )
}

diagnostics <- do.call(
    rbind,
    diag_list
)

balance_components <- do.call(
    rbind,
    balance_list
)

## =============================================================================
## 6. Full-pipeline bootstrap B=1000
## =============================================================================

set.seed(20260901)

B <- 1000

boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = 14
)

colnames(boot) <- c(
    paste0("cancer_", labs),
    paste0("other_", labs)
)

start_time <- Sys.time()

for (b in seq_len(B)) {
    db <- RAW[
        sample(
            seq_len(nrow(RAW)),
            size = nrow(RAW),
            replace = TRUE
        ),
        ,
        drop = FALSE
    ]

    boot[b, ] <- tryCatch(
        {
            dx <- pipeline(db)
            estimate_all(dx)
        },
        error = function(e) {
            rep(NA_real_, 14)
        }
    )

    if (b %% 100 == 0) {
        cat(
            "STEP 1B bootstrap ",
            b,
            "/",
            B,
            " | elapsed ",
            round(
                difftime(
                    Sys.time(),
                    start_time,
                    units = "mins"
                ),
                1
            ),
            " min\n",
            sep = ""
        )
    }
}

## =============================================================================
## 7. Cancer-specific final table
## =============================================================================

result_list <- list()

for (i in seq_along(labs)) {
    v <- boot[, i]

    valid <- is.finite(v)

    result_list[[i]] <- data.frame(
        analysis = labs[i],
        RD = point_cancer[i],
        RD_pp = point_cancer[i] * 100,
        lo = quantile(
            v[valid],
            .025,
            na.rm = TRUE
        ),
        hi = quantile(
            v[valid],
            .975,
            na.rm = TRUE
        ),
        effective_B = sum(valid),
        stringsAsFactors = FALSE
    )
}

results <- do.call(
    rbind,
    result_list
)

results$lo_pp <- results$lo * 100
results$hi_pp <- results$hi * 100

results$upper_possible_reduction_pp <- pmax(
    0,
    -results$lo_pp
)

lowrisk <- c(
    "Q1",
    "bottom_decile",
    "risk_lt_5pct"
)

results$criterion_2pp <- ifelse(
    results$analysis %in% lowrisk,
    results$upper_possible_reduction_pp < 2,
    NA
)

results$criterion_3pp <- ifelse(
    results$analysis %in% lowrisk,
    results$upper_possible_reduction_pp < 3,
    NA
)

results$criterion_5pp <- ifelse(
    results$analysis %in% lowrisk,
    results$upper_possible_reduction_pp < 5,
    NA
)

## Other-cause diagnostic CIs

other_list <- list()

for (i in seq_along(labs)) {
    v <- boot[, 7 + i]
    valid <- is.finite(v)

    other_list[[i]] <- data.frame(
        analysis = labs[i],
        other_RD_pp = point_other[i] * 100,
        other_lo_pp = quantile(
            v[valid],
            .025,
            na.rm = TRUE
        ) * 100,
        other_hi_pp = quantile(
            v[valid],
            .975,
            na.rm = TRUE
        ) * 100,
        effective_B = sum(valid),
        stringsAsFactors = FALSE
    )
}

other_results <- do.call(
    rbind,
    other_list
)

results <- merge(
    results,
    other_results,
    by = "analysis",
    all.x = TRUE,
    sort = FALSE
)

results <- results[
    match(
        labs,
        results$analysis
    ),
    ,
    drop = FALSE
]

## =============================================================================
## 8. Compare with Step 1A global no-stage OW
## =============================================================================

comparison <- NULL

if (file.exists(STEP1A_FILE)) {
    step1a <- read.csv(
        STEP1A_FILE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    step1a <- step1a[
        step1a$analysis %in% labs,
        ,
        drop = FALSE
    ]

    comparison <- merge(
        data.frame(
            analysis = results$analysis,
            subgroup_OW_RD_pp = results$RD_pp,
            subgroup_OW_lo_pp = results$lo_pp,
            subgroup_OW_hi_pp = results$hi_pp,
            stringsAsFactors = FALSE
        ),
        data.frame(
            analysis = step1a$analysis,
            global_no_stage_OW_RD_pp = step1a$rd_pp,
            global_no_stage_OW_lo_pp = step1a$lo_pp,
            global_no_stage_OW_hi_pp = step1a$hi_pp,
            stringsAsFactors = FALSE
        ),
        by = "analysis",
        all.x = TRUE,
        sort = FALSE
    )

    comparison <- comparison[
        match(
            labs,
            comparison$analysis
        ),
        ,
        drop = FALSE
    ]

    comparison$point_change_pp <- (
        comparison$subgroup_OW_RD_pp -
        comparison$global_no_stage_OW_RD_pp
    )
}

## =============================================================================
## 9. Save files
## =============================================================================

OUT_RESULT <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Primary_B1000.csv"
)

OUT_DIAG <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Diagnostics.csv"
)

OUT_BAL <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Balance_Components.csv"
)

OUT_BOOT <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Bootstrap_B1000.csv"
)

OUT_COMPARE <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Comparison.csv"
)

OUT_REPORT <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Report.txt"
)

write.csv(
    results,
    OUT_RESULT,
    row.names = FALSE
)

write.csv(
    diagnostics,
    OUT_DIAG,
    row.names = FALSE
)

write.csv(
    balance_components,
    OUT_BAL,
    row.names = FALSE
)

write.csv(
    data.frame(
        replicate = seq_len(B),
        boot,
        check.names = FALSE
    ),
    OUT_BOOT,
    row.names = FALSE
)

if (!is.null(comparison)) {
    write.csv(
        comparison,
        OUT_COMPARE,
        row.names = FALSE
    )
}

## =============================================================================
## 10. Human-readable report
## =============================================================================

report <- c(
    "======================================================================",
    "A6 RED-TEAM STEP 1B — RISK-STRATUM-SPECIFIC OVERLAP WEIGHTING",
    "======================================================================",
    paste0(
        "Removed PS term(s): ",
        paste(stage_terms, collapse = " | ")
    ),
    "Primary low-risk anchor remains bottom predicted-risk decile.",
    "risk<5% remains supportive.",
    "",
    "== CANCER-SPECIFIC MORTALITY =="
)

for (i in seq_len(nrow(results))) {
    z <- results[i, ]

    extra <- ""

    if (z$analysis %in% lowrisk) {
        extra <- sprintf(
            " | upper reduction %.1f pp | 2pp:%s 3pp:%s 5pp:%s",
            z$upper_possible_reduction_pp,
            ifelse(z$criterion_2pp, "MET", "NOT MET"),
            ifelse(z$criterion_3pp, "MET", "NOT MET"),
            ifelse(z$criterion_5pp, "MET", "NOT MET")
        )
    }

    report <- c(
        report,
        paste0(
            sprintf(
                "%-15s RD %+.1f pp (95%% CI %+.1f to %+.1f; B=%d)",
                z$analysis,
                z$RD_pp,
                z$lo_pp,
                z$hi_pp,
                z$effective_B.x
            ),
            extra
        )
    )
}

report <- c(
    report,
    "",
    "== OTHER-CAUSE MORTALITY RESIDUAL-SELECTION DIAGNOSTIC =="
)

for (i in seq_len(nrow(results))) {
    z <- results[i, ]

    report <- c(
        report,
        sprintf(
            "%-15s RD %+.1f pp (95%% CI %+.1f to %+.1f; B=%d)",
            z$analysis,
            z$other_RD_pp,
            z$other_lo_pp,
            z$other_hi_pp,
            z$effective_B.y
        )
    )
}

report <- c(
    report,
    "",
    "== SUBGROUP-SPECIFIC BALANCE / SUPPORT =="
)

for (i in seq_len(nrow(diagnostics))) {
    z <- diagnostics[i, ]

    report <- c(
        report,
        sprintf(
            paste0(
                "%-15s n=%d (%d/%d) events=%d/%d | ",
                "CSM CIF %.1f/%.1f%% RD %+.1f | ESS %.1f/%.1f | ",
                "max weight %.3f | max|SMD| %.3f [%s] | %s"
            ),
            z$analysis,
            z$n,
            z$limited,
            z$colectomy,
            z$cancer_events_limited,
            z$cancer_events_colectomy,
            z$cancer_CIF_limited_pct,
            z$cancer_CIF_colectomy_pct,
            z$cancer_RD_pp,
            z$ESS_limited,
            z$ESS_colectomy,
            z$max_weight,
            z$max_abs_weighted_SMD,
            z$worst_balance_variable,
            ifelse(
                z$balance_pass,
                "PASS",
                "FAIL"
            )
        )
    )
}

if (!is.null(comparison)) {
    report <- c(
        report,
        "",
        "== CHANGE VS STEP 1A GLOBAL NO-STAGE OW =="
    )

    for (i in seq_len(nrow(comparison))) {
        z <- comparison[i, ]

        report <- c(
            report,
            sprintf(
                "%-15s global %+.2f pp | subgroup-OW %+.2f pp | change %+.2f pp",
                z$analysis,
                z$global_no_stage_OW_RD_pp,
                z$subgroup_OW_RD_pp,
                z$point_change_pp
            )
        )
    }
}

report <- c(
    report,
    "",
    "== DECISION RULE ==",
    "If bottom-decile subgroup-specific weighting preserves balance and the 3-pp criterion, the core low-risk conclusion survives the strongest current exchangeability challenge.",
    "If the bottom-decile 3-pp criterion is lost, the manuscript's result-convincingness score must be downgraded and the headline reframed.",
    "======================================================================"
)

writeLines(
    report,
    OUT_REPORT,
    useBytes = TRUE
)

cat(
    paste(
        report,
        collapse = "\n"
    ),
    "\n"
)

cat(
    "\nREPORT: ",
    OUT_REPORT,
    "\n",
    sep = ""
)
