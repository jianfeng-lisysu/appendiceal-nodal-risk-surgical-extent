## =============================================================================
## A6 FINAL PRESUBMISSION — STEP 07 v2
## CORRECTED-FRAMEWORK LANDMARK (3/6/12 mo) + HISTOLOGY x RISK, B=1000
##
## PURPOSE
##   Final presubmission rerun of:
##     A) 3/6/12-month landmark sensitivity
##     B) histology x predicted-risk-half sensitivity
##
## LOCKED FRAMEWORK
##   - Summary Stage is excluded from every propensity-score model.
##   - Primary risk-stratified analyses refit PS/OW within each risk stratum.
##   - Histology subgroups use a histology-free, no-Summary-Stage PS.
##   - Full-pipeline bootstrap: raw resampling -> frozen pipeline / nodal-risk
##     model refit -> risk definitions rebuilt -> subgroup PS/OW refit ->
##     weighted Aalen-Johansen RD.
##   - B=1000; do not change B based on whether a decision bound is met.
##
## IMPORTANT CORRECTIONS VERSUS EARLIER STEP-07 DRAFT
##   1) Landmark analyses preserve the PRIMARY baseline risk definition:
##      percentile cut points/P10 are estimated in the full outcome-evaluable
##      bootstrap sample BEFORE restricting to survivors at L. The landmark
##      therefore evaluates the same baseline risk strata after delayed entry,
##      rather than redefining "bottom decile" among L-month survivors.
##   2) Histology risk halves preserve the legacy/frozen definition:
##      overall median predicted nodal risk is estimated in the full
##      risk-scored analytic cohort BEFORE outcome filtering.
##   3) Histology PS uses a full-rank design-matrix fitter to avoid artificial
##      failure from constant/aliased factor columns inside sparse subgroups.
##   4) Landmark results include n, arm counts, ESS, max weighted |SMD|,
##      2/3/5-pp decision-bound results, and effective bootstrap B.
##   5) Histology estimates require effective B >=900 to be labeled estimable.
##   6) Bootstrap checkpoint/resume is enabled.
##
## OUTPUTS (03_results/)
##   PRESUB_07_Landmark_B1000.csv
##   PRESUB_07_Landmark_Boot_B1000.csv
##   PRESUB_07_Histology_B1000.csv
##   PRESUB_07_Histology_Boot_B1000.csv
##   PRESUB_07_Report.txt
##   PRESUB_07_Audit.txt
##   PRESUB_07_checkpoints/*.rds
##
## INTERPRETATION
##   Landmark analyses assess sensitivity to early postdiagnosis guarantee-time
##   effects. Because SEER lacks exact definitive-surgery dates, they cannot
##   eliminate immortal-time bias.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

SCRIPT_VERSION <- "A6_PRESUB07_v2_1_20260906_SYNTAXFIX"
B_BOOT <- 1000L
MIN_EFFECTIVE_B <- 900L

ROOT_CANDIDATES <- c(
    ".",
    "."
)

ROOT <- ROOT_CANDIDATES[dir.exists(ROOT_CANDIDATES)][1]

if (length(ROOT) == 0 || is.na(ROOT)) {
    stop(
        "A6 root not found. Expected one of: ",
        paste(ROOT_CANDIDATES, collapse = " | ")
    )
}

DATA_DIR <- file.path(ROOT, "02_data")
RESULT_DIR <- file.path(ROOT, "03_results")
CHECK_DIR <- file.path(RESULT_DIR, "PRESUB_07_checkpoints")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)

PIPE_FILE <- file.path(DATA_DIR, "Aline_frozen_pipeline_v2.RData")
RAW_FILE <- file.path(DATA_DIR, "analytic_cohort_A.csv")
GATE_FILE <- file.path(
    RESULT_DIR,
    "REDTEAM_01B_SubgroupOW_Primary_B1000.csv"
)

for (f in c(PIPE_FILE, RAW_FILE, GATE_FILE)) {
    if (!file.exists(f)) {
        stop("Missing required file: ", f)
    }
}

load(PIPE_FILE)

if (!exists("pipeline", mode = "function")) {
    stop("pipeline() missing from frozen RData.")
}

if (!exists("PS_FORMULA")) {
    stop("PS_FORMULA missing from frozen RData.")
}

RAW <- read.csv(
    RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(RAW) != 6951) {
    stop(
        "Frozen cohort fingerprint failed: N must be 6951; current N=",
        nrow(RAW)
    )
}


## =============================================================================
## 1. Corrected PS specifications
## =============================================================================

PS_ORIG <- PS_FORMULA

terms_orig <- attr(
    terms(PS_ORIG),
    "term.labels"
)

stage_terms <- terms_orig[
    grepl(
        "stage",
        terms_orig,
        ignore.case = TRUE
    )
]

if (length(stage_terms) == 0) {
    stop("No stage-related term found in PS_FORMULA.")
}

PS_NOSTAGE <- reformulate(
    termlabels = setdiff(
        terms_orig,
        stage_terms
    ),
    response = all.vars(PS_ORIG)[1],
    env = environment(PS_ORIG)
)

pipeline_env <- environment(pipeline)

assign(
    "PS_FORMULA",
    PS_NOSTAGE,
    envir = .GlobalEnv
)

assign(
    "PS_FORMULA",
    PS_NOSTAGE,
    envir = pipeline_env
)

terms_ns <- attr(
    terms(PS_NOSTAGE),
    "term.labels"
)

hist_terms <- terms_ns[
    grepl(
        "hist",
        terms_ns,
        ignore.case = TRUE
    )
]

if (length(hist_terms) == 0) {
    stop("No histology-related term found in corrected PS formula.")
}

PS_NOHIST <- reformulate(
    termlabels = setdiff(
        terms_ns,
        hist_terms
    ),
    response = all.vars(PS_NOSTAGE)[1],
    env = environment(PS_NOSTAGE)
)


## =============================================================================
## 2. Common estimators — Step-1B conventions
## =============================================================================

aj_w <- function(
    time,
    ev1,
    ev2,
    wt,
    start = -1,
    tau = 60
) {

    ok <- (
        !is.na(time) &
        !is.na(ev1) &
        !is.na(ev2) &
        !is.na(wt) &
        is.finite(wt)
    )

    time <- time[ok]
    ev1 <- ev1[ok]
    ev2 <- ev2[ok]
    wt <- wt[ok]

    keep <- time > start

    time <- time[keep]
    ev1 <- ev1[keep]
    ev2 <- ev2[keep]
    wt <- wt[keep]

    if (
        length(time) == 0 ||
        sum(wt) <= 0
    ) {
        return(NA_real_)
    }

    dc <- rowsum(
        wt * ev1,
        time,
        reorder = FALSE
    )

    dr <- rowsum(
        wt * ev2,
        time,
        reorder = FALSE
    )

    dw <- rowsum(
        wt,
        time,
        reorder = FALSE
    )

    ut <- as.numeric(
        rownames(dc)
    )

    ord <- order(ut)

    ut <- ut[ord]
    dc <- dc[ord]
    dr <- dr[ord]
    dw <- dw[ord]

    S <- 1
    cif <- 0
    at_risk <- sum(wt)

    for (k in seq_along(ut)) {

        if (ut[k] > tau) {
            break
        }

        if (at_risk > 0) {

            cif <- cif +
                S *
                dc[k] /
                at_risk

            S <- S *
                (
                    1 -
                    (
                        dc[k] +
                        dr[k]
                    ) /
                    at_risk
                )
        }

        at_risk <- at_risk -
            dw[k]
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

    cuts <- as.numeric(
        quantile(
            x,
            probs = 0:5 / 5,
            na.rm = TRUE,
            names = FALSE
        )
    )

    cuts[1] <- cuts[1] - 1e-8
    cuts[6] <- cuts[6] + 1e-8

    if (
        any(!is.finite(cuts)) ||
        any(diff(cuts) <= 0)
    ) {
        stop("Invalid predicted-risk quintile cut points.")
    }

    cuts
}


make_thresholds <- function(reference_risk) {

    list(
        cuts = make_cuts(reference_risk),
        p10 = as.numeric(
            quantile(
                reference_risk,
                0.10,
                na.rm = TRUE
            )
        )
    )
}


make_masks_with_thresholds <- function(
    d,
    thresholds
) {

    grp <- cut(
        d$nodal_risk,
        breaks = thresholds$cuts,
        labels = FALSE,
        include.lowest = TRUE
    )

    list(
        Q1 = grp == 1 & !is.na(grp),
        Q2 = grp == 2 & !is.na(grp),
        Q3 = grp == 3 & !is.na(grp),
        Q4 = grp == 4 & !is.na(grp),
        Q5 = grp == 5 & !is.na(grp),
        bottom_decile = (
            !is.na(d$nodal_risk) &
            d$nodal_risk <= thresholds$p10
        ),
        risk_lt_5pct = (
            !is.na(d$nodal_risk) &
            d$nodal_risk < 0.05
        )
    )
}


make_masks <- function(d) {

    make_masks_with_thresholds(
        d,
        make_thresholds(
            d$nodal_risk
        )
    )
}


fit_subgroup_ow <- function(
    g,
    form = PS_NOSTAGE
) {

    if (nrow(g) < 40) {
        stop("Subgroup too small.")
    }

    if (
        sum(g$treat == 0, na.rm = TRUE) <
        15
    ) {
        stop("Too few limited-resection patients.")
    }

    if (
        sum(g$treat == 1, na.rm = TRUE) <
        15
    ) {
        stop("Too few colectomy patients.")
    }

    fit <- suppressWarnings(
        glm(
            form,
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

    if (
        length(ps) != nrow(g) ||
        any(!is.finite(ps))
    ) {
        stop("Non-finite or length-mismatched subgroup PS.")
    }

    ps <- pmin(
        pmax(
            ps,
            .Machine$double.eps
        ),
        1 - .Machine$double.eps
    )

    g$ps_sub <- ps

    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    if (
        any(!is.finite(g$ow_sub)) ||
        sum(g$ow_sub, na.rm = TRUE) <= 0
    ) {
        stop("Invalid subgroup overlap weights.")
    }

    g
}


estimate_rd <- function(
    g,
    form = PS_NOSTAGE,
    start = -1
) {

    gs <- fit_subgroup_ow(
        g,
        form
    )

    a1 <- gs[
        gs$treat == 1,
        ,
        drop = FALSE
    ]

    a0 <- gs[
        gs$treat == 0,
        ,
        drop = FALSE
    ]

    c1 <- aj_w(
        a1$os_months,
        a1$cancer,
        a1$otherd,
        a1$ow_sub,
        start = start,
        tau = 60
    )

    c0 <- aj_w(
        a0$os_months,
        a0$cancer,
        a0$otherd,
        a0$ow_sub,
        start = start,
        tau = 60
    )

    c1 - c0
}


ess <- function(w) {

    w <- w[
        is.finite(w) &
        w > 0
    ]

    if (length(w) == 0) {
        return(NA_real_)
    }

    sum(w)^2 /
        sum(w^2)
}


weighted_smd <- function(
    x,
    tr,
    w
) {

    ok <- (
        is.finite(x) &
        is.finite(tr) &
        is.finite(w)
    )

    x <- x[ok]
    tr <- tr[ok]
    w <- w[ok]

    if (
        sum(tr == 0) == 0 ||
        sum(tr == 1) == 0
    ) {
        return(NA_real_)
    }

    m0 <- weighted.mean(
        x[tr == 0],
        w[tr == 0]
    )

    m1 <- weighted.mean(
        x[tr == 1],
        w[tr == 1]
    )

    v0 <- weighted.mean(
        (x[tr == 0] - m0)^2,
        w[tr == 0]
    )

    v1 <- weighted.mean(
        (x[tr == 1] - m1)^2,
        w[tr == 1]
    )

    den <- sqrt(
        (v0 + v1) / 2
    )

    if (
        !is.finite(den) ||
        den <= 0
    ) {
        return(0)
    }

    (m1 - m0) /
        den
}


rhs_matrix <- function(
    form,
    data
) {

    X <- model.matrix(
        delete.response(
            terms(form)
        ),
        data = data,
        na.action = na.pass
    )

    if (
        "(Intercept)" %in%
        colnames(X)
    ) {
        X <- X[
            ,
            colnames(X) != "(Intercept)",
            drop = FALSE
        ]
    }

    X
}


balance_from_weighted_data <- function(
    gs,
    form,
    include_nodal_risk = TRUE
) {

    X <- rhs_matrix(
        form,
        gs
    )

    if (include_nodal_risk) {
        X <- cbind(
            X,
            nodal_risk = gs$nodal_risk
        )
    }

    if (ncol(X) == 0) {
        return(
            list(
                max_abs_smd = NA_real_,
                worst_variable = NA_character_
            )
        )
    }

    smd <- vapply(
        seq_len(ncol(X)),
        function(j) {
            weighted_smd(
                X[, j],
                gs$treat,
                gs$ow_sub
            )
        },
        numeric(1)
    )

    names(smd) <- colnames(X)

    finite <- is.finite(smd)

    if (!any(finite)) {
        return(
            list(
                max_abs_smd = NA_real_,
                worst_variable = NA_character_
            )
        )
    }

    abs_smd <- abs(smd)

    list(
        max_abs_smd = max(
            abs_smd[finite],
            na.rm = TRUE
        ),
        worst_variable = names(
            which.max(
                abs_smd
            )
        )
    )
}


diagnose_standard_subgroup <- function(
    g,
    form = PS_NOSTAGE
) {

    gs <- fit_subgroup_ow(
        g,
        form
    )

    a0 <- gs[
        gs$treat == 0,
        ,
        drop = FALSE
    ]

    a1 <- gs[
        gs$treat == 1,
        ,
        drop = FALSE
    ]

    bal <- balance_from_weighted_data(
        gs,
        form,
        include_nodal_risk = TRUE
    )

    data.frame(
        n = nrow(gs),
        n_limited = nrow(a0),
        n_colectomy = nrow(a1),
        ESS_limited = ess(
            a0$ow_sub
        ),
        ESS_colectomy = ess(
            a1$ow_sub
        ),
        max_weight = max(
            gs$ow_sub,
            na.rm = TRUE
        ),
        max_abs_wSMD = bal$max_abs_smd,
        worst_balance_variable =
            bal$worst_variable,
        balance_pass = (
            is.finite(
                bal$max_abs_smd
            ) &&
            bal$max_abs_smd < 0.10
        ),
        stringsAsFactors = FALSE
    )
}


## =============================================================================
## 3. Full-rank histology-subgroup PS
## =============================================================================

make_full_rank_matrix <- function(
    g,
    form
) {

    X0 <- model.matrix(
        delete.response(
            terms(form)
        ),
        data = g,
        na.action = na.pass
    )

    if (nrow(X0) != nrow(g)) {
        stop("Histology PS design matrix changed row count.")
    }

    keep <- rep(
        TRUE,
        ncol(X0)
    )

    names(keep) <- colnames(X0)

    for (j in seq_len(ncol(X0))) {

        nm <- colnames(X0)[j]

        if (nm == "(Intercept)") {
            next
        }

        x <- X0[, j]

        if (
            all(!is.finite(x)) ||
            var(
                x,
                na.rm = TRUE
            ) == 0
        ) {
            keep[j] <- FALSE
        }
    }

    X1 <- X0[
        ,
        keep,
        drop = FALSE
    ]

    selected <- integer(0)

    int_pos <- match(
        "(Intercept)",
        colnames(X1)
    )

    if (!is.na(int_pos)) {
        selected <- int_pos
    }

    candidates <- setdiff(
        seq_len(ncol(X1)),
        selected
    )

    for (j in candidates) {

        trial <- unique(
            c(
                selected,
                j
            )
        )

        if (
            qr(
                X1[
                    ,
                    trial,
                    drop = FALSE
                ],
                tol = 1e-10
            )$rank ==
            length(trial)
        ) {
            selected <- trial
        }
    }

    if (length(selected) < 2) {
        stop(
            "Histology subgroup PS matrix has insufficient rank."
        )
    }

    X1[
        ,
        selected,
        drop = FALSE
    ]
}


fit_histology_ow <- function(g) {

    if (nrow(g) < 40) {
        stop("Histology subgroup too small.")
    }

    if (
        sum(g$treat == 0, na.rm = TRUE) <
        15
    ) {
        stop("Too few limited-resection patients.")
    }

    if (
        sum(g$treat == 1, na.rm = TRUE) <
        15
    ) {
        stop("Too few colectomy patients.")
    }

    X <- make_full_rank_matrix(
        g,
        PS_NOHIST
    )

    fit <- suppressWarnings(
        glm.fit(
            x = X,
            y = g$treat,
            family = binomial()
        )
    )

    beta <- fit$coefficients

    if (
        length(beta) != ncol(X) ||
        any(!is.finite(beta))
    ) {
        stop("Non-finite histology-subgroup PS coefficient.")
    }

    ps <- plogis(
        as.numeric(
            X %*%
            beta
        )
    )

    if (
        length(ps) != nrow(g) ||
        any(!is.finite(ps))
    ) {
        stop("Non-finite histology-subgroup PS.")
    }

    ps <- pmin(
        pmax(
            ps,
            .Machine$double.eps
        ),
        1 - .Machine$double.eps
    )

    g$ps_sub <- ps

    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    if (
        any(!is.finite(g$ow_sub)) ||
        sum(g$ow_sub, na.rm = TRUE) <= 0
    ) {
        stop("Invalid histology-subgroup overlap weights.")
    }

    g
}


estimate_histology_rd <- function(g) {

    gs <- fit_histology_ow(
        g
    )

    a1 <- gs[
        gs$treat == 1,
        ,
        drop = FALSE
    ]

    a0 <- gs[
        gs$treat == 0,
        ,
        drop = FALSE
    ]

    c1 <- aj_w(
        a1$os_months,
        a1$cancer,
        a1$otherd,
        a1$ow_sub,
        start = -1,
        tau = 60
    )

    c0 <- aj_w(
        a0$os_months,
        a0$cancer,
        a0$otherd,
        a0$ow_sub,
        start = -1,
        tau = 60
    )

    c1 - c0
}


diagnose_histology_subgroup <- function(g) {

    gs <- fit_histology_ow(
        g
    )

    a0 <- gs[
        gs$treat == 0,
        ,
        drop = FALSE
    ]

    a1 <- gs[
        gs$treat == 1,
        ,
        drop = FALSE
    ]

    bal <- balance_from_weighted_data(
        gs,
        PS_NOHIST,
        include_nodal_risk = TRUE
    )

    data.frame(
        n = nrow(gs),
        n_limited = nrow(a0),
        n_colectomy = nrow(a1),
        ESS_limited = ess(
            a0$ow_sub
        ),
        ESS_colectomy = ess(
            a1$ow_sub
        ),
        max_weight = max(
            gs$ow_sub,
            na.rm = TRUE
        ),
        max_abs_wSMD = bal$max_abs_smd,
        worst_balance_variable =
            bal$worst_variable,
        balance_pass = (
            is.finite(
                bal$max_abs_smd
            ) &&
            bal$max_abs_smd < 0.10
        ),
        stringsAsFactors = FALSE
    )
}


## =============================================================================
## 4. Step-1B fingerprint gate
## =============================================================================

labs7 <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

d0 <- pipeline(
    RAW
)

z_gate <- make_outcome(
    d0
)

gate_masks <- make_masks(
    z_gate
)

gate_point <- vapply(
    labs7,
    function(nm) {

        tryCatch(
            estimate_rd(
                z_gate[
                    gate_masks[[nm]],
                    ,
                    drop = FALSE
                ],
                form = PS_NOSTAGE,
                start = -1
            ),
            error = function(e) {
                NA_real_
            }
        )
    },
    numeric(1)
)

gate_ref <- read.csv(
    GATE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (
    !all(
        c(
            "analysis",
            "RD"
        ) %in%
        names(gate_ref)
    )
) {
    stop(
        "Step-1B gate file must contain columns: analysis, RD."
    )
}

gate_ref_vec <- gate_ref[
    match(
        labs7,
        gate_ref$analysis
    ),
    "RD"
]

if (
    any(
        !is.finite(
            gate_ref_vec
        )
    )
) {
    stop(
        "Step-1B gate file does not contain finite RD values for all seven analyses."
    )
}

gate_diff <- max(
    abs(
        gate_point -
        gate_ref_vec
    )
)

cat(
    sprintf(
        "FINGERPRINT GATE: max |diff| vs Step 1B = %.3e\n",
        gate_diff
    )
)

if (
    !is.finite(gate_diff) ||
    gate_diff > 1e-9
) {
    stop(
        "FINGERPRINT GATE FAILED — do NOT proceed."
    )
}

cat(
    "FINGERPRINT GATE PASSED.\n\n"
)


## =============================================================================
## 5. Bootstrap checkpoint helper
## =============================================================================

run_bootstrap <- function(
    label,
    ncol_out,
    colnames_out,
    seed,
    estimator_fun
) {

    checkpoint_file <- file.path(
        CHECK_DIR,
        paste0(
            label,
            "_B1000.rds"
        )
    )

    boot <- matrix(
        NA_real_,
        nrow = B_BOOT,
        ncol = ncol_out,
        dimnames = list(
            NULL,
            colnames_out
        )
    )

    start_b <- 1L

    if (
        file.exists(
            checkpoint_file
        )
    ) {

        cp <- tryCatch(
            readRDS(
                checkpoint_file
            ),
            error = function(e) {
                NULL
            }
        )

        cp_ok <- (
            !is.null(cp) &&
            identical(
                cp$script_version,
                SCRIPT_VERSION
            ) &&
            identical(
                cp$label,
                label
            ) &&
            identical(
                cp$B,
                B_BOOT
            ) &&
            is.matrix(
                cp$boot
            ) &&
            identical(
                dim(cp$boot),
                dim(boot)
            ) &&
            identical(
                colnames(cp$boot),
                colnames(boot)
            )
        )

        if (cp_ok) {

            boot <- cp$boot

            start_b <- cp$completed_through +
                1L

            if (
                !is.null(
                    cp$rng_seed
                )
            ) {
                assign(
                    ".Random.seed",
                    cp$rng_seed,
                    envir = .GlobalEnv
                )
            }

            cat(
                label,
                ": checkpoint loaded; resume at ",
                start_b,
                "/",
                B_BOOT,
                "\n",
                sep = ""
            )
        }
    }

    if (start_b == 1L) {
        set.seed(
            seed
        )
    }

    if (start_b <= B_BOOT) {

        t0 <- Sys.time()

        for (
            b in seq.int(
                start_b,
                B_BOOT
            )
        ) {

            db <- RAW[
                sample(
                    seq_len(
                        nrow(RAW)
                    ),
                    size = nrow(RAW),
                    replace = TRUE
                ),
                ,
                drop = FALSE
            ]

            boot[
                b,
            ] <- tryCatch(
                estimator_fun(
                    db
                ),
                error = function(e) {
                    rep(
                        NA_real_,
                        ncol_out
                    )
                }
            )

            if (
                b %% 50 ==
                0 ||
                b ==
                B_BOOT
            ) {

                saveRDS(
                    list(
                        script_version =
                            SCRIPT_VERSION,
                        label = label,
                        B = B_BOOT,
                        completed_through = b,
                        boot = boot,
                        rng_seed = .Random.seed
                    ),
                    checkpoint_file
                )
            }

            if (
                b %% 100 ==
                0 ||
                b ==
                B_BOOT
            ) {

                cat(
                    label,
                    ": ",
                    b,
                    "/",
                    B_BOOT,
                    " | elapsed ",
                    sprintf(
                        "%.1f",
                        as.numeric(
                            difftime(
                                Sys.time(),
                                t0,
                                units = "mins"
                            )
                        )
                    ),
                    " min\n",
                    sep = ""
                )
            }
        }
    }

    boot
}


## =============================================================================
## 6. Report/audit initialization
## =============================================================================

audit <- c(
    "======================================================================",
    "A6 FINAL PRESUBMISSION STEP 07 v2 — AUDIT",
    "======================================================================",
    paste0(
        "Script version: ",
        SCRIPT_VERSION
    ),
    paste0(
        "Run time: ",
        format(
            Sys.time()
        )
    ),
    paste0(
        "ROOT: ",
        ROOT
    ),
    paste0(
        "R version: ",
        R.version.string
    ),
    sprintf(
        "Step-1B fingerprint gate: PASS; max |diff| = %.3e",
        gate_diff
    ),
    paste0(
        "Removed Summary-Stage term(s): ",
        paste(
            stage_terms,
            collapse = " | "
        )
    ),
    paste0(
        "Histology term(s) additionally removed in histology subgroup PS: ",
        paste(
            hist_terms,
            collapse = " | "
        )
    ),
    "Landmark risk definition: baseline/full outcome-evaluable bootstrap-sample quintiles and P10 are defined BEFORE landmark restriction.",
    "Landmark 0.05 threshold: fixed.",
    "Landmark horizon: conditional on surviving beyond L; outcome horizon remains 60 months post diagnosis.",
    "Histology risk-half definition: overall median predicted nodal risk in the full risk-scored analytic cohort, before outcome filtering.",
    paste0(
        "B=",
        B_BOOT,
        "; estimable-result minimum effective B=",
        MIN_EFFECTIVE_B,
        "."
    ),
    "Decision-bound sensitivity: 2/3/5 percentage points; 3 pp remains the prespecified primary clinical decision bound.",
    "======================================================================"
)

report <- c(
    "======================================================================",
    "A6 PRESUBMISSION STEP 07 v2 — LANDMARK + HISTOLOGY x RISK",
    "======================================================================",
    "Framework: no Summary Stage in PS; subgroup-specific OW; full-pipeline B=1000.",
    "RD = oncologic colectomy minus limited resection in 5-year cancer-specific mortality.",
    "Negative RD indicates lower cancer-specific mortality associated with colectomy.",
    ""
)


## =============================================================================
## 7. PART A — Landmark 3 / 6 / 12 months
## =============================================================================

landmark_masks <- function(
    d,
    L
) {

    z_full <- make_outcome(
        d
    )

    thresholds <- make_thresholds(
        z_full$nodal_risk
    )

    z_land <- z_full[
        z_full$os_months >
        L,
        ,
        drop = FALSE
    ]

    masks <- make_masks_with_thresholds(
        z_land,
        thresholds
    )

    list(
        z_full = z_full,
        z_land = z_land,
        thresholds = thresholds,
        masks = masks
    )
}


landmark_estimate_all <- function(
    d,
    L
) {

    obj <- landmark_masks(
        d,
        L
    )

    out <- setNames(
        rep(
            NA_real_,
            length(labs7)
        ),
        labs7
    )

    for (nm in labs7) {

        g <- obj$z_land[
            obj$masks[[nm]],
            ,
            drop = FALSE
        ]

        out[nm] <- tryCatch(
            estimate_rd(
                g,
                form = PS_NOSTAGE,
                start = L
            ),
            error = function(e) {
                NA_real_
            }
        )
    }

    out
}


landmark_point_diagnostics <- function(
    d,
    L
) {

    obj <- landmark_masks(
        d,
        L
    )

    rows <- list()

    for (nm in labs7) {

        g <- obj$z_land[
            obj$masks[[nm]],
            ,
            drop = FALSE
        ]

        dg <- tryCatch(
            diagnose_standard_subgroup(
                g,
                form = PS_NOSTAGE
            ),
            error = function(e) {
                data.frame(
                    n = nrow(g),
                    n_limited = sum(
                        g$treat == 0,
                        na.rm = TRUE
                    ),
                    n_colectomy = sum(
                        g$treat == 1,
                        na.rm = TRUE
                    ),
                    ESS_limited = NA_real_,
                    ESS_colectomy = NA_real_,
                    max_weight = NA_real_,
                    max_abs_wSMD = NA_real_,
                    worst_balance_variable =
                        NA_character_,
                    balance_pass = FALSE,
                    stringsAsFactors = FALSE
                )
            }
        )

        dg$landmark <- L
        dg$analysis <- nm

        rows[[nm]] <- dg
    }

    out <- do.call(
        rbind,
        rows
    )

    rownames(out) <- NULL

    out
}


landmark_rows <- list()
landmark_boot_list <- list()
landmark_row_counter <- 0L

LANDMARK_SEEDS <- c(
    `3` = 20260906L,
    `6` = 20260909L,
    `12` = 20260915L
)

for (L in c(3, 6, 12)) {

    cat(
        "\n============================================================\n"
    )

    cat(
        "LANDMARK ",
        L,
        " MONTHS\n",
        sep = ""
    )

    cat(
        "============================================================\n"
    )

    pt <- landmark_estimate_all(
        d0,
        L
    )

    dg <- landmark_point_diagnostics(
        d0,
        L
    )

    boot <- run_bootstrap(
        label = paste0(
            "Landmark_",
            L,
            "m"
        ),
        ncol_out = length(labs7),
        colnames_out = labs7,
        seed = LANDMARK_SEEDS[
            as.character(L)
        ],
        estimator_fun = function(db) {

            dx <- pipeline(
                db
            )

            landmark_estimate_all(
                dx,
                L
            )
        }
    )

    landmark_boot_list[[as.character(L)]] <- boot

    effB <- colSums(
        is.finite(
            boot
        )
    )

    low_anchor_eff <- effB[
        c(
            "bottom_decile",
            "risk_lt_5pct"
        )
    ]

    if (
        any(
            low_anchor_eff <
            MIN_EFFECTIVE_B
        )
    ) {
        stop(
            "Landmark ",
            L,
            " months: effective B < ",
            MIN_EFFECTIVE_B,
            " for a prespecified low-risk anchor."
        )
    }

    report <- c(
        report,
        paste0(
            "== LANDMARK ",
            L,
            " MONTHS ==",
            " full landmark population n=",
            nrow(
                landmark_masks(
                    d0,
                    L
                )$z_land
            )
        )
    )

    for (nm in labs7) {

        v <- boot[
            ,
            nm
        ]

        valid <- is.finite(
            v
        )

        lo <- if (
            sum(valid) >= 50
        ) {
            as.numeric(
                quantile(
                    v[valid],
                    0.025,
                    na.rm = TRUE
                )
            )
        } else {
            NA_real_
        }

        hi <- if (
            sum(valid) >= 50
        ) {
            as.numeric(
                quantile(
                    v[valid],
                    0.975,
                    na.rm = TRUE
                )
            )
        } else {
            NA_real_
        }

        upper_reduction <- if (
            is.finite(lo)
        ) {
            max(
                0,
                -lo * 100
            )
        } else {
            NA_real_
        }

        is_decision_anchor <- nm %in%
            c(
                "bottom_decile",
                "risk_lt_5pct"
            )

        criterion_2pp <- if (
            is_decision_anchor &&
            is.finite(
                upper_reduction
            )
        ) {
            upper_reduction < 2
        } else {
            NA
        }

        criterion_3pp <- if (
            is_decision_anchor &&
            is.finite(
                upper_reduction
            )
        ) {
            upper_reduction < 3
        } else {
            NA
        }

        criterion_5pp <- if (
            is_decision_anchor &&
            is.finite(
                upper_reduction
            )
        ) {
            upper_reduction < 5
        } else {
            NA
        }

        dgi <- dg[
            dg$analysis ==
            nm,
            ,
            drop = FALSE
        ]

        landmark_row_counter <-
            landmark_row_counter +
            1L

        landmark_rows[[landmark_row_counter]] <- data.frame(
            landmark_months = L,
            analysis = nm,
            n = dgi$n,
            n_limited = dgi$n_limited,
            n_colectomy = dgi$n_colectomy,
            ESS_limited = dgi$ESS_limited,
            ESS_colectomy = dgi$ESS_colectomy,
            max_weight = dgi$max_weight,
            max_abs_wSMD =
                dgi$max_abs_wSMD,
            worst_balance_variable =
                dgi$worst_balance_variable,
            balance_pass =
                dgi$balance_pass,
            RD = pt[nm],
            RD_pp = pt[nm] * 100,
            lo = lo,
            hi = hi,
            lo_pp = lo * 100,
            hi_pp = hi * 100,
            effective_B = effB[nm],
            upper_possible_reduction_pp =
                upper_reduction,
            criterion_2pp =
                criterion_2pp,
            criterion_3pp =
                criterion_3pp,
            criterion_5pp =
                criterion_5pp,
            stringsAsFactors = FALSE
        )

        criterion_text <- ""

        if (is_decision_anchor) {

            criterion_text <- paste0(
                " | upper reduction ",
                sprintf(
                    "%.2f",
                    upper_reduction
                ),
                " pp",
                " | 2pp:",
                ifelse(
                    isTRUE(
                        criterion_2pp
                    ),
                    "MET",
                    "NOT MET"
                ),
                " 3pp:",
                ifelse(
                    isTRUE(
                        criterion_3pp
                    ),
                    "MET",
                    "NOT MET"
                ),
                " 5pp:",
                ifelse(
                    isTRUE(
                        criterion_5pp
                    ),
                    "MET",
                    "NOT MET"
                )
            )
        }

        report <- c(
            report,
            paste0(
                sprintf(
                    "%-15s RD %+.2f pp (95%% CI %+.2f to %+.2f; B=%d)",
                    nm,
                    pt[nm] * 100,
                    lo * 100,
                    hi * 100,
                    effB[nm]
                ),
                sprintf(
                    " | n=%d (%d/%d) | max|SMD|=%.3f [%s]",
                    dgi$n,
                    dgi$n_limited,
                    dgi$n_colectomy,
                    dgi$max_abs_wSMD,
                    ifelse(
                        isTRUE(
                            dgi$balance_pass
                        ),
                        "PASS",
                        "FAIL"
                    )
                ),
                criterion_text
            )
        )
    }

    report <- c(
        report,
        ""
    )
}


LANDMARK_TABLE <- do.call(
    rbind,
    landmark_rows
)

write.csv(
    LANDMARK_TABLE,
    file.path(
        RESULT_DIR,
        "PRESUB_07_Landmark_B1000.csv"
    ),
    row.names = FALSE
)


landmark_boot_wide <- do.call(
    cbind,
    lapply(
        names(
            landmark_boot_list
        ),
        function(L) {

            m <- landmark_boot_list[[L]]

            colnames(m) <- paste0(
                "L",
                L,
                "_",
                colnames(m)
            )

            m
        }
    )
)

write.csv(
    landmark_boot_wide,
    file.path(
        RESULT_DIR,
        "PRESUB_07_Landmark_Boot_B1000.csv"
    ),
    row.names = FALSE
)

report <- c(
    report,
    "Landmark interpretation:",
    "Landmark analyses were used to assess sensitivity to early postdiagnosis guarantee-time effects; because exact timing of definitive surgery was unavailable, these analyses could not eliminate immortal-time bias.",
    ""
)


## =============================================================================
## 8. PART B — Histology x risk half
## =============================================================================

hist_levels <- c(
    "nonmucinous",
    "mucinous",
    "GCA",
    "SRCC"
)

hist_labs <- c(
    "nonmucinous_low",
    "nonmucinous_high",
    "mucinous_low",
    "mucinous_high",
    "GCA_low",
    "GCA_high",
    "SRCC_low",
    "SRCC_high"
)


histology_groups <- function(d) {

    ## Preserve the prior/frozen definition:
    ## median risk is defined in the full risk-scored analytic cohort,
    ## before outcome filtering.
    median_risk <- median(
        d$nodal_risk,
        na.rm = TRUE
    )

    z <- make_outcome(
        d
    )

    masks <- list()

    for (hg in hist_levels) {

        masks[[paste0(
            hg,
            "_low"
        )]] <- (
            z$hist_group ==
            hg &
            z$nodal_risk <
            median_risk
        )

        masks[[paste0(
            hg,
            "_high"
        )]] <- (
            z$hist_group ==
            hg &
            z$nodal_risk >=
            median_risk
        )
    }

    list(
        z = z,
        median_risk = median_risk,
        masks = masks
    )
}


histology_estimate_all <- function(d) {

    obj <- histology_groups(
        d
    )

    out <- setNames(
        rep(
            NA_real_,
            length(hist_labs)
        ),
        hist_labs
    )

    for (nm in hist_labs) {

        g <- obj$z[
            obj$masks[[nm]],
            ,
            drop = FALSE
        ]

        out[nm] <- tryCatch(
            estimate_histology_rd(
                g
            ),
            error = function(e) {
                NA_real_
            }
        )
    }

    out
}


histology_point_diagnostics <- function(d) {

    obj <- histology_groups(
        d
    )

    rows <- list()

    for (nm in hist_labs) {

        g <- obj$z[
            obj$masks[[nm]],
            ,
            drop = FALSE
        ]

        dg <- tryCatch(
            diagnose_histology_subgroup(
                g
            ),
            error = function(e) {

                data.frame(
                    n = nrow(g),
                    n_limited = sum(
                        g$treat == 0,
                        na.rm = TRUE
                    ),
                    n_colectomy = sum(
                        g$treat == 1,
                        na.rm = TRUE
                    ),
                    ESS_limited = NA_real_,
                    ESS_colectomy = NA_real_,
                    max_weight = NA_real_,
                    max_abs_wSMD = NA_real_,
                    worst_balance_variable =
                        NA_character_,
                    balance_pass = FALSE,
                    stringsAsFactors = FALSE
                )
            }
        )

        dg$subgroup <- nm

        dg$histology <- sub(
            "_(low|high)$",
            "",
            nm
        )

        dg$risk_half <- ifelse(
            grepl(
                "_low$",
                nm
            ),
            "below median",
            "at/above median"
        )

        dg$overall_median_nodal_risk <-
            obj$median_risk

        rows[[nm]] <- dg
    }

    out <- do.call(
        rbind,
        rows
    )

    rownames(out) <- NULL

    out
}


pt_hist <- histology_estimate_all(
    d0
)

diag_hist <- histology_point_diagnostics(
    d0
)

boot_hist <- run_bootstrap(
    label = "Histology_RiskHalf",
    ncol_out = length(hist_labs),
    colnames_out = hist_labs,
    seed = 20260904L,
    estimator_fun = function(db) {

        dx <- pipeline(
            db
        )

        histology_estimate_all(
            dx
        )
    }
)

hist_rows <- list()

report <- c(
    report,
    "== HISTOLOGY x RISK HALF ==",
    paste0(
        "Overall median predicted nodal risk = ",
        sprintf(
            "%.4f%%",
            histology_groups(
                d0
            )$median_risk *
            100
        )
    )
)

for (
    i in seq_along(
        hist_labs
    )
) {

    nm <- hist_labs[i]

    v <- boot_hist[
        ,
        nm
    ]

    valid <- is.finite(
        v
    )

    nB <- sum(
        valid
    )

    lo <- if (
        nB >= 50
    ) {
        as.numeric(
            quantile(
                v[valid],
                0.025,
                na.rm = TRUE
            )
        )
    } else {
        NA_real_
    }

    hi <- if (
        nB >= 50
    ) {
        as.numeric(
            quantile(
                v[valid],
                0.975,
                na.rm = TRUE
            )
        )
    } else {
        NA_real_
    }

    dg <- diag_hist[
        diag_hist$subgroup ==
        nm,
        ,
        drop = FALSE
    ]

    insufficient <- (
        !is.finite(
            pt_hist[nm]
        ) ||
        dg$n < 40 ||
        dg$n_limited < 15 ||
        dg$n_colectomy < 15 ||
        nB < MIN_EFFECTIVE_B ||
        !is.finite(
            dg$max_abs_wSMD
        )
    )

    hist_rows[[i]] <- data.frame(
        subgroup = nm,
        histology = dg$histology,
        risk_half = dg$risk_half,
        overall_median_nodal_risk =
            dg$overall_median_nodal_risk,
        n = dg$n,
        n_limited = dg$n_limited,
        n_colectomy = dg$n_colectomy,
        ESS_limited = dg$ESS_limited,
        ESS_colectomy = dg$ESS_colectomy,
        max_weight = dg$max_weight,
        max_abs_wSMD =
            dg$max_abs_wSMD,
        worst_balance_variable =
            dg$worst_balance_variable,
        balance_pass =
            dg$balance_pass,
        RD = pt_hist[nm],
        RD_pp = pt_hist[nm] * 100,
        lo = lo,
        hi = hi,
        lo_pp = lo * 100,
        hi_pp = hi * 100,
        effective_B = nB,
        insufficient =
            insufficient,
        stringsAsFactors = FALSE
    )

    if (insufficient) {

        report <- c(
            report,
            sprintf(
                "%-20s Insufficient/unstable | n=%d (%d/%d) | valid B=%d | max|SMD|=%s",
                nm,
                dg$n,
                dg$n_limited,
                dg$n_colectomy,
                nB,
                ifelse(
                    is.finite(
                        dg$max_abs_wSMD
                    ),
                    sprintf(
                        "%.3f",
                        dg$max_abs_wSMD
                    ),
                    "NA"
                )
            )
        )

    } else {

        report <- c(
            report,
            sprintf(
                "%-20s RD %+.2f pp (95%% CI %+.2f to %+.2f; B=%d) | n=%d (%d/%d) | ESS %.1f/%.1f | max|SMD|=%.3f [%s]",
                nm,
                pt_hist[nm] * 100,
                lo * 100,
                hi * 100,
                nB,
                dg$n,
                dg$n_limited,
                dg$n_colectomy,
                dg$ESS_limited,
                dg$ESS_colectomy,
                dg$max_abs_wSMD,
                ifelse(
                    isTRUE(
                        dg$balance_pass
                    ),
                    "PASS",
                    "FAIL"
                )
            )
        )
    }
}


HIST_TABLE <- do.call(
    rbind,
    hist_rows
)

write.csv(
    HIST_TABLE,
    file.path(
        RESULT_DIR,
        "PRESUB_07_Histology_B1000.csv"
    ),
    row.names = FALSE
)

write.csv(
    boot_hist,
    file.path(
        RESULT_DIR,
        "PRESUB_07_Histology_Boot_B1000.csv"
    ),
    row.names = FALSE
)

report <- c(
    report,
    "",
    "Histology positioning:",
    "Histology x risk is a secondary heterogeneity analysis. It must not displace the individualized-risk primary framework, and sparse SRCC estimates remain exploratory.",
    "",
    "======================================================================",
    "FINAL PRESUBMISSION STEP 07 COMPLETE",
    "======================================================================"
)


## =============================================================================
## 9. Save report and audit
## =============================================================================

writeLines(
    report,
    file.path(
        RESULT_DIR,
        "PRESUB_07_Report.txt"
    ),
    useBytes = TRUE
)

writeLines(
    c(
        audit,
        "",
        paste0(
            "Landmark seeds: 3m=",
            LANDMARK_SEEDS["3"],
            "; 6m=",
            LANDMARK_SEEDS["6"],
            "; 12m=",
            LANDMARK_SEEDS["12"],
            "."
        ),
        "Histology seed: 20260904.",
        "Checkpoint interval: every 50 replicates.",
        "Do not rerun at a larger B because of favorable/unfavorable decision-bound results."
    ),
    file.path(
        RESULT_DIR,
        "PRESUB_07_Audit.txt"
    ),
    useBytes = TRUE
)


cat(
    "\n============================================================\n"
)

cat(
    paste(
        report,
        collapse = "\n"
    ),
    "\n"
)

cat(
    "============================================================\n"
)

cat(
    "Outputs written to: ",
    RESULT_DIR,
    "\n",
    sep = ""
)

cat(
    "Please send back PRESUB_07_Report.txt and the two summary CSV files first.\n"
)
