## =============================================================================
## A6 RED-TEAM STEP 3 — FOLLOW-UP SUPPORT / CALENDAR-ELIGIBLE SENSITIVITY
##
## PURPOSE
##   1. Quantify follow-up support for the 5-year cancer-specific mortality estimand.
##   2. Report reverse-KM follow-up.
##   3. Report unweighted numbers at risk at 12, 36, and 60 months.
##   4. Quantify how many patients have unresolved 5-year status because they
##      were censored alive before 60 months.
##   5. Re-run the corrected full pipeline in patients diagnosed <=2017,
##      a conservative calendar-eligible cohort with >=5 years of potential
##      administrative follow-up through 2023.
##
## CORRECTED ANALYTIC FRAMEWORK FOR <=2017 SENSITIVITY
##   - Summary Stage is NOT included in the propensity-score model.
##   - Nodal-risk model is re-estimated in the <=2017 population.
##   - Predicted-risk cut points are re-estimated.
##   - PS is re-fitted separately WITHIN each risk stratum.
##   - Overlap weighting is stratum-specific.
##   - B=1000 full-pipeline bootstrap.
##   - Bottom predicted-risk decile remains PRIMARY low-risk anchor.
##   - risk<5% remains SUPPORTIVE.
##
## IMPORTANT
##   - This is a follow-up/calendar-support sensitivity analysis, not a new
##     primary analysis.
##   - It does NOT overwrite any prior FINAL / REDTEAM results.
##   - Checkpointing allows the B=1000 run to resume after interruption.
##
## OUTPUTS
##   03_results/REDTEAM_03_ReverseKM.csv
##   03_results/REDTEAM_03_NumbersAtRisk.csv
##   03_results/REDTEAM_03_FiveYearSupport_ByYear.csv
##   03_results/REDTEAM_03_FiveYearSupport_ByEra.csv
##   03_results/REDTEAM_03_Le2017_Primary_B1000.csv
##   03_results/REDTEAM_03_Le2017_Diagnostics.csv
##   03_results/REDTEAM_03_Le2017_Bootstrap_B1000.csv
##   03_results/REDTEAM_03_Followup_Support_Report.txt
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
RESULT_DIR <- file.path(ROOT, "03_results")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

B <- 1000L
SEED <- 20260903L

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

CHECKPOINT_FILE <- file.path(
    RESULT_DIR,
    "REDTEAM_03_Le2017_B1000_checkpoint.rds"
)

if (is.na(PIPE_FILE)) stop("Missing Aline_frozen_pipeline_v2.RData")
if (is.na(RAW_FILE)) stop("Missing analytic_cohort_A.csv")

if (!requireNamespace("survival", quietly = TRUE)) {
    install.packages("survival", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
}

load(PIPE_FILE)

if (!exists("pipeline", mode = "function")) stop("pipeline() missing")
if (!exists("PS_FORMULA")) stop("PS_FORMULA missing")

RAW <- read.csv(
    RAW_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(RAW) != 6951) {
    stop("Frozen cohort fingerprint failed: N must be 6951")
}

main <- function() {

    ## =========================================================================
    ## 1. Correct PS formula — remove Summary Stage
    ## =========================================================================

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
        stop("No stage-related term found in PS_FORMULA")
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

    had_env_formula <- exists(
        "PS_FORMULA",
        envir = pipeline_env,
        inherits = FALSE
    )

    old_env_formula <- if (had_env_formula) {
        get(
            "PS_FORMULA",
            envir = pipeline_env,
            inherits = FALSE
        )
    } else {
        NULL
    }

    restore_formula <- function() {
        assign(
            "PS_FORMULA",
            PS_ORIG,
            envir = .GlobalEnv
        )

        if (had_env_formula) {
            assign(
                "PS_FORMULA",
                old_env_formula,
                envir = pipeline_env
            )
        }
    }

    on.exit(
        restore_formula(),
        add = TRUE
    )

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

    ## =========================================================================
    ## 2. Common helpers
    ## =========================================================================

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

        m <- time > start

        time <- time[m]
        ev1 <- ev1[m]
        ev2 <- ev2[m]
        wt <- wt[m]

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

                S <- S * (
                    1 -
                    (dc[k] + dr[k]) / ar
                )
            }

            ar <- ar - dw[k]
        }

        as.numeric(cif)
    }

    make_outcome <- function(d) {
        z <- d[
            !is.na(d$css_event) &
            !is.na(d$os_months) &
            !is.na(d$os_event) &
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
        cuts <- make_cuts(
            d$nodal_risk
        )

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
        if (nrow(g) < 40) {
            stop("Subgroup too small")
        }

        if (sum(g$treat == 0) < 15) {
            stop("Too few limited-resection patients")
        }

        if (sum(g$treat == 1) < 15) {
            stop("Too few colectomy patients")
        }

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

        if (any(!is.finite(ps))) {
            stop("Non-finite subgroup propensity score")
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
            !is.finite(sum(g$ow_sub)) ||
            sum(g$ow_sub) <= 0
        ) {
            stop("Invalid subgroup overlap weights")
        }

        g
    }

    estimate_one <- function(g) {
        gs <- fit_subgroup_ow(g)

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
            a1$ow_sub
        )

        c0 <- aj_w(
            a0$os_months,
            a0$cancer,
            a0$otherd,
            a0$ow_sub
        )

        c(
            limited = c0,
            colectomy = c1,
            rd = c1 - c0
        )
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

    estimate_all <- function(d) {
        z <- make_outcome(d)
        masks <- make_masks(z)

        out <- rep(
            NA_real_,
            length(labs)
        )

        names(out) <- labs

        for (i in seq_along(labs)) {
            g <- z[
                masks[[labs[i]]],
                ,
                drop = FALSE
            ]

            out[i] <- tryCatch(
                estimate_one(g)["rd"],
                error = function(e) NA_real_
            )
        }

        out
    }

    ess <- function(w) {
        w <- w[
            is.finite(w) &
            w > 0
        ]

        if (length(w) == 0) {
            return(NA_real_)
        }

        sum(w)^2 / sum(w^2)
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
            sum(tr == 1) == 0 ||
            sum(tr == 0) == 0
        ) {
            return(NA_real_)
        }

        m1 <- weighted.mean(
            x[tr == 1],
            w[tr == 1]
        )

        m0 <- weighted.mean(
            x[tr == 0],
            w[tr == 0]
        )

        v1 <- weighted.mean(
            (x[tr == 1] - m1)^2,
            w[tr == 1]
        )

        v0 <- weighted.mean(
            (x[tr == 0] - m0)^2,
            w[tr == 0]
        )

        denom <- sqrt(
            (v1 + v0) / 2
        )

        if (
            !is.finite(denom) ||
            denom <= 0
        ) {
            return(0)
        }

        (m1 - m0) / denom
    }

    balance_design <- function(g) {
        X <- model.matrix(
            delete.response(
                terms(PS_NOSTAGE)
            ),
            data = g,
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

        cbind(
            X,
            nodal_risk = g$nodal_risk
        )
    }

    ## =========================================================================
    ## 3. Full-cohort outcome dataset for follow-up support
    ## =========================================================================

    d_full <- pipeline(RAW)
    z_full <- make_outcome(d_full)

    outcome_excluded_n <- (
        nrow(RAW) -
        nrow(z_full)
    )

    z_full$exposure_label <- ifelse(
        z_full$treat == 1,
        "Oncologic colectomy",
        "Limited resection"
    )

    z_full$era_follow <- ifelse(
        z_full$year <= 2009,
        "2004-2009",
        ifelse(
            z_full$year <= 2017,
            "2010-2017",
            "2018-2023"
        )
    )

    ## =========================================================================
    ## 4. Reverse-KM follow-up
    ## =========================================================================

    reverse_km_one <- function(
        d,
        label
    ) {
        if (nrow(d) < 5) {
            return(
                data.frame(
                    group = label,
                    n = nrow(d),
                    median_months = NA_real_,
                    lower95_months = NA_real_,
                    upper95_months = NA_real_,
                    stringsAsFactors = FALSE
                )
            )
        }

        fit <- survival::survfit(
            survival::Surv(
                d$os_months,
                1 - d$os_event
            ) ~ 1
        )

        tab <- summary(fit)$table

        get_stat <- function(nm) {
            if (nm %in% names(tab)) {
                as.numeric(tab[nm])
            } else {
                NA_real_
            }
        }

        data.frame(
            group = label,
            n = nrow(d),
            median_months = get_stat("median"),
            lower95_months = get_stat("0.95LCL"),
            upper95_months = get_stat("0.95UCL"),
            stringsAsFactors = FALSE
        )
    }

    rk_list <- list()

    rk_list[[length(rk_list) + 1]] <- reverse_km_one(
        z_full,
        "Overall"
    )

    for (x in c(
        "Limited resection",
        "Oncologic colectomy"
    )) {
        rk_list[[length(rk_list) + 1]] <- reverse_km_one(
            z_full[
                z_full$exposure_label == x,
                ,
                drop = FALSE
            ],
            x
        )
    }

    for (x in c(
        "2004-2009",
        "2010-2017",
        "2018-2023"
    )) {
        rk_list[[length(rk_list) + 1]] <- reverse_km_one(
            z_full[
                z_full$era_follow == x,
                ,
                drop = FALSE
            ],
            paste0("Era ", x)
        )
    }

    rk_list[[length(rk_list) + 1]] <- reverse_km_one(
        z_full[
            z_full$year <= 2017,
            ,
            drop = FALSE
        ],
        "Diagnosis <=2017"
    )

    rk_list[[length(rk_list) + 1]] <- reverse_km_one(
        z_full[
            z_full$year >= 2018,
            ,
            drop = FALSE
        ],
        "Diagnosis 2018-2023"
    )

    REVERSE_KM <- do.call(
        rbind,
        rk_list
    )

    ## =========================================================================
    ## 5. Numbers at risk at 12 / 36 / 60 months
    ## =========================================================================

    risk_times <- c(
        12,
        36,
        60
    )

    numbers_at_risk_one <- function(
        d,
        label
    ) {
        fit <- survival::survfit(
            survival::Surv(
                d$os_months,
                d$os_event
            ) ~ 1
        )

        ss <- summary(
            fit,
            times = risk_times,
            extend = TRUE
        )

        data.frame(
            group = label,
            time_months = ss$time,
            n_risk = ss$n.risk,
            stringsAsFactors = FALSE
        )
    }

    masks_full <- make_masks(
        z_full
    )

    nar_list <- list()

    nar_list[[length(nar_list) + 1]] <- numbers_at_risk_one(
        z_full,
        "Overall"
    )

    for (x in c(
        "Limited resection",
        "Oncologic colectomy"
    )) {
        nar_list[[length(nar_list) + 1]] <- numbers_at_risk_one(
            z_full[
                z_full$exposure_label == x,
                ,
                drop = FALSE
            ],
            x
        )
    }

    nar_list[[length(nar_list) + 1]] <- numbers_at_risk_one(
        z_full[
            masks_full$bottom_decile,
            ,
            drop = FALSE
        ],
        "Bottom predicted-risk decile"
    )

    nar_list[[length(nar_list) + 1]] <- numbers_at_risk_one(
        z_full[
            masks_full$risk_lt_5pct,
            ,
            drop = FALSE
        ],
        "Predicted nodal risk <5%"
    )

    NUMBERS_AT_RISK <- do.call(
        rbind,
        nar_list
    )

    ## =========================================================================
    ## 6. Direct 5-year support / unresolved censoring
    ## =========================================================================

    support_summary <- function(
        d,
        group_name
    ) {
        unresolved <- (
            d$os_event == 0 &
            d$os_months < 60
        )

        death_before60 <- (
            d$os_event == 1 &
            d$os_months < 60
        )

        followed60 <- (
            d$os_months >= 60
        )

        data.frame(
            group = group_name,
            n = nrow(d),
            death_before_60_n = sum(
                death_before60,
                na.rm = TRUE
            ),
            followed_at_least_60_n = sum(
                followed60,
                na.rm = TRUE
            ),
            censored_alive_before_60_n = sum(
                unresolved,
                na.rm = TRUE
            ),
            five_year_status_resolved_n = sum(
                !unresolved,
                na.rm = TRUE
            ),
            five_year_status_resolved_pct = 100 * mean(
                !unresolved,
                na.rm = TRUE
            ),
            stringsAsFactors = FALSE
        )
    }

    SUPPORT_BY_YEAR <- do.call(
        rbind,
        lapply(
            sort(
                unique(
                    z_full$year
                )
            ),
            function(y) {
                out <- support_summary(
                    z_full[
                        z_full$year == y,
                        ,
                        drop = FALSE
                    ],
                    as.character(y)
                )

                out$diagnosis_year <- y
                out$calendar_full_5yr_potential <- (
                    y <= 2017
                )

                out
            }
        )
    )

    SUPPORT_BY_YEAR <- SUPPORT_BY_YEAR[
        ,
        c(
            "diagnosis_year",
            "calendar_full_5yr_potential",
            setdiff(
                names(SUPPORT_BY_YEAR),
                c(
                    "diagnosis_year",
                    "calendar_full_5yr_potential"
                )
            )
        ),
        drop = FALSE
    ]

    era_levels <- c(
        "2004-2009",
        "2010-2017",
        "2018-2023"
    )

    SUPPORT_BY_ERA <- do.call(
        rbind,
        lapply(
            era_levels,
            function(er) {
                support_summary(
                    z_full[
                        z_full$era_follow == er,
                        ,
                        drop = FALSE
                    ],
                    er
                )
            }
        )
    )

    SUPPORT_OVERALL <- support_summary(
        z_full,
        "Overall"
    )

    ## =========================================================================
    ## 7. <=2017 point estimates and balance diagnostics
    ## =========================================================================

    RAW17 <- RAW[
        RAW$year <= 2017,
        ,
        drop = FALSE
    ]

    d17 <- pipeline(
        RAW17
    )

    z17 <- make_outcome(
        d17
    )

    masks17 <- make_masks(
        z17
    )

    point17 <- rep(
        NA_real_,
        length(labs)
    )

    names(point17) <- labs

    diag_list <- list()

    for (nm in labs) {
        g <- z17[
            masks17[[nm]],
            ,
            drop = FALSE
        ]

        gs <- fit_subgroup_ow(
            g
        )

        est <- estimate_one(
            g
        )

        point17[nm] <- est["rd"]

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

        Xbal <- balance_design(
            gs
        )

        smd_vec <- vapply(
            seq_len(
                ncol(Xbal)
            ),
            function(j) {
                weighted_smd(
                    Xbal[, j],
                    gs$treat,
                    gs$ow_sub
                )
            },
            numeric(1)
        )

        names(smd_vec) <- colnames(
            Xbal
        )

        diag_list[[nm]] <- data.frame(
            analysis = nm,
            n = nrow(gs),
            limited = nrow(a0),
            colectomy = nrow(a1),
            cancer_events_limited = sum(
                a0$cancer == 1,
                na.rm = TRUE
            ),
            cancer_events_colectomy = sum(
                a1$cancer == 1,
                na.rm = TRUE
            ),
            cancer_CIF_limited_pct = est["limited"] * 100,
            cancer_CIF_colectomy_pct = est["colectomy"] * 100,
            cancer_RD_pp = est["rd"] * 100,
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
            max_abs_weighted_SMD = max(
                abs(smd_vec),
                na.rm = TRUE
            ),
            worst_balance_variable = names(
                which.max(
                    abs(smd_vec)
                )
            ),
            balance_pass = (
                max(
                    abs(smd_vec),
                    na.rm = TRUE
                ) < 0.10
            ),
            stringsAsFactors = FALSE
        )
    }

    DIAGNOSTICS17 <- do.call(
        rbind,
        diag_list
    )

    ## =========================================================================
    ## 8. <=2017 B=1000 full-pipeline bootstrap with checkpoint/resume
    ## =========================================================================

    boot <- matrix(
        NA_real_,
        nrow = B,
        ncol = length(labs),
        dimnames = list(
            NULL,
            labs
        )
    )

    start_b <- 1L

    if (file.exists(CHECKPOINT_FILE)) {
        cp <- tryCatch(
            readRDS(
                CHECKPOINT_FILE
            ),
            error = function(e) NULL
        )

        if (
            !is.null(cp) &&
            is.matrix(cp$boot) &&
            identical(
                dim(cp$boot),
                dim(boot)
            ) &&
            identical(
                colnames(cp$boot),
                colnames(boot)
            )
        ) {
            boot <- cp$boot

            done <- which(
                apply(
                    boot,
                    1,
                    function(x) {
                        any(
                            is.finite(x)
                        )
                    }
                )
            )

            start_b <- if (
                length(done) == 0
            ) {
                1L
            } else {
                max(done) + 1L
            }

            if (
                !is.null(cp$rng_seed)
            ) {
                assign(
                    ".Random.seed",
                    cp$rng_seed,
                    envir = .GlobalEnv
                )
            }

            cat(
                "Checkpoint loaded. Resume at bootstrap ",
                start_b,
                "/",
                B,
                "\n",
                sep = ""
            )
        }
    }

    if (start_b == 1L) {
        set.seed(
            SEED
        )
    }

    if (start_b <= B) {
        t0 <- Sys.time()

        for (b in seq.int(
            start_b,
            B
        )) {
            db <- RAW17[
                sample(
                    seq_len(
                        nrow(RAW17)
                    ),
                    size = nrow(RAW17),
                    replace = TRUE
                ),
                ,
                drop = FALSE
            ]

            boot[b, ] <- tryCatch(
                {
                    dx <- pipeline(
                        db
                    )

                    estimate_all(
                        dx
                    )
                },
                error = function(e) {
                    rep(
                        NA_real_,
                        length(labs)
                    )
                }
            )

            if (
                b %% 50 == 0 ||
                b == B
            ) {
                saveRDS(
                    list(
                        boot = boot,
                        completed_through = b,
                        rng_seed = .Random.seed
                    ),
                    CHECKPOINT_FILE
                )

                cat(
                    "<=2017 bootstrap ",
                    b,
                    "/",
                    B,
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
    }

    ## =========================================================================
    ## 9. Final <=2017 bootstrap table
    ## =========================================================================

    result_list <- list()

    for (i in seq_along(labs)) {
        v <- boot[
            ,
            i
        ]

        valid <- is.finite(
            v
        )

        if (sum(valid) < 800) {
            warning(
                labs[i],
                ": effective B < 800 (",
                sum(valid),
                ")"
            )
        }

        lo <- quantile(
            v[valid],
            0.025,
            na.rm = TRUE
        )

        hi <- quantile(
            v[valid],
            0.975,
            na.rm = TRUE
        )

        result_list[[i]] <- data.frame(
            analysis = labs[i],
            RD = point17[i],
            RD_pp = point17[i] * 100,
            lo = lo,
            hi = hi,
            lo_pp = lo * 100,
            hi_pp = hi * 100,
            upper_possible_reduction_pp = max(
                0,
                -lo * 100
            ),
            effective_B = sum(
                valid
            ),
            stringsAsFactors = FALSE
        )
    }

    RESULTS17 <- do.call(
        rbind,
        result_list
    )

    lowrisk <- c(
        "Q1",
        "bottom_decile",
        "risk_lt_5pct"
    )

    RESULTS17$criterion_2pp <- ifelse(
        RESULTS17$analysis %in% lowrisk,
        RESULTS17$upper_possible_reduction_pp < 2,
        NA
    )

    RESULTS17$criterion_3pp <- ifelse(
        RESULTS17$analysis %in% lowrisk,
        RESULTS17$upper_possible_reduction_pp < 3,
        NA
    )

    RESULTS17$criterion_5pp <- ifelse(
        RESULTS17$analysis %in% lowrisk,
        RESULTS17$upper_possible_reduction_pp < 5,
        NA
    )

    ## =========================================================================
    ## 10. Compare <=2017 to corrected Step 1B full-cohort result
    ## =========================================================================

    if (file.exists(STEP1B_FILE)) {
        full_ref <- read.csv(
            STEP1B_FILE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )

        rd_col <- if (
            "RD_pp" %in% names(full_ref)
        ) {
            "RD_pp"
        } else if (
            "rd_pp" %in% names(full_ref)
        ) {
            "rd_pp"
        } else {
            NA_character_
        }

        if (!is.na(rd_col)) {
            ref <- data.frame(
                analysis = full_ref$analysis,
                full_cohort_RD_pp = full_ref[[rd_col]],
                stringsAsFactors = FALSE
            )

            RESULTS17 <- merge(
                RESULTS17,
                ref,
                by = "analysis",
                all.x = TRUE,
                sort = FALSE
            )

            RESULTS17 <- RESULTS17[
                match(
                    labs,
                    RESULTS17$analysis
                ),
                ,
                drop = FALSE
            ]

            RESULTS17$point_change_vs_full_pp <- (
                RESULTS17$RD_pp -
                RESULTS17$full_cohort_RD_pp
            )
        }
    }

    ## =========================================================================
    ## 11. Save machine-readable files
    ## =========================================================================

    OUT_RK <- file.path(
        RESULT_DIR,
        "REDTEAM_03_ReverseKM.csv"
    )

    OUT_NAR <- file.path(
        RESULT_DIR,
        "REDTEAM_03_NumbersAtRisk.csv"
    )

    OUT_YEAR <- file.path(
        RESULT_DIR,
        "REDTEAM_03_FiveYearSupport_ByYear.csv"
    )

    OUT_ERA <- file.path(
        RESULT_DIR,
        "REDTEAM_03_FiveYearSupport_ByEra.csv"
    )

    OUT_RES17 <- file.path(
        RESULT_DIR,
        "REDTEAM_03_Le2017_Primary_B1000.csv"
    )

    OUT_DIAG17 <- file.path(
        RESULT_DIR,
        "REDTEAM_03_Le2017_Diagnostics.csv"
    )

    OUT_BOOT17 <- file.path(
        RESULT_DIR,
        "REDTEAM_03_Le2017_Bootstrap_B1000.csv"
    )

    OUT_REPORT <- file.path(
        RESULT_DIR,
        "REDTEAM_03_Followup_Support_Report.txt"
    )

    write.csv(
        REVERSE_KM,
        OUT_RK,
        row.names = FALSE
    )

    write.csv(
        NUMBERS_AT_RISK,
        OUT_NAR,
        row.names = FALSE
    )

    write.csv(
        SUPPORT_BY_YEAR,
        OUT_YEAR,
        row.names = FALSE
    )

    write.csv(
        SUPPORT_BY_ERA,
        OUT_ERA,
        row.names = FALSE
    )

    write.csv(
        RESULTS17,
        OUT_RES17,
        row.names = FALSE
    )

    write.csv(
        DIAGNOSTICS17,
        OUT_DIAG17,
        row.names = FALSE
    )

    write.csv(
        data.frame(
            replicate = seq_len(B),
            boot,
            check.names = FALSE
        ),
        OUT_BOOT17,
        row.names = FALSE
    )

    ## =========================================================================
    ## 12. Human-readable report
    ## =========================================================================

    fmt_num <- function(x, digits = 1) {
        ifelse(
            is.finite(x),
            formatC(
                x,
                format = "f",
                digits = digits
            ),
            "NA"
        )
    }

    report <- c(
        "======================================================================",
        "A6 RED-TEAM STEP 3 — FOLLOW-UP SUPPORT / <=2017 SENSITIVITY",
        "======================================================================",
        paste0(
            "Frozen cohort N = ",
            nrow(RAW)
        ),
        paste0(
            "Outcome-evaluable N = ",
            nrow(z_full),
            " | excluded from outcome analysis = ",
            outcome_excluded_n
        ),
        paste0(
            "Removed PS term(s): ",
            paste(
                stage_terms,
                collapse = " | "
            )
        ),
        "Primary low-risk anchor remains bottom predicted-risk decile.",
        "risk<5% remains supportive.",
        "",
        "== REVERSE-KM FOLLOW-UP =="
    )

    for (i in seq_len(
        nrow(REVERSE_KM)
    )) {
        z <- REVERSE_KM[
            i,
            ,
            drop = FALSE
        ]

        report <- c(
            report,
            sprintf(
                "%-28s n=%d | median %.1f months (95%% CI %.1f to %.1f)",
                z$group,
                z$n,
                z$median_months,
                z$lower95_months,
                z$upper95_months
            )
        )
    }

    report <- c(
        report,
        "",
        "== UNWEIGHTED NUMBERS AT RISK =="
    )

    groups_order <- unique(
        NUMBERS_AT_RISK$group
    )

    for (g in groups_order) {
        z <- NUMBERS_AT_RISK[
            NUMBERS_AT_RISK$group == g,
            ,
            drop = FALSE
        ]

        vals <- paste(
            paste0(
                z$time_months,
                "m=",
                z$n_risk
            ),
            collapse = " | "
        )

        report <- c(
            report,
            paste0(
                g,
                " | ",
                vals
            )
        )
    }

    report <- c(
        report,
        "",
        "== 5-YEAR STATUS SUPPORT =="
    )

    report <- c(
        report,
        sprintf(
            "Overall: n=%d | death before 60m=%d | followed >=60m=%d | censored alive <60m=%d | 5y status resolved %.1f%%",
            SUPPORT_OVERALL$n,
            SUPPORT_OVERALL$death_before_60_n,
            SUPPORT_OVERALL$followed_at_least_60_n,
            SUPPORT_OVERALL$censored_alive_before_60_n,
            SUPPORT_OVERALL$five_year_status_resolved_pct
        )
    )

    for (i in seq_len(
        nrow(SUPPORT_BY_ERA)
    )) {
        z <- SUPPORT_BY_ERA[
            i,
            ,
            drop = FALSE
        ]

        report <- c(
            report,
            sprintf(
                "%-10s n=%d | censored alive <60m=%d | 5y status resolved %.1f%%",
                z$group,
                z$n,
                z$censored_alive_before_60_n,
                z$five_year_status_resolved_pct
            )
        )
    }

    report <- c(
        report,
        "",
        "Calendar-eligible sensitivity uses diagnosis year <=2017, conservatively ensuring at least 5 years of potential administrative follow-up through 2023.",
        paste0(
            "<=2017 raw N = ",
            nrow(RAW17),
            " | outcome-evaluable N = ",
            nrow(z17)
        ),
        "",
        "== <=2017 CORRECTED FULL-PIPELINE B=1000 =="
    )

    for (i in seq_len(
        nrow(RESULTS17)
    )) {
        z <- RESULTS17[
            i,
            ,
            drop = FALSE
        ]

        extra <- ""

        if (
            z$analysis %in%
            lowrisk
        ) {
            extra <- sprintf(
                " | upper reduction %.1f pp | 2pp:%s 3pp:%s 5pp:%s",
                z$upper_possible_reduction_pp,
                ifelse(
                    z$criterion_2pp,
                    "MET",
                    "NOT MET"
                ),
                ifelse(
                    z$criterion_3pp,
                    "MET",
                    "NOT MET"
                ),
                ifelse(
                    z$criterion_5pp,
                    "MET",
                    "NOT MET"
                )
            )
        }

        change_text <- if (
            "point_change_vs_full_pp" %in%
            names(RESULTS17)
        ) {
            sprintf(
                " | change vs full %+.2f pp",
                z$point_change_vs_full_pp
            )
        } else {
            ""
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
                    z$effective_B
                ),
                extra,
                change_text
            )
        )
    }

    report <- c(
        report,
        "",
        "== <=2017 BALANCE / SUPPORT =="
    )

    for (i in seq_len(
        nrow(DIAGNOSTICS17)
    )) {
        z <- DIAGNOSTICS17[
            i,
            ,
            drop = FALSE
        ]

        report <- c(
            report,
            sprintf(
                paste0(
                    "%-15s n=%d (%d/%d) events=%d/%d | ",
                    "CIF %.1f/%.1f%% RD %+.1f | ESS %.1f/%.1f | ",
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

    bottom <- RESULTS17[
        RESULTS17$analysis ==
        "bottom_decile",
        ,
        drop = FALSE
    ]

    risk5 <- RESULTS17[
        RESULTS17$analysis ==
        "risk_lt_5pct",
        ,
        drop = FALSE
    ]

    report <- c(
        report,
        "",
        "== DECISION / SCORE CHECK =="
    )

    if (
        nrow(bottom) == 1 &&
        is.finite(
            bottom$RD_pp
        )
    ) {
        if (
            "point_change_vs_full_pp" %in%
            names(bottom) &&
            is.finite(
                bottom$point_change_vs_full_pp
            ) &&
            abs(
                bottom$point_change_vs_full_pp
            ) > 2
        ) {
            report <- c(
                report,
                paste0(
                    "FOLLOW-UP SCORE ALERT: bottom-decile point estimate moved ",
                    sprintf(
                        "%+.2f",
                        bottom$point_change_vs_full_pp
                    ),
                    " pp vs the full corrected cohort (>2 pp)."
                )
            )
        } else {
            report <- c(
                report,
                "No >2-pp bottom-decile point-estimate shift versus the full corrected cohort."
            )
        }

        report <- c(
            report,
            paste0(
                "<=2017 bottom-decile 3-pp criterion: ",
                ifelse(
                    isTRUE(
                        bottom$criterion_3pp
                    ),
                    "MET",
                    "NOT MET"
                ),
                " | upper compatible reduction ",
                sprintf(
                    "%.2f",
                    bottom$upper_possible_reduction_pp
                ),
                " pp."
            )
        )
    }

    if (nrow(risk5) == 1) {
        report <- c(
            report,
            paste0(
                "<=2017 risk<5% 3-pp criterion: ",
                ifelse(
                    isTRUE(
                        risk5$criterion_3pp
                    ),
                    "MET",
                    "NOT MET"
                ),
                " | upper compatible reduction ",
                sprintf(
                    "%.2f",
                    risk5$upper_possible_reduction_pp
                ),
                " pp."
            )
        )
    }

    report <- c(
        report,
        "",
        "Interpretation rule:",
        "The <=2017 analysis is a calendar-eligibility / follow-up-support sensitivity analysis. Wider confidence intervals caused by the smaller historical cohort are not, by themselves, evidence that the main result is invalid.",
        "The key robustness question is whether the low-risk point estimate materially changes direction/magnitude when late-diagnosis patients are removed, while within-stratum balance remains acceptable.",
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
}

main()
