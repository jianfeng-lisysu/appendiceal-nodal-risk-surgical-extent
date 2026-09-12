## =============================================================================
## A6 - STEP 06 S2
## Alternative surgery-definition sensitivity analyses from ORIGINAL SEER source
## FINAL full-pipeline bootstrap B=1000 - ANALYTIC-CONTENT GATE VERSION
##
## Root:
##   <project root>
##
## Inputs:
##   02_data/SEER_appendix_malignant_2000_2023.txt
##   02_data/analytic_cohort_A.csv
##   02_data/Aline_frozen_pipeline_v2.RData
##
## Outputs:
##   03_results/06_S2_alternative_surgery_FINAL_B1000.txt
##   03_results/06_S2_alternative_surgery_FINAL_B1000.csv
##   03_results/06_S2_preE7_audit_FINAL.txt
##   03_results/06_S2_bootstrap_FINAL_B1000.csv
##
## S2 definitions reproduced from the prior frozen sensitivity analysis:
##
##   S2-A:
##     limited resection = surgery codes 20-39
##     oncologic colectomy = surgery codes 40-79
##     expected N = 7,450
##
##   S2-B:
##     main definition but exclude surgery code 32
##     limited resection = codes 30-31
##     oncologic colectomy = codes 40-41
##     expected N = 6,676
##
##   S2-C:
##     broaden colectomy to codes 50-79
##     limited resection = codes 30-32
##     oncologic colectomy = codes 40-79
##     expected N = 7,115
##
## IMPORTANT HARMONIZATION REPAIRS:
##
##   1. T stage includes:
##      Derived SEER Combined T (2016-2017)
##
##   2. 2023 surgery codes are alphanumeric:
##      A300 -> 30
##      A320 -> 32
##      A400 -> 40
##      A410 -> 41
##      etc.
##
## The script MUST reproduce the frozen main cohort N=6,951 and its main
## point estimates before S2 is allowed to run.
## =============================================================================


## =============================================================================
## 0. Setup
## =============================================================================

rm(list = ls())

options(
    stringsAsFactors = FALSE
)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
RESULT_DIR <- file.path(ROOT, "03_results")

dir.create(
    RESULT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

SEER_FILE <- file.path(
    DATA_DIR,
    "SEER_appendix_malignant_2000_2023.txt"
)

FROZEN_MAIN_FILE <- file.path(
    DATA_DIR,
    "analytic_cohort_A.csv"
)

PIPE_FILE <- file.path(
    DATA_DIR,
    "Aline_frozen_pipeline_v2.RData"
)

## Backward-compatible fallback
if (!file.exists(SEER_FILE)) {
    SEER_FILE <- file.path(
        ROOT,
        "SEER_appendix_malignant_2000_2023.txt"
    )
}

if (!file.exists(FROZEN_MAIN_FILE)) {
    FROZEN_MAIN_FILE <- file.path(
        ROOT,
        "analytic_cohort_A.csv"
    )
}

if (!file.exists(PIPE_FILE)) {
    PIPE_FILE <- file.path(
        ROOT,
        "Aline_frozen_pipeline_v2.RData"
    )
}

if (!file.exists(SEER_FILE)) {
    stop("Cannot find original SEER source file.")
}

if (!file.exists(FROZEN_MAIN_FILE)) {
    stop("Cannot find frozen analytic_cohort_A.csv.")
}

if (!file.exists(PIPE_FILE)) {
    stop("Cannot find frozen pipeline RData.")
}

load(
    PIPE_FILE
)

needed_objects <- c(
    "prep",
    "pipeline",
    "RISK_FORMULA",
    "PS_FORMULA"
)

missing_objects <- needed_objects[
    !vapply(
        needed_objects,
        exists,
        logical(1),
        inherits = TRUE
    )
]

if (length(missing_objects) > 0) {

    stop(
        "Frozen pipeline is missing: ",
        paste(
            missing_objects,
            collapse = ", "
        )
    )
}


## =============================================================================
## 1. Read original SEER source
## =============================================================================

cat("============================================================\n")
cat("STEP 06 S2 - reconstructing pre-E7 cohort\n")
cat("============================================================\n\n")

df <- read.delim(
    SEER_FILE,
    sep = "\t",
    quote = "\"",
    colClasses = "character",
    check.names = FALSE,
    na.strings = NULL,
    fileEncoding = "UTF-8"
)

names(df) <- trimws(
    names(df)
)


## =============================================================================
## 2. Exact frozen cohort construction through E6
## =============================================================================

FLOW <- data.frame(
    step = character(0),
    n = integer(0),
    stringsAsFactors = FALSE
)

add_flow <- function(
    step,
    n
) {

    FLOW <<- rbind(
        FLOW,
        data.frame(
            step = step,
            n = as.integer(n),
            stringsAsFactors = FALSE
        )
    )
}


add_flow(
    "E0",
    nrow(df)
)


## -----------------------------------------------------------------------------
## E1: diagnosis years 2004-2023
## -----------------------------------------------------------------------------

df$year <- suppressWarnings(
    as.numeric(
        df[["Year of diagnosis"]]
    )
)

df <- df[
    !is.na(df$year) &
    df$year >= 2004 &
    df$year <= 2023,
    ,
    drop = FALSE
]

add_flow(
    "E1",
    nrow(df)
)


## -----------------------------------------------------------------------------
## E2: prespecified appendiceal histology groups
## -----------------------------------------------------------------------------

H <- "Histologic Type ICD-O-3"

GROUPS <- list(
    nonmucinous = c(
        "8140",
        "8263",
        "8261",
        "8262",
        "8255",
        "8010",
        "8144",
        "8210",
        "8211"
    ),
    mucinous = c(
        "8480",
        "8481",
        "8470"
    ),
    GCA = c(
        "8243",
        "8244",
        "8245",
        "8249"
    ),
    SRCC = c(
        "8490"
    )
)

df$hist_group <- NA_character_

for (g in names(GROUPS)) {

    hist_codes <- getElement(
        GROUPS,
        g
    )

    df$hist_group[
        df[, H] %in%
        hist_codes
    ] <- g
}

df <- df[
    !is.na(
        df$hist_group
    ),
    ,
    drop = FALSE
]

add_flow(
    "E2",
    nrow(df)
)

cat(
    "E2 histology counts: ",
    paste(
        names(
            table(
                df$hist_group
            )
        ),
        as.integer(
            table(
                df$hist_group
            )
        ),
        sep = "=",
        collapse = " | "
    ),
    "\n",
    sep = ""
)


## -----------------------------------------------------------------------------
## E3: microscopic confirmation
## -----------------------------------------------------------------------------

OKDX <- c(
    "Positive histology",
    "Pos hist AND immunophenotyping AND/OR pos genetic studies",
    "Positive microscopic confirm, method not specified"
)

df <- df[
    df[["Diagnostic Confirmation"]] %in% OKDX,
    ,
    drop = FALSE
]

add_flow(
    "E3",
    nrow(df)
)


## -----------------------------------------------------------------------------
## E4: exclude death certificate/autopsy only
## -----------------------------------------------------------------------------

df <- df[
    !(
        df[["Type of Reporting Source"]] %in%
        c(
            "Death certificate only",
            "Autopsy only"
        )
    ),
    ,
    drop = FALSE
]

add_flow(
    "E4",
    nrow(df)
)


## -----------------------------------------------------------------------------
## E5: localized/regional summary stage
## -----------------------------------------------------------------------------

CSS <- "Combined Summary Stage with Expanded Regional Codes (2004+)"

df <- df[
    grepl(
        "^(Localized|Regional)",
        df[, CSS]
    ),
    ,
    drop = FALSE
]

add_flow(
    "E5",
    nrow(df)
)


## -----------------------------------------------------------------------------
## E6: harmonized T stage
##
## IMPORTANT:
## 2016-2017 uses Derived SEER Combined T, whose values are commonly p1A/p2/p3...
## rather than literal T1/T2/T3...
## -----------------------------------------------------------------------------

tnorm <- function(v) {

    s <- gsub(
        " ",
        "",
        toupper(
            as.character(v)
        ),
        fixed = TRUE
    )

    out <- rep(
        NA_character_,
        length(s)
    )

    for (k in 1:4) {

        mask <- (
            is.na(out) &
            grepl(
                paste0(
                    "^[TCP]?",
                    k
                ),
                s
            )
        )

        out[
            mask
        ] <- paste0(
            "T",
            k
        )
    }

    out
}


df$T_stage <- NA_character_

T_COLS <- c(
    "Derived EOD 2018 T Recode (2018+)",
    "Derived SEER Combined T (2016-2017)",
    "Derived AJCC T, 7th ed (2010-2015)",
    "Derived AJCC T, 6th ed (2004-2015)"
)

missing_t_cols <- T_COLS[
    !T_COLS %in%
    names(df)
]

if (length(missing_t_cols) > 0) {

    stop(
        "Missing T-stage field(s): ",
        paste(
            missing_t_cols,
            collapse = " | "
        )
    )
}

for (col in T_COLS) {

    z <- tnorm(
        df[, col]
    )

    mask <- is.na(
        df$T_stage
    )

    df$T_stage[
        mask
    ] <- z[
        mask
    ]
}

df <- df[
    !is.na(
        df$T_stage
    ),
    ,
    drop = FALSE
]

add_flow(
    "E6",
    nrow(df)
)


## =============================================================================
## 3. Exact 2023 surgery-code harmonization
## =============================================================================

sc_old <- suppressWarnings(
    as.numeric(
        df[["RX Summ--Surg Prim Site (1998-2022)"]]
    )
)

s23_raw <- as.character(
    df[["RX Summ--Surg Prim Site 2023 (2023+)"]]
)

s23_num <- suppressWarnings(
    as.numeric(
        sub(
            "^A",
            "",
            s23_raw
        )
    )
)

valid_2023_code <- grepl(
    "^A[0-9]{3}$",
    s23_raw
)

sc23 <- rep(
    NA_real_,
    length(
        s23_raw
    )
)

sc23[
    valid_2023_code
] <- s23_num[
    valid_2023_code
] / 10

sc <- sc_old

replace_2023 <- (
    is.na(sc) &
    !is.na(sc23)
)

sc[
    replace_2023
] <- sc23[
    replace_2023
]

df$surg_code <- sc


## =============================================================================
## 4. Create all non-exposure analytic variables BEFORE E7
## =============================================================================

## -----------------------------------------------------------------------------
## Survival / E8
## -----------------------------------------------------------------------------

df$os_months <- suppressWarnings(
    as.numeric(
        df[["Survival months"]]
    )
)

df <- df[
    !is.na(
        df$os_months
    ) &
    df[["Survival months flag"]] !=
    "Not calculated because a Death Certificate Only or Autopsy Only case",
    ,
    drop = FALSE
]


## -----------------------------------------------------------------------------
## Grade
## -----------------------------------------------------------------------------

G1 <- c(
    "Well differentiated; Grade I" = 1,
    "Moderately differentiated; Grade II" = 2,
    "Poorly differentiated; Grade III" = 3,
    "Undifferentiated; anaplastic; Grade IV" = 4
)

G2 <- c(
    "Site-specific grade system category (1)" = 1,
    "Site-specific grade system category (2)" = 2,
    "Site-specific grade system category (3)" = 3,
    "Well differentiated" = 1,
    "Low grade" = 1,
    "Moderately differentiated" = 2,
    "Poorly differentiated" = 3,
    "High grade" = 3,
    "Undifferentiated and anaplastic" = 4
)

ga <- unname(
    G1[
        df[["Grade Recode (thru 2017)"]]
    ]
)

gb <- unname(
    G2[
        df[["Derived Summary Grade 2018 (2018+)"]]
    ]
)

df$grade <- ifelse(
    df$year <= 2017,
    ifelse(
        is.na(ga),
        gb,
        ga
    ),
    ifelse(
        is.na(gb),
        ga,
        gb
    )
)


## -----------------------------------------------------------------------------
## Tumor size
## -----------------------------------------------------------------------------

s1 <- suppressWarnings(
    as.numeric(
        df[["CS tumor size (2004-2015)"]]
    )
)

s2 <- suppressWarnings(
    as.numeric(
        df[["Tumor Size Summary (2016+)"]]
    )
)

s1[
    is.na(s1) |
    s1 < 1 |
    s1 > 988
] <- NA_real_

s2[
    is.na(s2) |
    s2 < 1 |
    s2 > 988
] <- NA_real_

df$size_mm <- ifelse(
    is.na(s1),
    s2,
    s1
)


## -----------------------------------------------------------------------------
## Demographics
## -----------------------------------------------------------------------------

df$age <- suppressWarnings(
    as.numeric(
        sub(
            "^\\D*(\\d+).*$",
            "\\1",
            df[["Age recode with single ages and 90+"]]
        )
    )
)

df$female <- as.integer(
    df[["Sex"]] ==
    "Female"
)

df$race_eth <- df[[
    "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)"
]]

df$income <- df[[
    "Median household income inflation adj to 2024"
]]

df$rural_urban <- df[[
    "Rural-Urban Continuum Code"
]]

df$marital <- df[[
    "Marital status at diagnosis"
]]

df$era <- cut(
    df$year,
    breaks = c(
        2003,
        2009,
        2017,
        2023
    ),
    labels = c(
        "2004-2009",
        "2010-2017",
        "2018-2023"
    )
)


## -----------------------------------------------------------------------------
## Nodes
## -----------------------------------------------------------------------------

lx <- suppressWarnings(
    as.numeric(
        df[["Regional nodes examined (1988+)"]]
    )
)

lp <- suppressWarnings(
    as.numeric(
        df[["Regional nodes positive (1988+)"]]
    )
)

df$ln_examined <- ifelse(
    !is.na(lx) &
    lx < 90,
    lx,
    NA_real_
)

df$n_positive <- ifelse(
    !is.na(lp) &
    lp >= 1 &
    lp < 90,
    1,
    ifelse(
        !is.na(lp) &
        lp == 0,
        0,
        NA_real_
    )
)

df$adequate_ln <- as.integer(
    !is.na(
        df$ln_examined
    ) &
    df$ln_examined >= 12
)


## -----------------------------------------------------------------------------
## CEA / chemotherapy / stage / outcomes
## -----------------------------------------------------------------------------

CEA_MAP <- c(
    "CEA negative/normal; within normal limits" = "normal",
    "CEA positive/elevated" = "elevated",
    "Borderline" = "borderline"
)

df$cea <- unname(
    CEA_MAP[
        df[["CEA Pretreatment Interpretation Recode (2010+)"]]
    ]
)

df$chemo <- as.integer(
    df[["Chemotherapy recode (yes, no/unk)"]] ==
    "Yes"
)

df$stage_sum <- ifelse(
    grepl(
        "^Localized",
        df[, CSS]
    ),
    "Localized",
    "Regional"
)

csc <- df[[
    "SEER cause-specific death classification"
]]

df$css_event <- ifelse(
    csc ==
    "Dead (attributable to this cancer dx)",
    1,
    ifelse(
        csc ==
        "Alive or dead of other cause",
        0,
        NA_real_
    )
)

df$os_event <- as.integer(
    df[["Vital status recode (study cutoff used)"]] ==
    "Dead"
)

df$other_death <- as.integer(
    df$os_event == 1 &
    !is.na(
        df$css_event
    ) &
    df$css_event == 0
)


## =============================================================================
## 5. Function creating a 25-column analytic dataset for any surgery definition
## =============================================================================

make_analytic <- function(
    limited_mask,
    colectomy_mask
) {

    eligible <- (
        !is.na(limited_mask) &
        !is.na(colectomy_mask) &
        (
            limited_mask |
            colectomy_mask
        )
    )

    x <- df[
        eligible,
        ,
        drop = FALSE
    ]

    lm <- limited_mask[
        eligible
    ]

    cm <- colectomy_mask[
        eligible
    ]

    ## No observation may be assigned to both arms
    if (any(lm & cm)) {
        stop("Surgery definition has overlapping limited/colectomy codes.")
    }

    x$exposure <- ifelse(
        lm,
        "limited",
        "colectomy"
    )

    out <- data.frame(
        patient_id = x[["Patient ID"]],
        year = x$year,
        era = as.character(x$era),
        age = x$age,
        female = x$female,
        race_eth = x$race_eth,
        income = x$income,
        rural_urban = x$rural_urban,
        marital = x$marital,
        hist_group = x$hist_group,
        grade = x$grade,
        T = x$T_stage,
        size_mm = x$size_mm,
        stage_sum = x$stage_sum,
        surg_code = x$surg_code,
        exposure = x$exposure,
        ln_examined = x$ln_examined,
        n_positive = x$n_positive,
        adequate_ln = x$adequate_ln,
        cea = x$cea,
        chemo = x$chemo,
        os_months = x$os_months,
        os_event = x$os_event,
        css_event = x$css_event,
        other_death = x$other_death,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    out
}


## =============================================================================
## 6. Reconstruct MAIN cohort and hard fingerprint
## =============================================================================

MAIN <- make_analytic(
    limited_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 30 &
        df$surg_code <= 32
    ),
    colectomy_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 40 &
        df$surg_code <= 41
    )
)

FROZEN <- read.csv(
    FROZEN_MAIN_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

EXPECTED_FLOW <- c(
    E0 = 24632,
    E1 = 23021,
    E2 = 13972,
    E3 = 13877,
    E4 = 13868,
    E5 = 8223,
    E6 = 7689
)

FLOW_PASS <- TRUE

for (nm in names(EXPECTED_FLOW)) {

    observed <- FLOW$n[
        FLOW$step == nm
    ]

    if (
        length(observed) != 1 ||
        observed !=
        EXPECTED_FLOW[
            nm
        ]
    ) {

        FLOW_PASS <- FALSE
    }
}

MAIN_N_PASS <- (
    nrow(MAIN) ==
    6951
)

MAIN_ARM_PASS <- (
    sum(
        MAIN$exposure ==
        "limited"
    ) ==
    2877 &&
    sum(
        MAIN$exposure ==
        "colectomy"
    ) ==
    4074
)

## -------------------------------------------------------------------------
## Cohort identity audit
##
## Patient ID is retained as an informational check only because the original
## SEER export identifier and the frozen analytic CSV identifier are not
## guaranteed to be the same serialization.
##
## The hard cohort-identity gate instead compares the complete multiset of the
## 24 analytic fields OTHER THAN patient_id. This is order-independent and
## duplicate-aware.
## -------------------------------------------------------------------------

normalize_character <- function(x) {

    y <- trimws(
        as.character(x)
    )

    y[
        is.na(y) |
        y == ""
    ] <- "<NA>"

    y
}


normalize_numeric <- function(x) {

    y <- suppressWarnings(
        as.numeric(x)
    )

    out <- ifelse(
        is.na(y),
        "<NA>",
        format(
            y,
            scientific = FALSE,
            trim = TRUE,
            digits = 15
        )
    )

    out
}


make_cohort_signature <- function(d) {

    required_cols <- c(
        "year",
        "era",
        "age",
        "female",
        "race_eth",
        "income",
        "rural_urban",
        "marital",
        "hist_group",
        "grade",
        "T",
        "size_mm",
        "stage_sum",
        "surg_code",
        "exposure",
        "ln_examined",
        "n_positive",
        "adequate_ln",
        "cea",
        "chemo",
        "os_months",
        "os_event",
        "css_event",
        "other_death"
    )

    missing_cols <- required_cols[
        !required_cols %in%
        names(d)
    ]

    if (
        length(missing_cols) >
        0
    ) {

        stop(
            "Cohort signature missing columns: ",
            paste(
                missing_cols,
                collapse = " | "
            )
        )
    }

    numeric_cols <- c(
        "year",
        "age",
        "female",
        "grade",
        "size_mm",
        "surg_code",
        "ln_examined",
        "n_positive",
        "adequate_ln",
        "chemo",
        "os_months",
        "os_event",
        "css_event",
        "other_death"
    )

    z <- vector(
        "list",
        length(required_cols)
    )

    names(z) <- required_cols

    for (
        nm in required_cols
    ) {

        if (
            nm %in%
            numeric_cols
        ) {

            z[
                nm
            ] <- list(
                normalize_numeric(
                    d[, nm]
                )
            )

        } else {

            z[
                nm
            ] <- list(
                normalize_character(
                    d[, nm]
                )
            )
        }
    }

    do.call(
        paste,
        c(
            z,
            sep = "\u241F"
        )
    )
}


## Informational patient-ID comparison
main_id <- trimws(
    as.character(
        MAIN$patient_id
    )
)

frozen_id <- trimws(
    as.character(
        FROZEN$patient_id
    )
)

MAIN_DUPLICATE_IDS <- sum(
    duplicated(
        main_id
    )
)

FROZEN_DUPLICATE_IDS <- sum(
    duplicated(
        frozen_id
    )
)

MAIN_ID_SET_PASS <- (
    length(main_id) ==
    length(frozen_id) &&
    setequal(
        main_id,
        frozen_id
    )
)

MAIN_ID_ORDER_PASS <- (
    length(main_id) ==
    length(frozen_id) &&
    identical(
        main_id,
        frozen_id
    )
)


## Hard content-based cohort identity comparison
MAIN_SIGNATURE <- make_cohort_signature(
    MAIN
)

FROZEN_SIGNATURE <- make_cohort_signature(
    FROZEN
)

MAIN_CONTENT_PASS <- identical(
    sort(
        MAIN_SIGNATURE
    ),
    sort(
        FROZEN_SIGNATURE
    )
)


## Additional count audit by year x exposure x histology
MAIN_CROSS_TAB <- as.data.frame(
    xtabs(
        ~ year +
        exposure +
        hist_group,
        data = MAIN
    )
)

FROZEN_CROSS_TAB <- as.data.frame(
    xtabs(
        ~ year +
        exposure +
        hist_group,
        data = FROZEN
    )
)

MAIN_CROSS_TAB <- MAIN_CROSS_TAB[
    order(
        MAIN_CROSS_TAB$year,
        MAIN_CROSS_TAB$exposure,
        MAIN_CROSS_TAB$hist_group
    ),
    ,
    drop = FALSE
]

FROZEN_CROSS_TAB <- FROZEN_CROSS_TAB[
    order(
        FROZEN_CROSS_TAB$year,
        FROZEN_CROSS_TAB$exposure,
        FROZEN_CROSS_TAB$hist_group
    ),
    ,
    drop = FALSE
]

rownames(
    MAIN_CROSS_TAB
) <- NULL

rownames(
    FROZEN_CROSS_TAB
) <- NULL

MAIN_CROSSTAB_PASS <- identical(
    MAIN_CROSS_TAB,
    FROZEN_CROSS_TAB
)

cat("Frozen flow audit:\n")

for (i in seq_len(nrow(FLOW))) {

    cat(
        FLOW$step[i],
        " = ",
        FLOW$n[i],
        "\n",
        sep = ""
    )
}

cat(
    "\nMain reconstructed N = ",
    nrow(MAIN),
    " | limited=",
    sum(MAIN$exposure == "limited"),
    " | colectomy=",
    sum(MAIN$exposure == "colectomy"),
    "\n",
    sep = ""
)

cat(
    "Flow PASS = ",
    FLOW_PASS,
    "\n",
    sep = ""
)

cat(
    "Main N/arms PASS = ",
    MAIN_N_PASS &&
    MAIN_ARM_PASS,
    "\n",
    sep = ""
)

cat(
    "Patient ID informational check: set PASS = ",
    MAIN_ID_SET_PASS,
    " | order PASS = ",
    MAIN_ID_ORDER_PASS,
    " | duplicates MAIN/FROZEN = ",
    MAIN_DUPLICATE_IDS,
    "/",
    FROZEN_DUPLICATE_IDS,
    "\n",
    sep = ""
)

cat(
    "24-field analytic-content multiset PASS = ",
    MAIN_CONTENT_PASS,
    "\n",
    sep = ""
)

cat(
    "Year x exposure x histology crosstab PASS = ",
    MAIN_CROSSTAB_PASS,
    "\n\n",
    sep = ""
)


if (
    !FLOW_PASS ||
    !MAIN_N_PASS ||
    !MAIN_ARM_PASS ||
    !MAIN_CONTENT_PASS ||
    !MAIN_CROSSTAB_PASS
) {

    cat("\n")
    cat("============================================================\n")
    cat("PRE-E7 HARD GATE FAILED\n")
    cat("============================================================\n")
    cat("Flow PASS                  = ", FLOW_PASS, "\n", sep = "")
    cat("Main N PASS                = ", MAIN_N_PASS, "\n", sep = "")
    cat("Main arm-count PASS        = ", MAIN_ARM_PASS, "\n", sep = "")
    cat("24-field content PASS      = ", MAIN_CONTENT_PASS, "\n", sep = "")
    cat("Year/exposure/hist PASS    = ", MAIN_CROSSTAB_PASS, "\n", sep = "")
    cat("Patient-ID set (info only) = ", MAIN_ID_SET_PASS, "\n", sep = "")
    cat("============================================================\n")

    stop(
        "PRE-E7 reconstruction failed analytic-content fingerprint. S2 aborted."
    )
}


## =============================================================================
## 7. Corrected Step 6 framework
##    - Summary Stage removed from PS
##    - PS/OW refitted separately within each predicted-risk stratum
##    - weighted 5-year competing-risk CIF
##    - full-pipeline bootstrap B=1000
## =============================================================================

B <- 1000L
SCRIPT_VERSION <- "A6_REDTEAM06_CORRECTED_FROZEN_20260906_v4_E2FIX"

PS_ORIG <- PS_FORMULA

ps_terms <- attr(
    terms(PS_ORIG),
    "term.labels"
)

stage_terms <- ps_terms[
    grepl(
        "stage",
        ps_terms,
        ignore.case = TRUE
    )
]

if (length(stage_terms) == 0) {
    stop(
        "No stage-related term found in PS_FORMULA. ",
        "Audit the frozen pipeline before continuing."
    )
}

PS_NOSTAGE <- reformulate(
    setdiff(
        ps_terms,
        stage_terms
    ),
    response = all.vars(PS_ORIG)[1],
    env = environment(PS_ORIG)
)

penv <- environment(pipeline)

had_env_ps <- exists(
    "PS_FORMULA",
    envir = penv,
    inherits = FALSE
)

old_env_ps <- if (had_env_ps) {
    get(
        "PS_FORMULA",
        envir = penv,
        inherits = FALSE
    )
} else {
    NULL
}

assign(
    "PS_FORMULA",
    PS_NOSTAGE,
    envir = .GlobalEnv
)

assign(
    "PS_FORMULA",
    PS_NOSTAGE,
    envir = penv
)


## =============================================================================
## 8. Correct weighted competing-risk estimators
## =============================================================================

aj_w <- function(time, ev1, ev2, wt, start = -1, tau = 60) {

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
        !length(time) ||
        sum(wt) <= 0
    ) {
        return(NA_real_)
    }

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
            0:5 / 5,
            na.rm = TRUE,
            names = FALSE
        )
    )

    z[1] <- z[1] - 1e-8
    z[6] <- z[6] + 1e-8

    if (
        any(!is.finite(z)) ||
        any(diff(z) <= 0)
    ) {
        stop("Invalid quintile cut points.")
    }

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
        stop("Subgroup too small.")
    }

    if (
        sum(g$treat == 0) < 15 ||
        sum(g$treat == 1) < 15
    ) {
        stop("Insufficient treatment support.")
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

    if (
        length(ps) != nrow(g) ||
        any(!is.finite(ps))
    ) {
        stop("Non-finite subgroup PS.")
    }

    ps <- pmin(
        pmax(ps, 1e-6),
        1 - 1e-6
    )

    g$ps_sub <- ps

    g$ow_sub <- ifelse(
        g$treat == 1,
        1 - ps,
        ps
    )

    if (
        any(!is.finite(g$ow_sub)) ||
        sum(g$ow_sub) <= 0
    ) {
        stop("Invalid subgroup OW.")
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


LABS <- c(
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
        length(LABS)
    )

    names(out) <- LABS

    for (nm in LABS) {

        g <- z[
            masks[[nm]],
            ,
            drop = FALSE
        ]

        out[nm] <- tryCatch(
            estimate_one(g)["rd"],
            error = function(e) NA_real_
        )
    }

    out
}


ess <- function(w) {
    w <- w[is.finite(w) & w > 0]
    if (!length(w)) return(NA_real_)
    sum(w)^2 / sum(w^2)
}


weighted_smd <- function(x, tr, w) {

    ok <- is.finite(x) & is.finite(tr) & is.finite(w)

    x <- x[ok]
    tr <- tr[ok]
    w <- w[ok]

    if (!sum(tr == 0) || !sum(tr == 1)) {
        return(NA_real_)
    }

    m0 <- weighted.mean(x[tr == 0], w[tr == 0])
    m1 <- weighted.mean(x[tr == 1], w[tr == 1])

    v0 <- weighted.mean((x[tr == 0] - m0)^2, w[tr == 0])
    v1 <- weighted.mean((x[tr == 1] - m1)^2, w[tr == 1])

    den <- sqrt((v0 + v1) / 2)

    if (!is.finite(den) || den <= 0) {
        return(0)
    }

    (m1 - m0) / den
}


balance_design <- function(g) {

    X <- model.matrix(
        delete.response(
            terms(PS_NOSTAGE)
        ),
        data = g,
        na.action = na.pass
    )

    if ("(Intercept)" %in% colnames(X)) {
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


diagnose_all <- function(d) {

    z <- make_outcome(d)
    masks <- make_masks(z)

    rows <- list()

    for (nm in LABS) {

        g <- z[
            masks[[nm]],
            ,
            drop = FALSE
        ]

        gs <- fit_subgroup_ow(g)
        est <- estimate_one(g)

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

        X <- balance_design(gs)

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

        mx <- max(
            abs(smd),
            na.rm = TRUE
        )

        rows[[nm]] <- data.frame(
            analysis = nm,
            n = nrow(gs),
            limited = nrow(a0),
            colectomy = nrow(a1),
            cancer_events_limited =
                sum(a0$cancer == 1, na.rm = TRUE),
            cancer_events_colectomy =
                sum(a1$cancer == 1, na.rm = TRUE),
            CIF_limited_pct = est["limited"] * 100,
            CIF_colectomy_pct = est["colectomy"] * 100,
            RD_pp = est["rd"] * 100,
            ESS_limited = ess(a0$ow_sub),
            ESS_colectomy = ess(a1$ow_sub),
            max_weight = max(gs$ow_sub, na.rm = TRUE),
            max_abs_weighted_SMD = mx,
            worst_variable = names(
                which.max(abs(smd))
            ),
            balance_pass = is.finite(mx) && mx < 0.10,
            stringsAsFactors = FALSE
        )
    }

    do.call(
        rbind,
        rows
    )
}


## =============================================================================
## 9. Corrected main-analysis fingerprint
## =============================================================================

MAIN_ANALYSIS <- pipeline(
    MAIN
)

MAIN_CORRECTED <- estimate_all(
    MAIN_ANALYSIS
)

EXPECTED_CORRECTED_RD_PP <- c(
    Q1 = -0.170196,
    Q2 = -3.604407,
    Q3 = -0.204167,
    Q4 = -4.432059,
    Q5 = -5.066558,
    bottom_decile = 2.053016,
    risk_lt_5pct = 2.926840
)

MAIN_DIFF_PP <- abs(
    MAIN_CORRECTED * 100 -
    EXPECTED_CORRECTED_RD_PP[
        names(MAIN_CORRECTED)
    ]
)

cat("Corrected main-analysis fingerprint (pp):\n")
cat(
    paste(
        sprintf(
            "%s=%+.3f",
            names(MAIN_CORRECTED),
            MAIN_CORRECTED * 100
        ),
        collapse = " | "
    ),
    "\n"
)

cat(
    "Maximum difference from corrected frozen reference = ",
    sprintf(
        "%.3f pp",
        max(MAIN_DIFF_PP, na.rm = TRUE)
    ),
    "\n\n",
    sep = ""
)

if (
    max(MAIN_DIFF_PP, na.rm = TRUE) > 0.20
) {
    stop(
        "Corrected main-analysis fingerprint failed. ",
        "Do not run Step 6 until Step 1B is reproduced."
    )
}


## =============================================================================
## 10. Alternative surgery populations
## =============================================================================

ALT_A <- make_analytic(
    limited_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 20 &
        df$surg_code <= 39
    ),
    colectomy_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 40 &
        df$surg_code <= 79
    )
)

ALT_B <- make_analytic(
    limited_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 30 &
        df$surg_code <= 31
    ),
    colectomy_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 40 &
        df$surg_code <= 41
    )
)

ALT_C <- make_analytic(
    limited_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 30 &
        df$surg_code <= 32
    ),
    colectomy_mask = (
        !is.na(df$surg_code) &
        df$surg_code >= 40 &
        df$surg_code <= 79
    )
)


EXPECTED_N <- c(
    ALT_A = 7450,
    ALT_B = 6676,
    ALT_C = 7115
)

OBSERVED_N <- c(
    ALT_A = nrow(ALT_A),
    ALT_B = nrow(ALT_B),
    ALT_C = nrow(ALT_C)
)

cat("Alternative-surgery cohort sizes:\n")

for (nm in names(EXPECTED_N)) {
    cat(
        nm,
        ": N=",
        OBSERVED_N[nm],
        " | expected=",
        EXPECTED_N[nm],
        "\n",
        sep = ""
    )
}

if (any(OBSERVED_N != EXPECTED_N)) {
    stop(
        "Alternative-surgery cohort fingerprint failed. ",
        "Expected A=7450, B=6676, C=7115."
    )
}


ALT_LIST <- list(
    ALT_A_20_39_vs_40_79 = ALT_A,
    ALT_B_30_31_vs_40_41 = ALT_B,
    ALT_C_30_32_vs_40_79 = ALT_C
)

ALT_LABELS <- c(
    ALT_A_20_39_vs_40_79 =
        "limited 20-39 vs colectomy 40-79",
    ALT_B_30_31_vs_40_41 =
        "limited 30-31 vs colectomy 40-41",
    ALT_C_30_32_vs_40_79 =
        "limited 30-32 vs colectomy 40-79"
)


## =============================================================================
## 11. Corrected point estimates + balance
## =============================================================================

POINTS <- matrix(
    NA_real_,
    nrow = length(ALT_LIST),
    ncol = length(LABS),
    dimnames = list(
        names(ALT_LIST),
        LABS
    )
)

DIAG_LIST <- list()

for (nm in names(ALT_LIST)) {

    cat(
        "\nRunning point estimate: ",
        ALT_LABELS[nm],
        "\n",
        sep = ""
    )

    dx <- pipeline(
        ALT_LIST[[nm]]
    )

    if (nrow(dx) != nrow(ALT_LIST[[nm]])) {
        stop("pipeline() changed row count for ", nm)
    }

    POINTS[nm, ] <- estimate_all(dx)

    dg <- diagnose_all(dx)
    dg$definition <- nm
    dg$description <- unname(ALT_LABELS[nm])
    DIAG_LIST[[nm]] <- dg

    cat(
        paste(
            sprintf(
                "%s=%+.2f pp",
                LABS,
                POINTS[nm, ] * 100
            ),
            collapse = " | "
        ),
        "\n"
    )
}

DIAG <- do.call(
    rbind,
    DIAG_LIST
)


## =============================================================================
## 12. B=1000 full-pipeline bootstrap, with checkpoints
## =============================================================================

SEEDS <- c(
    ALT_A_20_39_vs_40_79 = 26090601L,
    ALT_B_30_31_vs_40_41 = 26090602L,
    ALT_C_30_32_vs_40_79 = 26090603L
)

BOOT_LIST <- list()

for (nm in names(ALT_LIST)) {

    dat_alt <- ALT_LIST[[nm]]

    checkpoint_file <- file.path(
        RESULT_DIR,
        paste0(
            "REDTEAM_06_Checkpoint_",
            nm,
            "_B1000.rds"
        )
    )

    boot <- matrix(
        NA_real_,
        nrow = B,
        ncol = length(LABS),
        dimnames = list(NULL, LABS)
    )

    start_b <- 1L

    if (file.exists(checkpoint_file)) {

        cp <- tryCatch(
            readRDS(checkpoint_file),
            error = function(e) NULL
        )

        checkpoint_ok <- (
            !is.null(cp) &&
            identical(cp$script_version, SCRIPT_VERSION) &&
            identical(cp$definition, nm) &&
            identical(cp$n, nrow(dat_alt)) &&
            is.matrix(cp$boot) &&
            identical(dim(cp$boot), dim(boot)) &&
            identical(colnames(cp$boot), colnames(boot))
        )

        if (checkpoint_ok) {

            boot <- cp$boot

            done <- which(
                apply(
                    boot,
                    1,
                    function(x) any(is.finite(x))
                )
            )

            start_b <- if (length(done)) {
                max(done) + 1L
            } else {
                1L
            }

            if (!is.null(cp$rng_seed)) {
                assign(
                    ".Random.seed",
                    cp$rng_seed,
                    envir = .GlobalEnv
                )
            }

            cat(
                "Checkpoint loaded for ",
                nm,
                ". Resume at ",
                start_b,
                "/",
                B,
                "\n",
                sep = ""
            )
        }
    }

    if (start_b == 1L) {
        set.seed(SEEDS[nm])
    }

    if (start_b <= B) {

        t0 <- Sys.time()

        for (b in seq.int(start_b, B)) {

            db <- dat_alt[
                sample(
                    seq_len(nrow(dat_alt)),
                    nrow(dat_alt),
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
                    rep(NA_real_, length(LABS))
                }
            )

            if (b %% 50 == 0 || b == B) {

                saveRDS(
                    list(
                        script_version = SCRIPT_VERSION,
                        definition = nm,
                        n = nrow(dat_alt),
                        completed_through = b,
                        boot = boot,
                        rng_seed = .Random.seed
                    ),
                    checkpoint_file
                )

                cat(
                    nm,
                    " bootstrap ",
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

    BOOT_LIST[[nm]] <- boot
}


## =============================================================================
## 13. Final summary
## =============================================================================

SUMMARY_ROWS <- list()
BOOT_ROWS <- list()

row_i <- 0L
boot_i <- 0L

for (nm in names(ALT_LIST)) {

    boot <- BOOT_LIST[[nm]]

    boot_i <- boot_i + 1L

    BOOT_ROWS[[boot_i]] <- data.frame(
        definition = nm,
        replicate = seq_len(B),
        boot,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    dg <- DIAG[
        DIAG$definition == nm,
        ,
        drop = FALSE
    ]

    for (j in seq_along(LABS)) {

        lab <- LABS[j]

        v <- boot[, j]
        valid <- is.finite(v)
        effB <- sum(valid)

        if (effB < 900) {
            warning(
                nm,
                " / ",
                lab,
                ": effective B <900 (",
                effB,
                ")."
            )
        }

        lo <- if (effB >= 50) {
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

        hi <- if (effB >= 50) {
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

        zdiag <- dg[
            dg$analysis == lab,
            ,
            drop = FALSE
        ]

        upper_possible_reduction_pp <- if (is.finite(lo)) {
            max(0, -lo * 100)
        } else {
            NA_real_
        }

        is_lowrisk <- lab %in%
            c(
                "Q1",
                "bottom_decile",
                "risk_lt_5pct"
            )

        row_i <- row_i + 1L

        SUMMARY_ROWS[[row_i]] <- data.frame(
            definition = nm,
            description = unname(ALT_LABELS[nm]),
            cohort_n = nrow(ALT_LIST[[nm]]),
            analysis = lab,
            n = zdiag$n,
            limited = zdiag$limited,
            colectomy = zdiag$colectomy,
            cancer_events_limited =
                zdiag$cancer_events_limited,
            cancer_events_colectomy =
                zdiag$cancer_events_colectomy,
            CIF_limited_pct =
                zdiag$CIF_limited_pct,
            CIF_colectomy_pct =
                zdiag$CIF_colectomy_pct,
            RD_pp = POINTS[nm, lab] * 100,
            CI_low_pp = lo * 100,
            CI_high_pp = hi * 100,
            upper_possible_reduction_pp =
                upper_possible_reduction_pp,
            criterion_3pp = if (is_lowrisk) {
                upper_possible_reduction_pp < 3
            } else {
                NA
            },
            effective_B = effB,
            ESS_limited = zdiag$ESS_limited,
            ESS_colectomy = zdiag$ESS_colectomy,
            max_weight = zdiag$max_weight,
            max_abs_weighted_SMD =
                zdiag$max_abs_weighted_SMD,
            worst_variable =
                zdiag$worst_variable,
            balance_pass =
                zdiag$balance_pass,
            stringsAsFactors = FALSE
        )
    }
}

SUMMARY <- do.call(
    rbind,
    SUMMARY_ROWS
)

BOOT_WIDE <- do.call(
    rbind,
    BOOT_ROWS
)


## =============================================================================
## 14. Save final files
## =============================================================================

OUT_SUMMARY <- file.path(
    RESULT_DIR,
    "REDTEAM_06_AlternativeSurgery_B1000.csv"
)

OUT_BOOT <- file.path(
    RESULT_DIR,
    "REDTEAM_06_AlternativeSurgery_Bootstrap_B1000.csv"
)

OUT_REPORT <- file.path(
    RESULT_DIR,
    "REDTEAM_06_AlternativeSurgery_Report.txt"
)

OUT_AUDIT <- file.path(
    RESULT_DIR,
    "REDTEAM_06_AlternativeSurgery_Audit.txt"
)

write.csv(
    SUMMARY,
    OUT_SUMMARY,
    row.names = FALSE
)

write.csv(
    BOOT_WIDE,
    OUT_BOOT,
    row.names = FALSE
)


audit_lines <- c(
    "======================================================================",
    "A6 RED-TEAM STEP 6 — CORRECTED ALTERNATIVE SURGERY DEFINITIONS",
    "======================================================================",
    paste0("Script version: ", SCRIPT_VERSION),
    paste0("Original source: ", SEER_FILE),
    paste0("Frozen analytic cohort: ", FROZEN_MAIN_FILE),
    paste0("Frozen pipeline: ", PIPE_FILE),
    "",
    paste0("Main N PASS: ", MAIN_N_PASS),
    paste0("Main arm-count PASS: ", MAIN_ARM_PASS),
    paste0("24-field analytic-content PASS: ", MAIN_CONTENT_PASS),
    paste0("Year/exposure/histology crosstab PASS: ", MAIN_CROSSTAB_PASS),
    paste0(
        "Main patient-ID membership PASS (diagnostic): ",
        MAIN_ID_SET_PASS
    ),
    paste0(
        "Main patient-ID row-order PASS: ",
        MAIN_ID_ORDER_PASS,
        " (informational only)"
    ),
    paste0(
        "Removed PS term(s): ",
        paste(stage_terms, collapse = " | ")
    ),
    paste0(
        "Corrected Step 1B fingerprint max difference: ",
        sprintf(
            "%.3f pp",
            max(MAIN_DIFF_PP, na.rm = TRUE)
        )
    ),
    "",
    paste0("Alternative A N=", nrow(ALT_A)),
    paste0("Alternative B N=", nrow(ALT_B)),
    paste0("Alternative C N=", nrow(ALT_C)),
    "",
    "Framework:",
    "- Summary Stage excluded from PS.",
    "- PS/OW refitted separately within every predicted-risk stratum.",
    "- Nodal-risk model re-estimated through pipeline().",
    "- Quintile/bottom-decile cut points re-estimated.",
    "- 5-year cancer-specific mortality uses weighted competing-risk CIF.",
    "- Other-cause death is the competing event.",
    "- Full-pipeline percentile bootstrap B=1000.",
    "======================================================================"
)

writeLines(
    audit_lines,
    OUT_AUDIT,
    useBytes = TRUE
)


report <- c(
    "======================================================================",
    "A6 RED-TEAM STEP 6 — ALTERNATIVE SURGERY DEFINITIONS, CORRECTED B=1000",
    "======================================================================",
    "RD = oncologic colectomy - limited resection.",
    "RD < 0 indicates lower 5-year cancer-specific mortality associated with colectomy.",
    "",
    paste0(
        "Patient-ID membership PASS: ",
        MAIN_ID_SET_PASS,
        " | row-order PASS: ",
        MAIN_ID_ORDER_PASS,
        " (row order is not a cohort criterion)"
    ),
    paste0(
        "Removed PS term(s): ",
        paste(stage_terms, collapse = " | ")
    ),
    ""
)

for (nm in names(ALT_LIST)) {

    report <- c(
        report,
        paste0(
            "-- ",
            ALT_LABELS[nm],
            " | cohort N=",
            nrow(ALT_LIST[[nm]]),
            " --"
        )
    )

    ss <- SUMMARY[
        SUMMARY$definition == nm,
        ,
        drop = FALSE
    ]

    for (i in seq_len(nrow(ss))) {

        crit_txt <- ""

        if (
            ss$analysis[i] %in%
            c(
                "Q1",
                "bottom_decile",
                "risk_lt_5pct"
            )
        ) {
            crit_txt <- paste0(
                " | upper compatible colectomy reduction ",
                sprintf(
                    "%.2f",
                    ss$upper_possible_reduction_pp[i]
                ),
                " pp | 3-pp criterion: ",
                ifelse(
                    isTRUE(ss$criterion_3pp[i]),
                    "MET",
                    "NOT MET"
                )
            )
        }

        report <- c(
            report,
            paste0(
                sprintf(
                    "  %-15s RD %+.2f pp (95%% CI %+.2f to %+.2f; B=%d)",
                    ss$analysis[i],
                    ss$RD_pp[i],
                    ss$CI_low_pp[i],
                    ss$CI_high_pp[i],
                    ss$effective_B[i]
                ),
                sprintf(
                    " | CIF %.2f%% vs %.2f%% | max|SMD| %.3f [%s]",
                    ss$CIF_limited_pct[i],
                    ss$CIF_colectomy_pct[i],
                    ss$max_abs_weighted_SMD[i],
                    ifelse(
                        ss$balance_pass[i],
                        "PASS",
                        "FAIL"
                    )
                ),
                crit_txt
            )
        )
    }

    report <- c(
        report,
        ""
    )
}

report <- c(
    report,
    "Interpretation:",
    "These sensitivity analyses change only the operational surgery-code definition.",
    "They preserve the corrected Step 1B framework and the frozen 5-year competing-risk estimand.",
    "",
    "IMPORTANT:",
    "Do not use the temporary +4.18/+4.67 pp event-proportion results in manuscript text or figures.",
    "======================================================================"
)

writeLines(
    report,
    OUT_REPORT,
    useBytes = TRUE
)


## =============================================================================
## 15. Restore original formula and console output
## =============================================================================

assign(
    "PS_FORMULA",
    PS_ORIG,
    envir = .GlobalEnv
)

if (had_env_ps) {
    assign(
        "PS_FORMULA",
        old_env_ps,
        envir = penv
    )
}

cat("\n============================================================\n")
cat("A6 RED-TEAM STEP 6 CORRECTED COMPLETE\n")
cat("============================================================\n")
cat(
    paste(
        report,
        collapse = "\n"
    ),
    "\n"
)

cat("\nFiles written:\n")
cat(OUT_SUMMARY, "\n")
cat(OUT_BOOT, "\n")
cat(OUT_REPORT, "\n")
cat(OUT_AUDIT, "\n")
cat("============================================================\n")
