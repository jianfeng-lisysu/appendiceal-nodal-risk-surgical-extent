## =============================================================================
## A6 - GLOBAL CONSISTENCY AUDIT
## Figures 1-4 + Main Tables 1-2 + Supplementary Tables S1-S9
## READ-ONLY AUDIT — no file is modified, moved, renamed, or deleted
##
## Root:
##   <project root>
##
## OUTPUTS:
##   06_audit/00_GLOBAL_CONSISTENCY_AUDIT_FINAL.txt
##   06_audit/00_GLOBAL_FILE_MANIFEST_FINAL.csv
##   06_audit/00_GLOBAL_NUMERIC_CHECKS_FINAL.csv
##
## PURPOSE:
##   1. Check required final analysis-source files
##   2. Check Figure 1-4 / Table 1-2 / Supplementary Table S1-S9 existence
##   3. Check the frozen cohort/model/main-result fingerprints
##   4. Check Supplementary Table source consistency
##   5. Check PNG/TIFF figure-pair availability
##   6. If 'magick' is installed, audit image dimensions / density metadata
##
## IMPORTANT:
##   This script is intentionally READ-ONLY.
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

if (!dir.exists(ROOT)) {
    stop("A6 root does not exist: ", ROOT)
}

AUDIT_DIR <- file.path(
    ROOT,
    "06_audit"
)

dir.create(
    AUDIT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

REPORT_FILE <- file.path(
    AUDIT_DIR,
    "00_GLOBAL_CONSISTENCY_AUDIT_FINAL.txt"
)

MANIFEST_FILE <- file.path(
    AUDIT_DIR,
    "00_GLOBAL_FILE_MANIFEST_FINAL.csv"
)

CHECK_FILE <- file.path(
    AUDIT_DIR,
    "00_GLOBAL_NUMERIC_CHECKS_FINAL.csv"
)


## =============================================================================
## 1. Recursive file inventory
## =============================================================================

all_files <- list.files(
    ROOT,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    include.dirs = FALSE
)

all_files <- unique(
    normalizePath(
        all_files,
        winslash = "/",
        mustWork = FALSE
    )
)

rel_path <- function(x) {

    root_norm <- normalizePath(
        ROOT,
        winslash = "/",
        mustWork = FALSE
    )

    sub(
        paste0(
            "^",
            gsub(
                "([][{}()+*^$|\\\\?.])",
                "\\\\\\1",
                root_norm
            ),
            "/?"
        ),
        "",
        x
    )
}

manifest <- data.frame(
    path = all_files,
    relative_path = vapply(
        all_files,
        rel_path,
        character(1)
    ),
    filename = basename(
        all_files
    ),
    extension = tolower(
        tools::file_ext(
            all_files
        )
    ),
    size_bytes = file.info(
        all_files
    )$size,
    modified_time = as.character(
        file.info(
            all_files
        )$mtime
    ),
    stringsAsFactors = FALSE
)

write.csv(
    manifest,
    MANIFEST_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 2. Helper functions
## =============================================================================

find_files <- function(
    pattern,
    ignore_case = TRUE
) {

    hit <- grepl(
        pattern,
        manifest$filename,
        ignore.case = ignore_case,
        perl = TRUE
    )

    manifest$path[hit]
}


find_exact <- function(
    filename
) {

    hit <- manifest$filename == filename

    manifest$path[hit]
}


pick_preferred <- function(
    filenames
) {

    found <- character(0)

    for (nm in filenames) {

        z <- find_exact(
            nm
        )

        if (length(z) > 0) {
            found <- c(
                found,
                z
            )
        }
    }

    unique(
        found
    )
}


safe_read <- function(path) {

    if (
        length(path) == 0 ||
        !file.exists(path[1])
    ) {
        return("")
    }

    x <- tryCatch(
        readLines(
            path[1],
            warn = FALSE,
            encoding = "UTF-8"
        ),
        error = function(e) {
            tryCatch(
                readLines(
                    path[1],
                    warn = FALSE
                ),
                error = function(e2) character(0)
            )
        }
    )

    paste(
        x,
        collapse = "\n"
    )
}


safe_csv <- function(path) {

    if (
        length(path) == 0 ||
        !file.exists(path[1])
    ) {
        return(NULL)
    }

    tryCatch(
        read.csv(
            path[1],
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        error = function(e) NULL
    )
}


check_close <- function(
    observed,
    expected,
    tolerance,
    label,
    source
) {

    pass <- (
        length(observed) == 1 &&
        is.finite(observed) &&
        abs(observed - expected) <= tolerance
    )

    data.frame(
        category = "numeric",
        label = label,
        observed = ifelse(
            length(observed) == 1 && is.finite(observed),
            format(
                observed,
                digits = 12,
                scientific = FALSE,
                trim = TRUE
            ),
            "NA"
        ),
        expected = format(
            expected,
            digits = 12,
            scientific = FALSE,
            trim = TRUE
        ),
        tolerance = tolerance,
        status = ifelse(
            pass,
            "PASS",
            "FAIL"
        ),
        source = source,
        stringsAsFactors = FALSE
    )
}


check_equal <- function(
    observed,
    expected,
    label,
    source
) {

    pass <- identical(
        as.character(observed),
        as.character(expected)
    )

    data.frame(
        category = "exact",
        label = label,
        observed = paste(
            observed,
            collapse = " | "
        ),
        expected = paste(
            expected,
            collapse = " | "
        ),
        tolerance = NA_real_,
        status = ifelse(
            pass,
            "PASS",
            "FAIL"
        ),
        source = source,
        stringsAsFactors = FALSE
    )
}


check_contains <- function(
    text,
    token,
    label,
    source
) {

    pass <- grepl(
        token,
        text,
        fixed = TRUE
    )

    data.frame(
        category = "text",
        label = label,
        observed = ifelse(
            pass,
            token,
            "<not found>"
        ),
        expected = token,
        tolerance = NA_real_,
        status = ifelse(
            pass,
            "PASS",
            "FAIL"
        ),
        source = source,
        stringsAsFactors = FALSE
    )
}


file_check <- function(
    label,
    patterns
) {

    hits <- character(0)

    for (pat in patterns) {

        hits <- c(
            hits,
            find_files(
                pat
            )
        )
    }

    hits <- unique(
        hits
    )

    data.frame(
        category = "file",
        label = label,
        observed = ifelse(
            length(hits) > 0,
            paste(
                vapply(
                    hits,
                    rel_path,
                    character(1)
                ),
                collapse = " | "
            ),
            "<missing>"
        ),
        expected = "file exists",
        tolerance = NA_real_,
        status = ifelse(
            length(hits) > 0,
            "PASS",
            "MISSING"
        ),
        source = ROOT,
        stringsAsFactors = FALSE
    )
}


append_check <- function(x) {

    checks[
        length(checks) + 1
    ] <<- list(x)
}


## =============================================================================
## 3. Locate core frozen source files
## =============================================================================

RAW_FILE <- pick_preferred(
    c(
        "analytic_cohort_A.csv"
    )
)

RISK_FILE <- pick_preferred(
    c(
        "analytic_cohort_A_riskscored.csv"
    )
)

WT_FILE <- pick_preferred(
    c(
        "analytic_cohort_A_weighted.csv"
    )
)

PIPE_FILE <- pick_preferred(
    c(
        "Aline_frozen_pipeline_v2.RData"
    )
)

MAIN_REPORT <- pick_preferred(
    c(
        "04_main_report_FINAL_B1000.txt"
    )
)

CONT_REPORT <- pick_preferred(
    c(
        "04b_continuous_report_FINAL_B1000.txt"
    )
)

CONT_CURVE <- pick_preferred(
    c(
        "04b_rd_curve_FINAL_B1000.csv"
    )
)

LANDMARK_REPORT <- pick_preferred(
    c(
        "05_landmark_report_FINAL_B1000.txt"
    )
)

CORE_SENS_REPORT <- pick_preferred(
    c(
        "06_sensitivity_core_report_FINAL.txt"
    )
)

NEGCTRL_REPORT <- pick_preferred(
    c(
        "06_negative_control_othercause_FINAL_B1000.txt"
    )
)

HIST_REPORT <- pick_preferred(
    c(
        "06_histology_subgroup_reweight_FINAL_B1000.txt"
    )
)

ALT_SURGERY_REPORT <- pick_preferred(
    c(
        "06_S2_alternative_surgery_FINAL_B1000.txt"
    )
)

MIXED_CSV <- pick_preferred(
    c(
        "06_negative_control_histology_era_FINAL_B1000.csv"
    )
)


## =============================================================================
## 4. File-presence audit
## =============================================================================

checks <- list()

append_check(
    file_check(
        "Frozen analytic cohort",
        c(
            "^analytic_cohort_A\\.csv$"
        )
    )
)

append_check(
    file_check(
        "Frozen risk-scored cohort",
        c(
            "^analytic_cohort_A_riskscored\\.csv$"
        )
    )
)

append_check(
    file_check(
        "Frozen weighted cohort",
        c(
            "^analytic_cohort_A_weighted\\.csv$"
        )
    )
)

append_check(
    file_check(
        "Frozen pipeline RData",
        c(
            "^Aline_frozen_pipeline_v2\\.RData$"
        )
    )
)

append_check(
    file_check(
        "Final main B1000 report",
        c(
            "^04_main_report_FINAL_B1000.*\\.txt$"
        )
    )
)

append_check(
    file_check(
        "Final continuous-risk B1000 report",
        c(
            "^04b_continuous_report_FINAL_B1000.*\\.txt$"
        )
    )
)

append_check(
    file_check(
        "Final landmark B1000 report",
        c(
            "^05_landmark_report_FINAL_B1000.*\\.txt$"
        )
    )
)

append_check(
    file_check(
        "Final sensitivity core report",
        c(
            "^06_sensitivity_core_report_FINAL.*\\.txt$"
        )
    )
)


## =============================================================================
## 5. Frozen cohort numeric fingerprints
## =============================================================================

raw <- safe_csv(
    RAW_FILE
)

if (!is.null(raw)) {

    append_check(
        check_close(
            nrow(raw),
            6951,
            0,
            "Final analytic cohort N",
            rel_path(RAW_FILE[1])
        )
    )

    if ("exposure" %in% names(raw)) {

        append_check(
            check_close(
                sum(
                    raw$exposure == "limited",
                    na.rm = TRUE
                ),
                2877,
                0,
                "Limited-resection N",
                rel_path(RAW_FILE[1])
            )
        )

        append_check(
            check_close(
                sum(
                    raw$exposure == "colectomy",
                    na.rm = TRUE
                ),
                4074,
                0,
                "Oncologic-colectomy N",
                rel_path(RAW_FILE[1])
            )
        )
    }

    if ("css_event" %in% names(raw)) {

        append_check(
            check_close(
                sum(
                    raw$css_event == 1,
                    na.rm = TRUE
                ),
                1212,
                0,
                "Cancer-specific death count",
                rel_path(RAW_FILE[1])
            )
        )
    }

    if ("os_event" %in% names(raw)) {

        append_check(
            check_close(
                sum(
                    raw$os_event == 1,
                    na.rm = TRUE
                ),
                2111,
                0,
                "Overall death count",
                rel_path(RAW_FILE[1])
            )
        )
    }
}


## =============================================================================
## 6. Main-result text fingerprints
## =============================================================================

main_text <- safe_read(
    MAIN_REPORT
)

if (nzchar(main_text)) {

    main_tokens <- c(
        "Q1             : +0.3 pp (95% CI -3.4 to +3.0)",
        "Q2             : -3.6 pp (95% CI -7.8 to +1.6)",
        "Q3             : -0.7 pp (95% CI -7.4 to +3.8)",
        "Q4             : -5.1 pp (95% CI -9.6 to +2.1)",
        "Q5             : -6.0 pp (95% CI -12.7 to +0.5)",
        "bottom_decile  : +3.1 pp (95% CI -2.2 to +7.0)",
        "risk_lt_5pct   : +4.2 pp (95% CI -2.0 to +7.8)",
        "Q5-Q1 contrast: -6.3 pp (95% CI -13.5 to +1.6)",
        "bottom-decile upper 2.2 [MET]",
        "risk<5% upper 2.0 [MET]",
        "Q5 cancer-specific mortality RR = 0.847",
        "E-value(point)=1.64",
        "E-value(CI-limit)=1.00"
    )

    for (tok in main_tokens) {

        append_check(
            check_contains(
                main_text,
                tok,
                paste0(
                    "Main result token: ",
                    tok
                ),
                rel_path(MAIN_REPORT[1])
            )
        )
    }
}


## =============================================================================
## 7. Continuous-risk fingerprints
## =============================================================================

cont_text <- safe_read(
    CONT_REPORT
)

if (nzchar(cont_text)) {

    possible_p_tokens <- c(
        "0.226",
        "0.2264"
    )

    p_pass <- any(
        vapply(
            possible_p_tokens,
            function(z) grepl(
                z,
                cont_text,
                fixed = TRUE
            ),
            logical(1)
        )
    )

    append_check(
        data.frame(
            category = "text",
            label = "Continuous-risk global interaction P",
            observed = ifelse(
                p_pass,
                "0.226 / 0.2264 found",
                "<not found>"
            ),
            expected = "approximately 0.226",
            tolerance = NA_real_,
            status = ifelse(
                p_pass,
                "PASS",
                "FAIL"
            ),
            source = ifelse(
                length(CONT_REPORT) > 0,
                rel_path(CONT_REPORT[1]),
                "<missing>"
            ),
            stringsAsFactors = FALSE
        )
    )
}


curve <- safe_csv(
    CONT_CURVE
)

if (!is.null(curve)) {

    possible_risk_names <- c(
        "risk",
        "nodal_risk",
        "grid",
        "risk_grid"
    )

    possible_rd_names <- c(
        "rd",
        "estimate",
        "rd_point"
    )

    risk_name <- possible_risk_names[
        possible_risk_names %in%
        names(curve)
    ]

    rd_name <- possible_rd_names[
        possible_rd_names %in%
        names(curve)
    ]

    if (
        length(risk_name) > 0 &&
        length(rd_name) > 0
    ) {

        risk_values <- curve[
            ,
            risk_name[1]
        ]

        rd_values <- curve[
            ,
            rd_name[1]
        ]

        idx5 <- which.min(
            abs(
                risk_values -
                0.05
            )
        )

        if (length(idx5) == 1) {

            observed_rd5 <- rd_values[idx5]

            if (
                max(
                    abs(
                        rd_values
                    ),
                    na.rm = TRUE
                ) <
                1
            ) {
                observed_rd5 <- observed_rd5 * 100
            }

            append_check(
                check_close(
                    observed_rd5,
                    0.3,
                    0.5,
                    "Continuous curve RD near 5% nodal risk, pp",
                    rel_path(CONT_CURVE[1])
                )
            )
        }
    }
}


## =============================================================================
## 8. Landmark fingerprints
## =============================================================================

landmark_text <- safe_read(
    LANDMARK_REPORT
)

if (nzchar(landmark_text)) {

    landmark_tokens <- c(
        "Landmark 3 months (effective B=1000)",
        "bottom_decile  : +2.2 pp (95% CI -2.3 to +5.8)",
        "risk_lt_5pct   : +3.2 pp (95% CI -2.3 to +6.7)",
        "Landmark 6 months (effective B=1000)",
        "bottom_decile  : +2.1 pp (95% CI -2.3 to +5.5)",
        "risk_lt_5pct   : +3.1 pp (95% CI -2.5 to +6.6)",
        "Landmark 12 months (effective B=1000)",
        "bottom_decile  : +1.9 pp (95% CI -2.7 to +5.3)",
        "risk_lt_5pct   : +2.8 pp (95% CI -2.5 to +6.5)"
    )

    for (tok in landmark_tokens) {

        append_check(
            check_contains(
                landmark_text,
                tok,
                paste0(
                    "Landmark token: ",
                    tok
                ),
                rel_path(LANDMARK_REPORT[1])
            )
        )
    }
}


## =============================================================================
## 9. Adequate nodal assessment sensitivity fingerprints
## =============================================================================

sens_text <- safe_read(
    CORE_SENS_REPORT
)

if (nzchar(sens_text)) {

    sens_tokens <- c(
        "colectomy arm restricted to LN>=12 (input n=6130; effective B=500)",
        "Q1: -0.6 pp (95% CI -4.0 to +2.7)",
        "Q2: -5.0 pp (95% CI -9.3 to -0.2)",
        "Q3: -1.5 pp (95% CI -7.9 to +3.2)",
        "Q4: -6.6 pp (95% CI -10.7 to +1.2)",
        "Q5: -7.0 pp (95% CI -13.7 to -0.1)"
    )

    for (tok in sens_tokens) {

        append_check(
            check_contains(
                sens_text,
                tok,
                paste0(
                    "Adequate-node sensitivity: ",
                    tok
                ),
                rel_path(CORE_SENS_REPORT[1])
            )
        )
    }
}


## =============================================================================
## 10. Supplementary table files S1-S9
## =============================================================================

for (s in 1:9) {

    label <- paste0(
        "Supplementary Table S",
        s,
        " DOCX"
    )

    pattern <- paste0(
        "Supplementary_Table_S",
        s,
        ".*FINAL.*\\.docx$"
    )

    append_check(
        file_check(
            label,
            c(
                pattern
            )
        )
    )
}


## =============================================================================
## 11. Main Table 1-2 file audit
## =============================================================================

append_check(
    file_check(
        "Main Table 1",
        c(
            "^Table_?1.*\\.docx$",
            "^Table1.*\\.docx$",
            "^Main_Table_?1.*\\.docx$"
        )
    )
)

append_check(
    file_check(
        "Main Table 2",
        c(
            "^Table_?2.*\\.docx$",
            "^Table2.*\\.docx$",
            "^Main_Table_?2.*\\.docx$"
        )
    )
)


## =============================================================================
## 12. Figure 1-4 file audit
## =============================================================================

for (f in 1:4) {

    append_check(
        file_check(
            paste0(
                "Figure ",
                f,
                " PNG"
            ),
            c(
                paste0(
                    "^Figure_?",
                    f,
                    ".*\\.png$"
                ),
                paste0(
                    "^Fig_?",
                    f,
                    ".*\\.png$"
                )
            )
        )
    )

    append_check(
        file_check(
            paste0(
                "Figure ",
                f,
                " TIFF"
            ),
            c(
                paste0(
                    "^Figure_?",
                    f,
                    ".*\\.tif$"
                ),
                paste0(
                    "^Figure_?",
                    f,
                    ".*\\.tiff$"
                ),
                paste0(
                    "^Fig_?",
                    f,
                    ".*\\.tif$"
                ),
                paste0(
                    "^Fig_?",
                    f,
                    ".*\\.tiff$"
                )
            )
        )
    )
}


## =============================================================================
## 13. Optional image metadata audit
## =============================================================================

image_audit_lines <- c(
    "== IMAGE METADATA AUDIT =="
)

has_magick <- requireNamespace(
    "magick",
    quietly = TRUE
)

if (!has_magick) {

    image_audit_lines <- c(
        image_audit_lines,
        "Package 'magick' is not installed; image density metadata was not checked.",
        "File-pair existence checks above remain valid."
    )

} else {

    image_files <- manifest$path[
        manifest$extension %in%
        c(
            "png",
            "tif",
            "tiff"
        )
    ]

    fig_image_files <- image_files[
        grepl(
            "^(Figure|Fig)_?[1-4]",
            basename(
                image_files
            ),
            ignore.case = TRUE
        )
    ]

    if (length(fig_image_files) == 0) {

        image_audit_lines <- c(
            image_audit_lines,
            "No Figure 1-4 PNG/TIFF files found for metadata inspection."
        )

    } else {

        for (img_path in fig_image_files) {

            info <- tryCatch(
                magick::image_info(
                    magick::image_read(
                        img_path
                    )
                ),
                error = function(e) NULL
            )

            if (is.null(info)) {

                image_audit_lines <- c(
                    image_audit_lines,
                    paste0(
                        rel_path(
                            img_path
                        ),
                        " | metadata read FAILED"
                    )
                )

            } else {

                density_value <- ""

                if ("density" %in% names(info)) {
                    density_value <- as.character(
                        info$density[1]
                    )
                }

                image_audit_lines <- c(
                    image_audit_lines,
                    paste0(
                        rel_path(
                            img_path
                        ),
                        " | ",
                        info$width[1],
                        "x",
                        info$height[1],
                        " px",
                        ifelse(
                            nzchar(
                                density_value
                            ),
                            paste0(
                                " | density=",
                                density_value
                            ),
                            ""
                        )
                    )
                )
            }
        }
    }
}


## =============================================================================
## 14. Optional DOCX title/text audit using officer
## =============================================================================

docx_audit_lines <- c(
    "== DOCX TEXT AUDIT =="
)

has_officer <- requireNamespace(
    "officer",
    quietly = TRUE
)

if (!has_officer) {

    docx_audit_lines <- c(
        docx_audit_lines,
        "Package 'officer' is not installed; DOCX text was not parsed."
    )

} else {

    supp_docx <- character(0)

    for (s in 1:9) {

        z <- find_files(
            paste0(
                "^Supplementary_Table_S",
                s,
                ".*FINAL.*\\.docx$"
            )
        )

        if (length(z) > 0) {

            ## Prefer names containing FINAL_LOCK or COMPACT_FINAL.
            ord <- order(
                !grepl(
                    "FINAL_LOCK|COMPACT_FINAL",
                    basename(
                        z
                    ),
                    ignore.case = TRUE
                ),
                basename(
                    z
                )
            )

            supp_docx <- c(
                supp_docx,
                z[ord][1]
            )
        }
    }

    if (length(supp_docx) == 0) {

        docx_audit_lines <- c(
            docx_audit_lines,
            "No Supplementary Table S1-S9 DOCX files found."
        )

    } else {

        for (docx_path in supp_docx) {

            sx <- tryCatch(
                officer::docx_summary(
                    officer::read_docx(
                        docx_path
                    )
                ),
                error = function(e) NULL
            )

            if (is.null(sx)) {

                docx_audit_lines <- c(
                    docx_audit_lines,
                    paste0(
                        rel_path(
                            docx_path
                        ),
                        " | parse FAILED"
                    )
                )

            } else {

                text_all <- paste(
                    sx$text,
                    collapse = " "
                )

                title_guess <- sub(
                    "^(.*?Supplementary Table S[0-9]+\\..{0,180}).*$",
                    "\\1",
                    text_all,
                    perl = TRUE
                )

                docx_audit_lines <- c(
                    docx_audit_lines,
                    paste0(
                        rel_path(
                            docx_path
                        ),
                        " | chars=",
                        nchar(
                            text_all
                        ),
                        " | ",
                        substr(
                            title_guess,
                            1,
                            220
                        )
                    )
                )
            }
        }
    }
}


## =============================================================================
## 15. Build combined numeric/file check table
## =============================================================================

check_table <- do.call(
    rbind,
    checks
)

write.csv(
    check_table,
    CHECK_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)


## =============================================================================
## 16. Summary counts
## =============================================================================

n_pass <- sum(
    check_table$status == "PASS",
    na.rm = TRUE
)

n_fail <- sum(
    check_table$status == "FAIL",
    na.rm = TRUE
)

n_missing <- sum(
    check_table$status == "MISSING",
    na.rm = TRUE
)

GLOBAL_PASS <- (
    n_fail == 0 &&
    n_missing == 0
)


## =============================================================================
## 17. Human-readable report
## =============================================================================

report <- c(
    "======================================================================",
    "A6 GLOBAL CONSISTENCY AUDIT — FINAL",
    "======================================================================",
    paste0(
        "Timestamp: ",
        Sys.time()
    ),
    paste0(
        "Root: ",
        ROOT
    ),
    "",
    paste0(
        "PASS = ",
        n_pass
    ),
    paste0(
        "FAIL = ",
        n_fail
    ),
    paste0(
        "MISSING = ",
        n_missing
    ),
    paste0(
        "GLOBAL STATUS = ",
        ifelse(
            GLOBAL_PASS,
            "PASS",
            "REVIEW REQUIRED"
        )
    ),
    "",
    "== FAILED / MISSING CHECKS =="
)

bad <- check_table[
    check_table$status %in%
    c(
        "FAIL",
        "MISSING"
    ),
    ,
    drop = FALSE
]

if (nrow(bad) == 0) {

    report <- c(
        report,
        "None."
    )

} else {

    for (i in seq_len(nrow(bad))) {

        report <- c(
            report,
            paste0(
                "[",
                bad$status[i],
                "] ",
                bad$label[i],
                " | observed=",
                bad$observed[i],
                " | expected=",
                bad$expected[i],
                " | source=",
                bad$source[i]
            )
        )
    }
}


report <- c(
    report,
    "",
    "== ALL CHECKS =="
)

for (i in seq_len(nrow(check_table))) {

    report <- c(
        report,
        paste0(
            "[",
            check_table$status[i],
            "] ",
            check_table$label[i],
            " | observed=",
            check_table$observed[i],
            " | expected=",
            check_table$expected[i]
        )
    )
}


report <- c(
    report,
    "",
    image_audit_lines,
    "",
    docx_audit_lines,
    "",
    "== INTERPRETATION ==",
    paste0(
        "A PASS here means the machine-readable frozen sources and expected final ",
        "artifacts are internally consistent on the checks implemented by this script."
    ),
    paste0(
        "It does not replace a final human visual review of figure labels, legends, ",
        "font sizes, cropping, line breaks, and journal-specific formatting."
    ),
    "",
    "======================================================================"
)


writeLines(
    report,
    REPORT_FILE,
    useBytes = TRUE
)


## =============================================================================
## 18. Console
## =============================================================================

cat(
    paste(
        report[
            seq_len(
                min(
                    length(
                        report
                    ),
                    80
                )
            )
        ],
        collapse = "\n"
    ),
    "\n"
)

if (length(report) > 80) {
    cat(
        "\n... full report saved to:\n",
        REPORT_FILE,
        "\n",
        sep = ""
    )
}

cat("\n")
cat("Manifest: ", MANIFEST_FILE, "\n", sep = "")
cat("Checks  : ", CHECK_FILE, "\n", sep = "")
cat("Report  : ", REPORT_FILE, "\n", sep = "")

cat("\n")
cat("============================================================\n")
cat(
    "GLOBAL STATUS: ",
    ifelse(
        GLOBAL_PASS,
        "PASS",
        "REVIEW REQUIRED"
    ),
    "\n",
    sep = ""
)
cat("============================================================\n")
