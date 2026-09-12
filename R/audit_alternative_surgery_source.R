## =============================================================================
## A6 - S7 SOURCE AUDIT ONLY
## Locate existing alternative-surgery-definition sensitivity results
##
## Root:
##   <project root>
##
## PURPOSE
##   This is NOT an analysis script.
##   It does NOT reconstruct PRE-E7 data and does NOT re-run any model.
##   It only searches existing result/report files for the already-completed
##   alternative surgery-definition sensitivity analyses.
##
## OUTPUT
##   03_results/S7_alternative_surgery_source_audit.txt
##
## EXPECTED historical sensitivity definitions
##   A. Limited 20–39 vs colectomy >=40
##   B. Excluding surgery code 32
##   C. Colectomy broadened to include 50–79
##
## IMPORTANT
##   These populations cannot be reconstructed from analytic_cohort_A.csv
##   (N=6,951), because some alternative surgery codes lie outside the
##   frozen primary E7 cohort.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
RESULT_DIR <- file.path(ROOT, "03_results")

dir.create(
    RESULT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

OUT_FILE <- file.path(
    RESULT_DIR,
    "S7_alternative_surgery_source_audit.txt"
)

## =============================================================================
## 1. Candidate directories
## =============================================================================

dirs <- unique(
    c(
        RESULT_DIR,
        ROOT
    )
)

all_files <- character(0)

for (dd in dirs) {

    if (dir.exists(dd)) {

        ff <- list.files(
            dd,
            full.names = TRUE,
            recursive = FALSE,
            all.files = FALSE
        )

        all_files <- c(
            all_files,
            ff
        )
    }
}

all_files <- unique(all_files)

## =============================================================================
## 2. Candidate file-name filter
## =============================================================================

name_pat <- paste(
    c(
        "06",
        "sensitivity",
        "surgery",
        "alternative",
        "definition",
        "robust",
        "code32",
        "code_32",
        "50-79",
        "50_79",
        "20-39",
        "20_39"
    ),
    collapse = "|"
)

candidate_files <- all_files[
    grepl(
        name_pat,
        basename(all_files),
        ignore.case = TRUE
    )
]

## Keep readable text-like sources first.
text_ext <- "\\.(txt|csv|tsv|md|log)$"

text_candidates <- candidate_files[
    grepl(
        text_ext,
        candidate_files,
        ignore.case = TRUE
    )
]

## =============================================================================
## 3. Search keywords inside readable files
## =============================================================================

content_keywords <- c(
    "20–39",
    "20-39",
    ">=40",
    "≥40",
    "code 32",
    "code32",
    "50–79",
    "50-79",
    "alternative surgery",
    "surgery definition",
    "exclude code 32",
    "excluding code 32",
    "colectomy includes",
    "broaden",
    "7450",
    "6676",
    "7115"
)

safe_read_lines <- function(path) {

    out <- tryCatch(
        readLines(
            path,
            warn = FALSE,
            encoding = "UTF-8"
        ),
        error = function(e) {
            character(0)
        }
    )

    if (length(out) == 0) {

        out <- tryCatch(
            readLines(
                path,
                warn = FALSE
            ),
            error = function(e) {
                character(0)
            }
        )
    }

    out
}

hits <- list()
hit_counter <- 0

for (path in text_candidates) {

    x <- safe_read_lines(path)

    if (length(x) == 0) {
        next
    }

    match_idx <- integer(0)

    for (kw in content_keywords) {

        match_idx <- union(
            match_idx,
            grep(
                kw,
                x,
                ignore.case = TRUE,
                fixed = TRUE
            )
        )
    }

    if (length(match_idx) > 0) {

        ## Add context ±3 lines.
        context_idx <- sort(
            unique(
                unlist(
                    lapply(
                        match_idx,
                        function(i) {
                            seq(
                                max(1, i - 3),
                                min(length(x), i + 3)
                            )
                        }
                    )
                )
            )
        )

        hit_counter <- hit_counter + 1

        hits[[hit_counter]] <- list(
            file = path,
            lines = x[context_idx],
            line_numbers = context_idx
        )
    }
}

## =============================================================================
## 4. Build audit report
## =============================================================================

report <- c(
    "== S7 alternative surgery-definition source audit ==",
    paste0("Root: ", ROOT),
    "",
    paste0(
        "Candidate files by filename: ",
        length(candidate_files)
    ),
    paste0(
        "Readable text candidates: ",
        length(text_candidates)
    ),
    paste0(
        "Files with content hits: ",
        length(hits)
    ),
    ""
)

if (length(candidate_files) > 0) {

    report <- c(
        report,
        "== Candidate files ==",
        candidate_files,
        ""
    )
}

if (length(hits) > 0) {

    report <- c(
        report,
        "== Content hits with context =="
    )

    for (h in hits) {

        report <- c(
            report,
            "",
            paste0(
                "FILE: ",
                h$file
            )
        )

        for (j in seq_along(h$lines)) {

            report <- c(
                report,
                sprintf(
                    "L%-5d %s",
                    h$line_numbers[j],
                    h$lines[j]
                )
            )
        }
    }

} else {

    report <- c(
        report,
        "",
        "No alternative-surgery result text was found automatically.",
        "",
        "DO NOT reconstruct S7 from analytic_cohort_A.csv.",
        "If the old sensitivity results exist only in a Word/Excel file,",
        "place or identify that file and inspect it separately."
    )
}

## Historical checksum only — NOT used as analysis output.
report <- c(
    report,
    "",
    "== Historical checksum for identification only ==",
    "Expected population sizes previously recorded:",
    "  20–39 vs >=40: n = 7,450",
    "  Exclude code 32: n = 6,676",
    "  Colectomy broadened to 50–79: n = 7,115",
    "",
    "Expected quintile point estimates previously recorded (pp):",
    "  20–39 vs >=40:          -0.6, -2.6, -1.7, -3.6, -6.5",
    "  Exclude code 32:        -0.1, -3.6, -0.9, -5.2, -5.2",
    "  Colectomy incl 50–79:   -0.4, -2.3, -1.8, -3.5, -6.5",
    "",
    "These checksum values are for source identification only.",
    "They must not be promoted to a final supplementary table unless",
    "the underlying existing result source is located and verified."
)

writeLines(
    report,
    OUT_FILE,
    useBytes = TRUE
)

## =============================================================================
## 5. Console
## =============================================================================

cat(
    paste(
        report,
        collapse = "\n"
    ),
    "\n"
)

cat("\n============================================================\n")
cat("S7 SOURCE AUDIT COMPLETE\n")
cat("============================================================\n")
cat("Saved: ", OUT_FILE, "\n", sep = "")
cat("No statistical analysis was re-run.\n")
cat("============================================================\n")
