## =============================================================================
## A6 - SUPPLEMENTARY TABLE S6 ONLY
## Histology- and risk-stratified 5-year cancer-specific mortality differences — COMPACT FINAL LOCK
## after subgroup-specific overlap weighting
##
## Root:
##   <project root>
##
## INPUT
##   03_results/06_histology_subgroup_reweight_FINAL_B1000.csv
##
## OUTPUTS
##   05_tables/Supplementary_Table_S6_Histology_Risk_Subgroups_FINAL.docx
##   05_tables/Supplementary_Table_S6_Histology_Risk_Subgroups_FINAL.csv
##
## FINAL B=1000 RESULTS ALREADY COMPLETED
##   - Histology: nonmucinous, mucinous, GCA, SRCC
##   - Risk half: below vs at/above overall median predicted nodal risk
##   - Subgroup-specific propensity score and overlap weighting
##   - Full-pipeline bootstrap B=1000
##
## IMPORTANT
##   This script formats the already-completed final sensitivity analysis.
##   It does NOT re-run bootstrap or alter any estimate.
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
RESULT_DIR <- file.path(ROOT, "03_results")
TABLE_DIR <- file.path(ROOT, "05_tables")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 1. Packages
## =============================================================================

need_pkg <- c("officer", "flextable")

miss_pkg <- need_pkg[
    !vapply(
        need_pkg,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(miss_pkg) > 0) {
    stop(
        "请先安装：",
        paste(miss_pkg, collapse = ", ")
    )
}

library(officer)
library(flextable)

## =============================================================================
## 2. Locate final subgroup file
## =============================================================================

INPUT_FILE <- file.path(
    RESULT_DIR,
    "06_histology_subgroup_reweight_FINAL_B1000.csv"
)

if (!file.exists(INPUT_FILE)) {
    INPUT_FILE <- file.path(
        ROOT,
        "06_histology_subgroup_reweight_FINAL_B1000.csv"
    )
}

if (!file.exists(INPUT_FILE)) {
    stop(
        "找不到 06_histology_subgroup_reweight_FINAL_B1000.csv"
    )
}

d <- read.csv(
    INPUT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

## =============================================================================
## 3. Required columns
## =============================================================================

needed <- c(
    "subgroup",
    "histology",
    "risk_half",
    "n",
    "limited_n",
    "colectomy_n",
    "old_global_weight_rd",
    "subgroup_reweighted_rd",
    "rd_lo",
    "rd_hi",
    "bootstrap_effective_B",
    "ess_limited",
    "ess_colectomy",
    "max_abs_smd_all",
    "n_smd_gt_010_all",
    "status"
)

miss <- setdiff(needed, names(d))

if (length(miss) > 0) {
    stop(
        "最终 subgroup CSV 缺少字段：",
        paste(miss, collapse = ", ")
    )
}

## =============================================================================
## 4. Freeze expected order
## =============================================================================

expected_order <- c(
    "nonmucinous_low",
    "nonmucinous_high",
    "mucinous_low",
    "mucinous_high",
    "GCA_low",
    "GCA_high",
    "SRCC_low",
    "SRCC_high"
)

if (!all(expected_order %in% d$subgroup)) {
    stop("8个最终 histology × risk-half subgroup 不完整。")
}

d <- d[
    match(expected_order, d$subgroup),
    ,
    drop = FALSE
]

## =============================================================================
## 5. Hard fingerprints for final B=1000 point estimates
## =============================================================================

expected_new_rd_pp <- c(
    -2.9,
    -8.2,
    -2.3,
    -3.6,
    -1.0,
    +1.5,
    NA,
    +3.3
)

expected_lo_pp <- c(
    -9.5,
    -14.1,
    -7.2,
    -11.0,
    -4.0,
    -5.1,
    NA,
    -9.9
)

expected_hi_pp <- c(
    +2.4,
    -2.5,
    +2.4,
    +3.5,
    +1.6,
    +8.6,
    NA,
    +16.5
)

for (i in seq_len(nrow(d))) {

    if (is.finite(expected_new_rd_pp[i])) {

        if (
            !is.finite(d$subgroup_reweighted_rd[i]) ||
            abs(d$subgroup_reweighted_rd[i] * 100 - expected_new_rd_pp[i]) > 0.15
        ) {
            stop(
                "Final subgroup RD fingerprint failed for ",
                d$subgroup[i]
            )
        }

        if (
            !is.finite(d$rd_lo[i]) ||
            abs(d$rd_lo[i] * 100 - expected_lo_pp[i]) > 0.15
        ) {
            stop(
                "Final subgroup lower CI fingerprint failed for ",
                d$subgroup[i]
            )
        }

        if (
            !is.finite(d$rd_hi[i]) ||
            abs(d$rd_hi[i] * 100 - expected_hi_pp[i]) > 0.15
        ) {
            stop(
                "Final subgroup upper CI fingerprint failed for ",
                d$subgroup[i]
            )
        }
    }
}

## =============================================================================
## 6. Publication labels
## =============================================================================

hist_label <- c(
    nonmucinous = "Nonmucinous adenocarcinoma",
    mucinous = "Mucinous adenocarcinoma",
    GCA = "Goblet cell adenocarcinoma",
    SRCC = "Signet-ring cell carcinoma"
)

risk_label <- c(
    low = "Below median risk",
    high = "At/above median risk"
)

## =============================================================================
## 7. Formatting helpers
## =============================================================================

fmt_pp <- function(x) {
    if (!is.finite(x)) return("—")
    sprintf("%+.1f", x * 100)
}

fmt_ci_pp <- function(lo, hi) {
    if (!is.finite(lo) || !is.finite(hi)) return("—")
    sprintf("%+.1f to %+.1f", lo * 100, hi * 100)
}

fmt_ess <- function(x) {
    if (!is.finite(x)) return("—")
    sprintf("%.1f", x)
}

fmt_smd <- function(x) {
    if (!is.finite(x)) return("—")
    if (abs(x) < 0.0005) x <- 0
    sprintf("%.3f", x)
}

fmt_b <- function(x) {
    if (!is.finite(x)) return("—")
    format(as.integer(x), big.mark = ",", scientific = FALSE)
}

## =============================================================================
## 8. Build final display table — compact journal version
## =============================================================================

estimate_status <- ifelse(
    grepl("^FAILED", d$status),
    "Insufficient",
    "Estimable"
)

table_s6 <- data.frame(
    `Histology` = unname(hist_label[d$histology]),
    `Predicted nodal-risk half` = unname(risk_label[d$risk_half]),
    `n (limited / colectomy)` = paste0(
        d$n,
        " (",
        d$limited_n,
        " / ",
        d$colectomy_n,
        ")"
    ),
    `Subgroup-OW RD, pp` = vapply(
        d$subgroup_reweighted_rd,
        fmt_pp,
        character(1)
    ),
    `95% CI, pp` = mapply(
        fmt_ci_pp,
        d$rd_lo,
        d$rd_hi,
        USE.NAMES = FALSE
    ),
    `ESS limited / colectomy` = ifelse(
        is.finite(d$ess_limited) & is.finite(d$ess_colectomy),
        paste0(
            sprintf("%.1f", d$ess_limited),
            " / ",
            sprintf("%.1f", d$ess_colectomy)
        ),
        "—"
    ),
    `Maximum weighted |SMD|` = vapply(
        d$max_abs_smd_all,
        fmt_smd,
        character(1)
    ),
    `Effective bootstrap B` = vapply(
        d$bootstrap_effective_B,
        fmt_b,
        character(1)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

## Sparse SRCC-low subgroup: report insufficiency explicitly in estimate columns.
sparse_row <- d$subgroup == "SRCC_low"

table_s6[sparse_row, "Subgroup-OW RD, pp"] <- "Insufficient"
table_s6[sparse_row, "95% CI, pp"] <- "—"
table_s6[sparse_row, "ESS limited / colectomy"] <- "—"
table_s6[sparse_row, "Maximum weighted |SMD|"] <- "—"
table_s6[sparse_row, "Effective bootstrap B"] <- "—"

## =============================================================================
## 9. Save CSV audit copy
## =============================================================================

CSV_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S6_Histology_Risk_Subgroups_FINAL.csv"
)

write.csv(
    table_s6,
    CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
)

## =============================================================================
## 10. Build journal-ready three-line table
## =============================================================================

ft <- flextable(table_s6)

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
    j = c(1, 2),
    align = "left",
    part = "all"
)

ft <- align(
    ft,
    j = 3:ncol(table_s6),
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

## Small separation at the start of each histology block.
ft <- padding(
    ft,
    i = c(1, 3, 5, 7),
    padding.top = 2.4,
    part = "body"
)

## Width-safe layout.
ft <- width(ft, j = 1, width = 2.05)
ft <- width(ft, j = 2, width = 1.45)
ft <- width(ft, j = 3, width = 1.35)
ft <- width(ft, j = 4, width = 1.15)
ft <- width(ft, j = 5, width = 1.25)
ft <- width(ft, j = 6, width = 1.45)
ft <- width(ft, j = 7, width = 1.20)
ft <- width(ft, j = 8, width = 1.10)

ft <- set_table_properties(
    ft,
    layout = "fixed",
    width = 1
)

## =============================================================================
## 11. Word document
## =============================================================================

DOCX_FILE <- file.path(
    TABLE_DIR,
    "Supplementary_Table_S6_Histology_Risk_Subgroups_FINAL.docx"
)

doc <- read_docx()

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            "Supplementary Table S6. Histology- and risk-stratified 5-year cancer-specific mortality differences after subgroup-specific overlap weighting",
            fp_text(
                font.family = "Times New Roman",
                font.size = 9.2,
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
    "RD is the absolute 5-year cancer-specific mortality risk difference, defined as oncologic colectomy minus limited resection; ",
    "negative values indicate lower cancer-specific mortality associated with oncologic colectomy. ",
    "Risk halves were defined by the overall median predicted nodal risk, which was re-estimated within each full-pipeline bootstrap replicate. ",
    "Within each histology-by-risk subgroup, the propensity-score model was refitted using the prespecified baseline covariates except histology, ",
    "which was constant by design, and subgroup-specific overlap weights were then applied. ",
    "Subgroup-specific overlap-weighted estimates are shown because these were the preferred sensitivity estimates; ",
    "the corresponding primary cohort-wide overlap-weighted estimates are available in the analysis audit file. ",
    "Balance was considered acceptable when all modeled covariates had absolute standardized mean difference <0.10; ",
    "all estimable subgroups met this criterion, so only the maximum weighted |SMD| is displayed. ",
    "The signet-ring cell carcinoma low-risk subgroup was too sparse for stable estimation; ",
    "the high-risk signet-ring cell carcinoma estimate remains exploratory because of limited sample size."
)

doc <- body_add_fpar(
    doc,
    fpar(
        ftext(
            note_text,
            fp_text(
                font.family = "Times New Roman",
                font.size = 7.2
            )
        )
    )
)

doc <- body_end_section_landscape(doc)

print(
    doc,
    target = DOCX_FILE
)

## =============================================================================
## 12. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY TABLE S6 COMPLETE\n")
cat("============================================================\n")

for (i in seq_len(nrow(d))) {

    if (d$subgroup[i] == "SRCC_low") {
        cat(
            sprintf(
                "%-20s | n=%d | insufficient for stable estimation\n",
                d$subgroup[i],
                d$n[i]
            )
        )
    } else {
        cat(
            sprintf(
                "%-20s | subgroup-OW RD %+.1f pp (%+.1f to %+.1f) | max|SMD| %.3f | B=%d\n",
                d$subgroup[i],
                d$subgroup_reweighted_rd[i] * 100,
                d$rd_lo[i] * 100,
                d$rd_hi[i] * 100,
                d$max_abs_smd_all[i],
                d$bootstrap_effective_B[i]
            )
        )
    }
}

cat("\nDOCX: ", DOCX_FILE, "\n", sep = "")
cat("CSV : ", CSV_FILE, "\n", sep = "")

cat("\n")
cat("Formatting only; final B=1000 subgroup estimates were not re-estimated.\n")
cat("Standard three-line table; no vertical rules.\n")
cat("============================================================\n")
