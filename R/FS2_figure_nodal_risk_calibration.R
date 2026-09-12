## =============================================================================
## A6 - FIGURE 2 ONLY
## Calibration of predicted nodal metastasis risk
## FINAL LOCK VERSION
##
## Root:
##   <project root>
##
## INPUT
##   02_data/analytic_cohort_A_weighted.csv
##
## OUTPUT
##   04_figures/Figure2_Nodal_Risk_Calibration_FINAL.png
##   04_figures/Figure2_Nodal_Risk_Calibration_FINAL.tiff
##
## STANDARD
##   - Color figure
##   - White background
##   - PNG 600 dpi
##   - TIFF 600 dpi + LZW compression
##   - Development curve = apparent calibration
##   - Bootstrap-corrected development performance shown in subtitle
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
DATA_DIR <- file.path(ROOT, "02_data")
FIG_DIR <- file.path(ROOT, "04_figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先运行：install.packages('ggplot2')")
}

library(ggplot2)

## =============================================================================
## 1. Locate and read frozen weighted cohort
## =============================================================================

WT_FILE <- file.path(DATA_DIR, "analytic_cohort_A_weighted.csv")
if (!file.exists(WT_FILE)) WT_FILE <- file.path(ROOT, "analytic_cohort_A_weighted.csv")
if (!file.exists(WT_FILE)) stop("找不到 analytic_cohort_A_weighted.csv")

WT <- read.csv(WT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(WT) != 6951) {
    stop("Frozen cohort fingerprint failed: N should be 6951.")
}

## =============================================================================
## 2. Define cohorts
## =============================================================================

dev <- WT[
    WT$exposure == "colectomy" &
    !is.na(WT$ln_examined) &
    WT$ln_examined >= 12 &
    !is.na(WT$n_positive),
    ,
    drop = FALSE
]

transport <- WT[
    WT$exposure == "limited" &
    !is.na(WT$ln_examined) &
    WT$ln_examined >= 12 &
    !is.na(WT$n_positive),
    ,
    drop = FALSE
]

if (nrow(dev) != 3252) stop("Development cohort fingerprint failed: expected n=3252.")
if (nrow(transport) != 879) stop("Transport cohort fingerprint failed: expected n=879.")

## =============================================================================
## 3. Helper functions
## =============================================================================

calibration_bins <- function(d, n_bins = 10) {
    cuts <- unique(as.numeric(
        quantile(d$nodal_risk, probs = 0:n_bins / n_bins, na.rm = TRUE, names = FALSE)
    ))
    if (length(cuts) < 4) stop("Insufficient unique calibration cut points.")

    cuts[1] <- cuts[1] - 1e-10
    cuts[length(cuts)] <- cuts[length(cuts)] + 1e-10

    bin <- cut(
        d$nodal_risk,
        breaks = cuts,
        include.lowest = TRUE,
        labels = FALSE
    )

    groups <- sort(unique(bin[!is.na(bin)]))

    out <- data.frame(
        bin = groups,
        predicted = NA_real_,
        observed = NA_real_,
        n = NA_integer_,
        stringsAsFactors = FALSE
    )

    for (i in seq_along(groups)) {
        idx <- bin == groups[i]
        out$predicted[i] <- mean(d$nodal_risk[idx], na.rm = TRUE)
        out$observed[i] <- mean(d$n_positive[idx], na.rm = TRUE)
        out$n[i] <- sum(idx)
    }

    out
}

auc_rank <- function(y, score) {
    ok <- !is.na(y) & !is.na(score)
    y <- y[ok]
    score <- score[ok]

    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)

    r <- rank(score, ties.method = "average")
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calibration_slope <- function(y, p) {
    ok <- !is.na(y) & !is.na(p) & p > 0 & p < 1
    fit <- glm(y[ok] ~ qlogis(p[ok]), family = binomial())
    unname(coef(fit)[2])
}

## =============================================================================
## 4. Calibration data
## =============================================================================

dev_cal <- calibration_bins(dev, 10)
dev_cal$Cohort <- "Development (apparent)"

tr_cal <- calibration_bins(transport, 10)
tr_cal$Cohort <- "Adequate-node transport"

cal_dat <- rbind(dev_cal, tr_cal)

## =============================================================================
## 5. Frozen performance fingerprints
## =============================================================================

AUC_DEV_APPARENT <- auc_rank(dev$n_positive, dev$nodal_risk)
AUC_TRANSPORT <- auc_rank(transport$n_positive, transport$nodal_risk)
SLOPE_TRANSPORT <- calibration_slope(transport$n_positive, transport$nodal_risk)

AUC_DEV_BOOT_CORRECTED <- 0.732
SLOPE_DEV_BOOT_CORRECTED <- 0.963

if (abs(AUC_DEV_APPARENT - 0.738) > 0.01) stop("Development apparent AUC fingerprint failed.")
if (abs(AUC_TRANSPORT - 0.750) > 0.02) stop("Transport AUC fingerprint failed.")
if (abs(SLOPE_TRANSPORT - 1.014) > 0.03) stop("Transport slope fingerprint failed.")

## =============================================================================
## 6. Axis range
## =============================================================================

cal_max <- max(c(cal_dat$predicted, cal_dat$observed), na.rm = TRUE)
cal_lim <- max(ceiling(cal_max * 10) / 10, 0.60)

## =============================================================================
## 7. Colors
## =============================================================================

COL_DEV <- "#2166AC"
COL_TRANSPORT <- "#D95F0E"
COL_REF <- "#666666"

## =============================================================================
## 8. Plot
## =============================================================================

p <- ggplot(
    cal_dat,
    aes(
        x = predicted,
        y = observed,
        group = Cohort,
        color = Cohort,
        shape = Cohort,
        linetype = Cohort
    )
) +
    geom_abline(
        intercept = 0,
        slope = 1,
        linewidth = 0.55,
        linetype = 3,
        color = COL_REF
    ) +
    geom_line(linewidth = 0.90) +
    geom_point(size = 3.0, stroke = 0.8) +
    scale_color_manual(
        values = c(
            "Development (apparent)" = COL_DEV,
            "Adequate-node transport" = COL_TRANSPORT
        ),
        breaks = c("Development (apparent)", "Adequate-node transport")
    ) +
    scale_shape_manual(
        values = c(
            "Development (apparent)" = 16,
            "Adequate-node transport" = 17
        ),
        breaks = c("Development (apparent)", "Adequate-node transport")
    ) +
    scale_linetype_manual(
        values = c(
            "Development (apparent)" = 1,
            "Adequate-node transport" = 2
        ),
        breaks = c("Development (apparent)", "Adequate-node transport")
    ) +
    scale_x_continuous(
        limits = c(0, cal_lim),
        breaks = seq(0, cal_lim, by = 0.1),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
        limits = c(0, cal_lim),
        breaks = seq(0, cal_lim, by = 0.1),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0, 0.02))
    ) +
    coord_equal() +
    labs(
        title = "Calibration of predicted nodal metastasis risk",
        subtitle = "Development (bootstrap-corrected): AUC 0.732, slope 0.963; transport: AUC 0.750, slope 1.014",
        x = "Mean predicted nodal risk",
        y = "Observed node-positive proportion",
        color = NULL,
        shape = NULL,
        linetype = NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
        plot.title = element_text(
            family = "Arial",
            face = "bold",
            size = 12.5,
            hjust = 0
        ),
        plot.subtitle = element_text(
            family = "Arial",
            size = 8.8,
            hjust = 0,
            margin = margin(b = 8)
        ),
        axis.title = element_text(
            family = "Arial",
            size = 10.5
        ),
        axis.text = element_text(
            family = "Arial",
            size = 9.5,
            color = "black"
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(
            family = "Arial",
            size = 9.3
        ),
        legend.key.width = unit(0.55, "inches"),
        plot.margin = margin(8, 10, 8, 8)
    )

## =============================================================================
## 9. Save
## =============================================================================

PNG_FILE <- file.path(FIG_DIR, "Figure2_Nodal_Risk_Calibration_FINAL.png")
TIFF_FILE <- file.path(FIG_DIR, "Figure2_Nodal_Risk_Calibration_FINAL.tiff")

ggsave(
    filename = PNG_FILE,
    plot = p,
    width = 7.2,
    height = 6.3,
    units = "in",
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = TIFF_FILE,
    plot = p,
    width = 7.2,
    height = 6.3,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

## =============================================================================
## 10. Console
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 2 FINAL LOCK VERSION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Development apparent AUC = %.3f\n", AUC_DEV_APPARENT))
cat(sprintf("Bootstrap-corrected development AUC = %.3f\n", AUC_DEV_BOOT_CORRECTED))
cat(sprintf("Bootstrap-corrected development slope = %.3f\n", SLOPE_DEV_BOOT_CORRECTED))
cat(sprintf("Adequate-node transport AUC = %.3f\n", AUC_TRANSPORT))
cat(sprintf("Adequate-node transport slope = %.3f\n", SLOPE_TRANSPORT))
cat("\nPNG : ", PNG_FILE, "\n", sep = "")
cat("TIFF: ", TIFF_FILE, "\n", sep = "")
cat("\nPNG = 600 dpi.\n")
cat("TIFF = 600 dpi + LZW.\n")
cat("============================================================\n")
