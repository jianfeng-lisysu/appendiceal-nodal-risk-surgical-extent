## =============================================================================
## A6 - FIGURE 4 ONLY
## Continuous association between predicted nodal risk and 5-year
## cancer-specific mortality risk difference
## FINAL LOCK VERSION
##
## Root:
##   <project root>
##
## INPUTS
##   03_results/04b_rd_curve_FINAL_B1000.csv
##   02_data/analytic_cohort_A_weighted.csv
##
## OUTPUTS
##   04_figures/Figure4_Continuous_RD_Curve_FINAL.png
##   04_figures/Figure4_Continuous_RD_Curve_FINAL.tiff
##
## STANDARD
##   - Color figure
##   - White background
##   - PNG 600 dpi
##   - TIFF 600 dpi + LZW compression
##   - Frozen continuous full-pipeline bootstrap results only
##   - RD = oncologic colectomy - limited resection
##   - Negative RD favors oncologic colectomy
##   - Pointwise 95% CI ribbon
##   - P10/P90 vertical reference lines
##   - Global interaction P = 0.2264
## =============================================================================

rm(list = ls())

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
RESULT_DIR <- file.path(ROOT, "03_results")
DATA_DIR <- file.path(ROOT, "02_data")
FIG_DIR <- file.path(ROOT, "04_figures")

dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先运行：install.packages('ggplot2')")
}

library(ggplot2)


## =============================================================================
## 1. Locate files
## =============================================================================

CURVE_FILE <- file.path(
    RESULT_DIR,
    "04b_rd_curve_FINAL_B1000.csv"
)

if (!file.exists(CURVE_FILE)) {
    CURVE_FILE <- file.path(
        ROOT,
        "04b_rd_curve_FINAL_B1000.csv"
    )
}

if (!file.exists(CURVE_FILE)) {
    stop("找不到 04b_rd_curve_FINAL_B1000.csv")
}


WT_FILE <- file.path(
    DATA_DIR,
    "analytic_cohort_A_weighted.csv"
)

if (!file.exists(WT_FILE)) {
    WT_FILE <- file.path(
        ROOT,
        "analytic_cohort_A_weighted.csv"
    )
}

if (!file.exists(WT_FILE)) {
    stop("找不到 analytic_cohort_A_weighted.csv")
}


## =============================================================================
## 2. Read frozen data
## =============================================================================

CURVE <- read.csv(
    CURVE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

WT <- read.csv(
    WT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

needed_cols <- c(
    "risk",
    "rd",
    "lo",
    "hi"
)

if (!all(needed_cols %in% names(CURVE))) {
    stop(
        "04b_rd_curve_FINAL_B1000.csv 缺少字段：",
        paste(
            setdiff(needed_cols, names(CURVE)),
            collapse = ", "
        )
    )
}

if (nrow(WT) != 6951) {
    stop(
        "Frozen cohort fingerprint failed: expected N=6951."
    )
}


## =============================================================================
## 3. Prepare plotting values
## =============================================================================

CURVE <- CURVE[
    order(CURVE$risk),
    ,
    drop = FALSE
]

CURVE$risk_pct <- CURVE$risk * 100
CURVE$rd_pp <- CURVE$rd * 100
CURVE$lo_pp <- CURVE$lo * 100
CURVE$hi_pp <- CURVE$hi * 100


P10 <- as.numeric(
    quantile(
        WT$nodal_risk,
        probs = 0.10,
        na.rm = TRUE,
        names = FALSE
    )
) * 100

P90 <- as.numeric(
    quantile(
        WT$nodal_risk,
        probs = 0.90,
        na.rm = TRUE,
        names = FALSE
    )
) * 100


## =============================================================================
## 4. Frozen sanity checks
## =============================================================================

if (abs(P10 - 5.37) > 0.5) {
    warning(
        sprintf(
            "P10 is %.3f%%; expected approximately 5.4%%.",
            P10
        )
    )
}

if (abs(P90 - 42.81) > 1.0) {
    warning(
        sprintf(
            "P90 is %.3f%%; expected approximately 42.8%%.",
            P90
        )
    )
}


## =============================================================================
## 5. Journal colors
## =============================================================================

COL_CURVE <- "#6A51A3"
COL_RIBBON <- "#D9EAF7"
COL_ZERO <- "#4D4D4D"
COL_REF <- "#777777"
COL_TEXT <- "#333333"


## =============================================================================
## 6. Axis limits
## =============================================================================

x_min <- floor(
    min(
        CURVE$risk_pct,
        na.rm = TRUE
    )
)

x_max <- ceiling(
    max(
        CURVE$risk_pct,
        na.rm = TRUE
    )
)

y_min <- floor(
    min(
        CURVE$lo_pp,
        na.rm = TRUE
    ) - 0.5
)

y_max <- ceiling(
    max(
        CURVE$hi_pp,
        na.rm = TRUE
    ) + 0.5
)


## =============================================================================
## 7. Annotation positions
##
## Put P10 label to the RIGHT of the P10 line.
## Put P90 label to the LEFT of the P90 line.
## Global-interaction annotation sits below the very top so nothing collides.
## =============================================================================

p10_label_x <- P10 + 0.7
p90_label_x <- P90 - 0.7
reference_label_y <- y_max - 0.55

interaction_x <- x_max - 1.0
interaction_y <- y_max - 1.45


## =============================================================================
## 8. Plot
## =============================================================================

p <- ggplot(
    CURVE,
    aes(
        x = risk_pct,
        y = rd_pp
    )
) +

    geom_ribbon(
        aes(
            ymin = lo_pp,
            ymax = hi_pp
        ),
        fill = COL_RIBBON,
        alpha = 0.85,
        linewidth = 0
    ) +

    geom_hline(
        yintercept = 0,
        linewidth = 0.65,
        color = COL_ZERO
    ) +

    geom_vline(
        xintercept = P10,
        linetype = 3,
        linewidth = 0.55,
        color = COL_REF
    ) +

    geom_vline(
        xintercept = P90,
        linetype = 3,
        linewidth = 0.55,
        color = COL_REF
    ) +

    geom_line(
        linewidth = 1.0,
        color = COL_CURVE
    ) +

    geom_point(
        size = 2.7,
        color = COL_CURVE
    ) +

    annotate(
        "text",
        x = p10_label_x,
        y = reference_label_y,
        label = sprintf(
            "P10 = %.1f%%",
            P10
        ),
        hjust = 0,
        vjust = 1,
        family = "Arial",
        size = 3.0,
        color = COL_REF
    ) +

    annotate(
        "text",
        x = p90_label_x,
        y = reference_label_y,
        label = sprintf(
            "P90 = %.1f%%",
            P90
        ),
        hjust = 1,
        vjust = 1,
        family = "Arial",
        size = 3.0,
        color = COL_REF
    ) +

    annotate(
        "label",
        x = interaction_x,
        y = interaction_y,
        label = "Global interaction P = 0.226",
        hjust = 1,
        vjust = 1,
        family = "Arial",
        size = 3.05,
        color = COL_TEXT,
        fill = "white",
        label.size = 0.25,
        label.padding = grid::unit(
            0.12,
            "lines"
        )
    ) +

    scale_x_continuous(
        breaks = seq(
            0,
            50,
            by = 10
        ),
        limits = c(
            x_min,
            x_max
        ),
        expand = expansion(
            mult = c(
                0.01,
                0.01
            )
        )
    ) +

    scale_y_continuous(
        breaks = seq(
            -15,
            10,
            by = 5
        ),
        limits = c(
            y_min,
            y_max
        ),
        expand = expansion(
            mult = c(
                0.01,
                0.01
            )
        )
    ) +

    labs(
        title = "Predicted nodal risk and 5-year cancer-specific mortality risk difference",
        x = "Predicted nodal metastasis risk, %",
        y = paste0(
            "Absolute 5-year cancer-specific mortality risk difference, pp\n",
            "(oncologic colectomy - limited resection)"
        )
    ) +

    theme_classic(
        base_family = "Arial",
        base_size = 11
    ) +

    theme(
        plot.title = element_text(
            family = "Arial",
            face = "bold",
            size = 12.5,
            hjust = 0
        ),
        axis.title.x = element_text(
            family = "Arial",
            size = 10.5,
            margin = margin(
                t = 8
            )
        ),
        axis.title.y = element_text(
            family = "Arial",
            size = 10.5,
            margin = margin(
                r = 8
            )
        ),
        axis.text = element_text(
            family = "Arial",
            size = 9.5,
            color = "black"
        ),
        plot.margin = margin(
            10,
            12,
            8,
            8
        )
    )


## =============================================================================
## 9. Save true 600-dpi PNG + TIFF
## =============================================================================

PNG_FILE <- file.path(
    FIG_DIR,
    "Figure4_Continuous_RD_Curve_FINAL.png"
)

TIFF_FILE <- file.path(
    FIG_DIR,
    "Figure4_Continuous_RD_Curve_FINAL.tiff"
)


ggsave(
    filename = PNG_FILE,
    plot = p,
    width = 8.4,
    height = 5.9,
    units = "in",
    dpi = 600,
    bg = "white"
)


ggsave(
    filename = TIFF_FILE,
    plot = p,
    width = 8.4,
    height = 5.9,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)


## =============================================================================
## 10. Console audit
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 4 FINAL LOCK VERSION COMPLETE\n")
cat("============================================================\n")

cat(
    sprintf(
        "P10 = %.3f%%\n",
        P10
    )
)

cat(
    sprintf(
        "P90 = %.3f%%\n",
        P90
    )
)

cat(
    "Global interaction P = 0.2264\n"
)

cat("\nCurve grid:\n")

for (i in seq_len(nrow(CURVE))) {
    cat(
        sprintf(
            "Risk %5.1f%% | RD %+.2f pp | 95%% CI %+.2f to %+.2f\n",
            CURVE$risk_pct[i],
            CURVE$rd_pp[i],
            CURVE$lo_pp[i],
            CURVE$hi_pp[i]
        )
    )
}

cat("\nPNG : ", PNG_FILE, "\n", sep = "")
cat("TIFF: ", TIFF_FILE, "\n", sep = "")

cat("\nPNG = 600 dpi.\n")
cat("TIFF = 600 dpi + LZW.\n")
cat("============================================================\n")
