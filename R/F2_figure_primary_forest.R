## =============================================================================
## A6 - FIGURE 3 ONLY
## Five-year cancer-specific mortality differences by predicted nodal risk
## JOURNAL-READY FINAL CANDIDATE
##
## Root:
##   <project root>
##
## INPUT
##   03_results/04_curve_data_FINAL_B1000.csv
##
## OUTPUT
##   04_figures/Figure3_Primary_RD_Forest_FINAL.png
##   04_figures/Figure3_Primary_RD_Forest_FINAL.tiff
##
## STANDARD
##   - Color figure
##   - White background
##   - PNG 600 dpi
##   - TIFF 600 dpi + LZW compression
##   - Frozen B=1000 estimates only
##   - RD = oncologic colectomy - limited resection
##   - Negative RD favors oncologic colectomy
##   - Red dashed line = prespecified -3 pp clinically important margin
## =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("A6_ROOT", unset = ".")  # set A6_ROOT to your project directory
RESULT_DIR <- file.path(ROOT, "03_results")
FIG_DIR <- file.path(ROOT, "04_figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先运行：install.packages('ggplot2')")
}

library(ggplot2)

## =============================================================================
## 1. Locate and read frozen result
## =============================================================================

INPUT_FILE <- file.path(RESULT_DIR, "04_curve_data_FINAL_B1000.csv")

if (!file.exists(INPUT_FILE)) {
    INPUT_FILE <- file.path(ROOT, "04_curve_data_FINAL_B1000.csv")
}

if (!file.exists(INPUT_FILE)) {
    stop("找不到 04_curve_data_FINAL_B1000.csv")
}

d <- read.csv(
    INPUT_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

needed <- c("stratum", "rd", "lo", "hi")

if (!all(needed %in% names(d))) {
    stop(
        "04_curve_data_FINAL_B1000.csv 缺少必要字段：",
        paste(setdiff(needed, names(d)), collapse = ", ")
    )
}

## =============================================================================
## 2. Frozen fingerprint
## =============================================================================

expected_order <- c(
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "bottom_decile",
    "risk_lt_5pct"
)

if (!all(expected_order %in% d$stratum)) {
    stop("主分析7个冻结strata不完整。")
}

d <- d[
    match(expected_order, d$stratum),
    ,
    drop = FALSE
]

expected_point_pp <- c(
    0.3,
    -3.6,
    -0.7,
    -5.1,
    -6.0,
    3.1,
    4.2
)

observed_point_pp <- d$rd * 100

if (
    max(
        abs(
            observed_point_pp -
            expected_point_pp
        )
    ) >
    0.15
) {
    stop("Figure 3 frozen point-estimate fingerprint failed.")
}

## =============================================================================
## 3. Prepare display data
## =============================================================================

label_map <- c(
    Q1 = "Q1 (lowest risk)",
    Q2 = "Q2",
    Q3 = "Q3",
    Q4 = "Q4",
    Q5 = "Q5 (highest risk)",
    bottom_decile = "Bottom predicted-risk decile",
    risk_lt_5pct = "Predicted nodal risk <5%"
)

d$label <- unname(
    label_map[
        d$stratum
    ]
)

## Numeric y positions create a deliberate gap between Q5 and the two anchors.
d$y <- c(
    8,
    7,
    6,
    5,
    4,
    2.5,
    1.5
)

d$rd_pp <- d$rd * 100
d$lo_pp <- d$lo * 100
d$hi_pp <- d$hi * 100

d$analysis_group <- ifelse(
    d$stratum %in% c(
        "bottom_decile",
        "risk_lt_5pct"
    ),
    "Prespecified extreme low-risk anchors",
    "Risk quintiles"
)

d$estimate_text <- sprintf(
    "%+.1f (%+.1f to %+.1f)",
    d$rd_pp,
    d$lo_pp,
    d$hi_pp
)

## Prespecified 3-pp criterion only applies to Q1 and two low-risk anchors.
d$criterion <- ""

criterion_rows <- d$stratum %in% c(
    "Q1",
    "bottom_decile",
    "risk_lt_5pct"
)

d$criterion[criterion_rows] <- ifelse(
    d$lo_pp[criterion_rows] > -3,
    "Met",
    "Not met"
)

## =============================================================================
## 4. Journal colors
## =============================================================================

COL_QUINTILE <- "#2166AC"
COL_ANCHOR <- "#D95F0E"
COL_ZERO <- "#444444"
COL_MARGIN <- "#B2182B"
COL_SEPARATOR <- "#BDBDBD"

## =============================================================================
## 5. Figure geometry
## =============================================================================

x_left <- floor(
    min(
        d$lo_pp,
        na.rm = TRUE
    ) -
    1
)

estimate_x <- ceiling(
    max(
        d$hi_pp,
        na.rm = TRUE
    ) +
    1.2
)

criterion_x <- estimate_x + 6.3

x_right <- criterion_x + 2.4

## =============================================================================
## 6. Plot
## =============================================================================

p <- ggplot(
    d,
    aes(
        x = rd_pp,
        y = y
    )
) +

    geom_vline(
        xintercept = 0,
        linewidth = 0.65,
        color = COL_ZERO
    ) +

    geom_vline(
        xintercept = -3,
        linewidth = 0.65,
        linetype = 2,
        color = COL_MARGIN
    ) +

    annotate(
        "segment",
        x = x_left,
        xend = x_right,
        y = 3.25,
        yend = 3.25,
        linewidth = 0.45,
        color = COL_SEPARATOR
    ) +

    geom_errorbar(
        aes(
            xmin = lo_pp,
            xmax = hi_pp,
            color = analysis_group
        ),
        width = 0.16,
        linewidth = 0.78
    ) +

    geom_point(
        aes(
            color = analysis_group
        ),
        size = 3.0
    ) +

    geom_text(
        aes(
            x = estimate_x,
            label = estimate_text
        ),
        hjust = 0,
        family = "Arial",
        size = 3.15,
        color = "black"
    ) +

    geom_text(
        data = d[
            criterion_rows,
            ,
            drop = FALSE
        ],
        aes(
            x = criterion_x,
            label = criterion
        ),
        hjust = 0,
        family = "Arial",
        size = 3.15,
        color = "black"
    ) +

    annotate(
        "text",
        x = estimate_x,
        y = 8.75,
        label = "RD (95% CI), pp",
        hjust = 0,
        family = "Arial",
        fontface = "bold",
        size = 3.25
    ) +

    annotate(
        "text",
        x = criterion_x,
        y = 8.75,
        label = "3-pp criterion",
        hjust = 0,
        family = "Arial",
        fontface = "bold",
        size = 3.25
    ) +

    annotate(
        "text",
        x = -3.15,
        y = 8.75,
        label = "Prespecified margin: -3 pp",
        hjust = 1,
        family = "Arial",
        size = 3.0,
        color = COL_MARGIN
    ) +

    scale_y_continuous(
        breaks = d$y,
        labels = d$label,
        limits = c(
            1.0,
            9.0
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_x_continuous(
        breaks = seq(
            -15,
            5,
            by = 5
        ),
        limits = c(
            x_left,
            x_right
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_color_manual(
        values = c(
            "Risk quintiles" = COL_QUINTILE,
            "Prespecified extreme low-risk anchors" = COL_ANCHOR
        ),
        breaks = c(
            "Risk quintiles",
            "Prespecified extreme low-risk anchors"
        )
    ) +

    labs(
        title = "Five-year cancer-specific mortality differences by predicted nodal risk",
        x = "Absolute risk difference, percentage points\n(oncologic colectomy - limited resection)",
        y = NULL,
        color = NULL
    ) +

    coord_cartesian(
        clip = "off"
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
        axis.text.x = element_text(
            family = "Arial",
            size = 9.5,
            color = "black"
        ),
        axis.text.y = element_text(
            family = "Arial",
            size = 9.7,
            color = "black"
        ),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(
            family = "Arial",
            size = 9.2
        ),
        legend.key.width = grid::unit(
            0.45,
            "inches"
        ),
        plot.margin = margin(
            12,
            20,
            8,
            8
        )
    )

## =============================================================================
## 7. Save true 600-dpi PNG + TIFF
## =============================================================================

PNG_FILE <- file.path(
    FIG_DIR,
    "Figure3_Primary_RD_Forest_FINAL.png"
)

TIFF_FILE <- file.path(
    FIG_DIR,
    "Figure3_Primary_RD_Forest_FINAL.tiff"
)

ggsave(
    filename = PNG_FILE,
    plot = p,
    width = 10.5,
    height = 6.4,
    units = "in",
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = TIFF_FILE,
    plot = p,
    width = 10.5,
    height = 6.4,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

## =============================================================================
## 8. Console
## =============================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 3 COMPLETE\n")
cat("============================================================\n")

for (i in seq_len(nrow(d))) {
    cat(
        sprintf(
            "%-30s %+.1f (%+.1f to %+.1f) pp %s\n",
            d$label[i],
            d$rd_pp[i],
            d$lo_pp[i],
            d$hi_pp[i],
            ifelse(
                nzchar(d$criterion[i]),
                paste0(
                    "| 3-pp criterion: ",
                    d$criterion[i]
                ),
                ""
            )
        )
    )
}

cat("\nPNG : ", PNG_FILE, "\n", sep = "")
cat("TIFF: ", TIFF_FILE, "\n", sep = "")

cat("\nPNG = 600 dpi.\n")
cat("TIFF = 600 dpi + LZW.\n")
cat("============================================================\n")
