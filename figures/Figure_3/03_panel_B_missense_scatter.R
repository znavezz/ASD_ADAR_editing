root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Combined Figure — Panel B  (missense optimization scatter).
#
# x = variant SIFT (pre-edit, deleterious < 0.05)
# y = post-edit SIFT via best neighbor-edit  (tolerated >= 0.05)
# Log10 on BOTH axes with a small pseudo-log floor for x = 0 rows.
# Style copied from ver-1 Figure 2 (no marginals): ocean-teal pass points,
# grey "other" points, dashed 0.05 threshold cross-hairs, quadrant n(%),
# arrow strips ("more severe" bottom, "less deleterious" left).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales); library(patchwork); library(tibble)
})

VER1     <- asd_source_v1()
VER2     <- asd_source()
INT_DIR <- asd_intermediate_dir()
OUT_DIR <- asd_figure_dir(2)
OUT_MAIN <- asd_panels_dir()
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)
source(asd_code("ver1", "_helpers.R"))   # theme_publication

# --- reload the FULL missense set (no SIFT<0.05 / improves filters) --------
paneldata <- readRDS(file.path(INT_DIR, "MissenseOptimization_full_flat.rds")) %>%
  filter(!is.na(variant_sift), !is.na(post_sift)) %>%
  mutate(hl = ifelse(!sift_tolerated(variant_sift) & sift_tolerated(post_sift),
                     "pass", "other"))

n_total <- nrow(paneldata)
n_pass  <- sum(paneldata$hl == "pass")
pct     <- 100 * n_pass / n_total
cat("plotted:", n_total, " pass (del→tol):", n_pass,
    sprintf(" (%.1f%%)", pct), "\n", sep = "")

# --- axis transform: pseudo-log10 (handles variant_sift == 0) ---------------
# sigma controls where linearity switches to log; use 0.001 so anything
# above 1e-3 sits on the log ramp.  Ticks placed at 0.001/0.01/0.05/0.1/0.5/1.
PSEUDO_SIG <- 1e-3
axis_trans  <- scales::pseudo_log_trans(sigma = PSEUDO_SIG, base = 10)
axis_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 0.5, 1)
axis_labels <- c("0", "0.001", "0.01", "0.05", "0.1", "0.5", "1")
axis_lims   <- c(0, 1)

# --- Colours (ver-1 palette) ------------------------------------------------
OCEAN_TEAL <- "#0d6e8c"
GREY_OFF   <- "grey75"

# --- Quadrant counts + percentages ------------------------------------------
qd <- paneldata %>%
  mutate(qx = ifelse(!sift_tolerated(variant_sift), "lo", "hi"),
         qy = ifelse(sift_tolerated(post_sift),    "hi", "lo")) %>%
  count(qx, qy, name = "n") %>%
  mutate(pct = 100 * n / sum(n),
         lbl = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","), pct))
get_lbl <- function(qx, qy) {
  v <- qd$lbl[qd$qx == qx & qd$qy == qy]
  if (length(v) == 0) "0\n(0.0%)" else v
}

# Anchor positions in original (non-transformed) coords.  Because pseudo-log
# is symmetric around 0, we can place labels close to 0 on the variant side.
x_lo <- 0.005     # inside deleterious column (variant_sift < 0.05)
x_hi <- 0.30      # inside "tolerated" column (unused — all variants < 0.05)
y_lo <- 0.005     # below post_sift = 0.05  (unused — all >= 0.05)
y_hi <- 0.30      # above 0.05

# --- Main scatter -----------------------------------------------------------
p_main <- ggplot(paneldata,
                 aes(x = variant_sift, y = post_sift, colour = hl)) +
  # y = x (no-change) diagonal — above the line the rescue lifts SIFT,
  # below the line the neighbor edit lowers it.
  geom_abline(intercept = 0, slope = 1,
              linetype = "dotted", colour = "grey35", linewidth = 0.5) +
  geom_hline(yintercept = SIFT_CUT, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_vline(xintercept = SIFT_CUT, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_point(alpha = 0.65, size = 1.9) +
  scale_colour_manual(values = c("pass" = OCEAN_TEAL, "other" = GREY_OFF),
                      guide = "none") +
  annotate("label", x = x_lo, y = y_hi, label = get_lbl("lo", "hi"),
           family = "Nimbus Sans", fontface = "bold", size = 4.4,
           colour = OCEAN_TEAL, fill = "white",
           label.size = 0.3, alpha = 0.85,
           label.padding = unit(0.30, "lines"), hjust = 0) +
  annotate("label", x = x_hi, y = y_hi, label = get_lbl("hi", "hi"),
           family = "Nimbus Sans", fontface = "bold", size = 4.0,
           colour = "grey25", fill = "white",
           label.size = 0, alpha = 0.65,
           label.padding = unit(0.25, "lines")) +
  annotate("label", x = x_lo, y = y_lo, label = get_lbl("lo", "lo"),
           family = "Nimbus Sans", fontface = "bold", size = 4.0,
           colour = "grey25", fill = "white",
           label.size = 0, alpha = 0.65,
           label.padding = unit(0.25, "lines"), hjust = 0) +
  annotate("label", x = x_hi, y = y_lo, label = get_lbl("hi", "lo"),
           family = "Nimbus Sans", fontface = "bold", size = 4.0,
           colour = "grey25", fill = "white",
           label.size = 0, alpha = 0.65,
           label.padding = unit(0.25, "lines")) +
  scale_x_continuous(trans  = axis_trans,
                     breaks = axis_breaks,
                     labels = axis_labels,
                     limits = axis_lims,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(trans  = axis_trans,
                     breaks = axis_breaks,
                     labels = axis_labels,
                     limits = axis_lims,
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Variant SIFT (pre-edit)",
       y = "Optimization SIFT (post-edit)",
       tag = "B") +
  theme_publication(base_size = 12, x_angle = 0, y_bold = TRUE) +
  theme(axis.title.x       = element_text(face = "bold", colour = "black"),
        plot.tag           = element_text(face = "bold", size = 16,
                                          family = "Nimbus Sans"),
        plot.tag.position  = c(0, 1),
        plot.margin        = margin(8, 2, 2, 2))

# --- Bottom arrow strip: single arrow pointing LEFT to "deleterious" -------
arrow_x_strip <- ggplot() +
  annotate("segment",
           x = 0.95, xend = 0.05, y = 0.5, yend = 0.5,
           arrow    = arrow(length = unit(0.35, "cm"), type = "closed"),
           linewidth = 0.7, colour = "grey20") +
  annotate("label",
           x = 0.5, y = 0.5, label = "deleterious",
           family = "Nimbus Sans", fontface = "italic",
           colour = "grey20", size = 4,
           fill = "white", label.size = 0,
           label.padding = unit(0.15, "lines")) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
  theme_void() +
  theme(plot.margin = margin(2, 0, 2, 0))

arrow_y_strip <- ggplot() +
  annotate("segment",
           x = 0.5, xend = 0.5, y = 0.05, yend = 0.95,
           arrow    = arrow(length = unit(0.35, "cm"), type = "closed"),
           linewidth = 0.7, colour = "grey20") +
  annotate("label",
           x = 0.5, y = 0.5, label = "tolerated",
           family = "Nimbus Sans", fontface = "italic",
           colour = "grey20", size = 4, angle = 90,
           fill = "white", label.size = 0,
           label.padding = unit(0.15, "lines")) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
  theme_void() +
  theme(plot.margin = margin(0, 2, 0, 0))

# --- Compose ---------------------------------------------------------------
spacer <- patchwork::plot_spacer()
panel <- (arrow_y_strip + p_main +
          spacer        + arrow_x_strip) +
  plot_layout(ncol = 2,
              widths  = c(0.10, 5),
              heights = c(5, 0.30))

# --- Save ------------------------------------------------------------------
out_pdf <- file.path(OUT_MAIN, "Figure_B_missense_scatter_v7.pdf")
out_png <- file.path(OUT_MAIN, "Figure_B_missense_scatter_v7.png")
out_rds <- file.path(OUT_MAIN, "Figure_B_missense_scatter_v7.rds")
ggsave(out_pdf, panel, width = 8, height = 6,
       device = cairo_pdf, bg = "transparent")
ggsave(out_png, panel, width = 8, height = 6,
       dpi = 150, bg = "white")
saveRDS(panel, out_rds)
cat("wrote:\n  ", out_pdf, "\n  ", out_png, "\n  ", out_rds, "\n", sep = "")
cat("\nplotted:", n_total,
    " del→tol pass:", n_pass,
    sprintf(" (%.1f%%)", pct), "\n")
