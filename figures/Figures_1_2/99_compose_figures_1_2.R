root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Figure 1 (ver-2) — composite v2, now with A illustration on top.
#
# Layout:
#   ┌─────────────────────────────┐
#   │             A               │  ADAR-mechanism illustration
#   ├─────┬─────┬─────┤
#   │     │  C  │  E  │
#   │  B  ├─────┼─────┤
#   │     │  D  │  F  │
#   └─────┴─────┴─────┘
#
# Sources:
#   A — Results_Plots_Tables/Figure1A_illustration.png   (external PNG)
#   B — ver-1 funnel
#   C — ver-2 ClinVar donut (project palette)
#   D — ver-1 G>A-only bar
#   E — ver-2 DNA off-targets gradient (project palette)
#   F — ver-2 bystander CADD25 gradient (project palette)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

VER1_INT <- asd_source_v1("Intermediate_tables/Figure_1")
VER2_INT <- asd_source("Intermediate_tables/Figure_1")
OUT_DIR  <- asd_figure_dir(1)
ILL_PNG  <- asd_code("Figure1A_illustration.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Panel A: raster the illustration, preserving its aspect ratio ---------
# We deliberately let the A-row container be WIDER than needed (equal to the
# full composite width) and then centre the image at native aspect inside
# that container.  Downstream, plot_layout heights control the A-row height.
img      <- readPNG(ILL_PNG)
img_h_px <- dim(img)[1]
img_w_px <- dim(img)[2]
img_asp  <- img_h_px / img_w_px    # height / width
cat(sprintf("illustration: %d × %d px, aspect %.3f\n",
            img_w_px, img_h_px, img_asp))
img_grob <- rasterGrob(img, interpolate = TRUE)   # aspect-preserving default

# Placeholders that get resolved in the "compose" block below — we need the
# container aspect (H_A/W) to know how to inset the image.
build_p_A <- function(W_container, H_container) {
  container_asp <- H_container / W_container
  if (img_asp > container_asp) {
    img_w_frac <- container_asp / img_asp
    img_h_frac <- 1
  } else {
    img_w_frac <- 1
    img_h_frac <- img_asp / container_asp
  }
  xmin <- (1 - img_w_frac) / 2
  xmax <- 1 - xmin
  ymin <- (1 - img_h_frac) / 2
  ymax <- 1 - ymin

  ggplot() +
    annotation_custom(img_grob,
                      xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    # No tag: a figure of one panel carries no letter.
    theme_void() +
    theme(plot.tag          = element_text(face = "bold", size = 18,
                                           family = "Nimbus Sans"),
          plot.tag.position = c(0, 1),
          plot.margin       = margin(4, 4, 4, 4))
}

# --- Panels B..F: existing ver-2 outputs -----------------------------------
p_B <- readRDS(asd_panel("Figure_1A_funnel_ggplot.rds", VER1_INT))
p_C <- readRDS(asd_panel("Figure_1B_clinvar_donut_ver2_ggplot.rds", VER2_INT))
p_D <- readRDS(asd_panel("Figure_1B_mismatch_GA_only_ggplot.rds", VER1_INT))
p_E <- readRDS(asd_panel("Figure_1D_offtargets_dna_gradient_ver2_ggplot.rds", VER2_INT))
p_F <- readRDS(asd_panel("Figure_1E_bystander_cadd25_gradient_ver2_ggplot.rds", VER2_INT))

p_B <- p_B + labs(tag = "A")   # funnel
p_C <- p_C + labs(tag = "B")   # ClinVar donut
p_D <- p_D + labs(tag = "C")   # consequence bar
p_E <- p_E + labs(tag = "D")   # off-targets
p_F <- p_F + labs(tag = "E")   # bystander CADD

ef_col <- (p_E / p_F)
cd_col <- (p_C / p_D) + plot_layout(heights = c(0.72, 0.28))

# Squeeze the funnel column (B) narrower + strip its right-side plot margin so
# it sits visually closer to C (the ClinVar donut column).
p_B <- p_B + theme(plot.margin = margin(5, 0, 5, 8))
p_C <- p_C + theme(plot.margin = margin(5, 8, 0, 0))

bottom_block <- (p_B | cd_col | ef_col) +
  plot_layout(widths = c(1.4, 3.0, 3.0))   # was 2.0 — B closer to C

# Full composite: illustration on top, panels below.
# 2026-07-08 update: shrink A row to ~2/3 of native height (was W*img_asp).
# Image is centred inside its (now shorter) container so aspect is preserved
# but B-F row occupies more of the total figure.
W          <- 18                   # composite width, inches
H_BOTTOM   <- 9                    # B-F row height (v1 shape preserved)
H_A_NATIVE <- W * img_asp          # A spans full width at native aspect
# 2026-07-08 (revert): the collaborators asked for A back at full size — image spans the full
# 18-in width at native aspect (12 in tall).
H_A        <- H_A_NATIVE
H          <- H_A + H_BOTTOM

p_A <- build_p_A(W_container = W, H_container = H_A)

fig1 <- p_A
fig2 <- bottom_block

cat(sprintf("Figure 1: %.1f x %.1f in\nFigure 2: %.1f x %.1f in, panels A-E\n",
            W, H_A, W, H_BOTTOM))

save_figure <- function(plot, n, width, height) {
  stem <- paste0("Figure_", n)
  dir  <- asd_figure_dir(n)
  pdf  <- file.path(dir, paste0(stem, ".pdf"))
  png  <- file.path(dir, paste0(stem, ".png"))
  rds  <- file.path(dir, paste0(stem, ".rds"))
  ggsave(pdf, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  saveRDS(plot, rds)
  ggsave(png, plot, width = width, height = height, dpi = 150, bg = "white")
  cat("wrote:\n  ", pdf, "\n  ", png, "\n  ", rds, "\n", sep = "")
}

save_figure(fig1, 1, W, H_A)
save_figure(fig2, 2, W, H_BOTTOM)
