root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Ver-2 Figure 2 composite — vertical: A illustration on top, B scatter below.
#
# A = Figure2A_illustration.png (Nonsense rescue vs Missense optimization diagram)
# B = Ver-2 missense-optimization scatter (variant SIFT × post-edit SIFT), v7
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
  library(png);     library(grid)
})

VER2    <- asd_source()
OUT_DIR <- asd_figure_dir(3)
ILL_PNG <- asd_code("Figure2A_illustration.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Panel A: raster illustration, aspect-preserved -----------------------
img      <- readPNG(ILL_PNG)
img_h_px <- dim(img)[1]
img_w_px <- dim(img)[2]
img_asp  <- img_h_px / img_w_px
cat(sprintf("illustration: %d × %d px, aspect %.3f\n",
            img_w_px, img_h_px, img_asp))
img_grob <- rasterGrob(img, interpolate = TRUE)

build_p_A <- function(W_container, H_container) {
  container_asp <- H_container / W_container
  if (img_asp > container_asp) {
    img_w_frac <- container_asp / img_asp; img_h_frac <- 1
  } else {
    img_w_frac <- 1;                       img_h_frac <- img_asp / container_asp
  }
  xmin <- (1 - img_w_frac) / 2; xmax <- 1 - xmin
  ymin <- (1 - img_h_frac) / 2; ymax <- 1 - ymin

  ggplot() +
    annotation_custom(img_grob,
                      xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    labs(tag = "A") +
    theme_void() +
    theme(plot.tag          = element_text(face = "bold", size = 18,
                                           family = "Nimbus Sans"),
          plot.tag.position = c(0, 1),
          plot.margin       = margin(4, 4, 4, 4))
}

# --- Panel B: missense-optimization scatter (v7 has the "B" tag baked in) ---
# Override the tag so it matches Panel A's size/position exactly.
# Read the panel from where 03_panel_B_missense_scatter.R writes it. It used to
# read the recovered July snapshot under asd_source(), which is an *inputs* tree
# -- so regenerating the panel changed nothing here and the composite silently
# kept the stale render. Fall back to the snapshot only if the panel has not
# been built yet, so a fresh checkout still composes.
p_B_panel <- file.path(asd_panels_dir(), "Figure_B_missense_scatter_v7.rds")
p_B_snap  <- file.path(VER2, "Results_Plots_Tables/Main",
                       "Figure_B_missense_scatter_v7.rds")
p_B_src <- if (file.exists(p_B_panel)) p_B_panel else p_B_snap
cat("panel B from: ", p_B_src, "\n", sep = "")
p_B <- readRDS(p_B_src)
p_B <- p_B +
  labs(tag = "B") +
  theme(plot.tag          = element_text(face = "bold", size = 18,
                                         family = "Nimbus Sans"),
        plot.tag.position = c(0, 1),
        plot.margin       = margin(4, 4, 4, 4))

# --- Layout: composite width chosen so A sits at native aspect --------------
W    <- 10                    # composite width, inches
H_A  <- W * img_asp           # A row height at native aspect (~5.1 in)
H_B  <- 4.5                   # smaller than A row on purpose (per the collaborators, 2026-07-08)
H    <- H_A + H_B

p_A <- build_p_A(W_container = W, H_container = H_A)

# Centre B inside its row with white padding on both sides so it visually
# reads smaller than A. Widths c(1.5, 7, 1.5) → B is 70% of composite width.
p_B_row <- (plot_spacer() | p_B | plot_spacer()) +
  plot_layout(widths = c(1.5, 7, 1.5))

composite <- (p_A / p_B_row) + plot_layout(heights = c(H_A, H_B))

cat(sprintf("composite: %.1f × %.1f in, A row = %.1f in\n", W, H, H_A))

out_pdf <- file.path(OUT_DIR, "Figure_3.pdf")
out_png <- file.path(OUT_DIR, "Figure_3.png")
out_rds <- file.path(OUT_DIR, "Figure_3.rds")
ggsave(out_pdf, composite, width = W, height = H,
       device = cairo_pdf, bg = "transparent")
ggsave(out_png, composite, width = W, height = H,
       dpi = 150, bg = "white")
saveRDS(composite, out_rds)
cat("wrote:\n  ", out_pdf, "\n  ", out_png, "\n  ", out_rds, "\n", sep = "")
