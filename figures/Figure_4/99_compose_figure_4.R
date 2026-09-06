root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Ver-2 Figure 3 — full composite
#   ┌─────────────┬─────────────┐
#   │   A heat    │  B treatab. │
#   ├─────────────┴─────────────┤
#   │      C  CHD8 lollipop     │
#   └───────────────────────────┘
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

VER1    <- asd_source_v1()
VER2    <- asd_source()
V1_RDS  <- file.path(VER1, "Intermediate_tables/Figure_3")
V2_RDS  <- asd_intermediate_dir()
OUT_DIR <- asd_figure_dir(4)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Panel A: brain heatmap (ver-1 reused) ---------------------------------
p_A <- readRDS(asd_panel("Figure_3_brain_heatmap_ggplot.rds", V1_RDS)) +
  labs(tag = "A") +
  theme(legend.position       = "top",
        legend.direction      = "horizontal",
        legend.justification  = "left",
        legend.title          = element_text(size = 10),
        legend.text           = element_text(size = 9),
        legend.key.width      = unit(1.4, "cm"),
        legend.key.height     = unit(0.4, "cm"),
        legend.box.margin     = margin(0, 0, 2, 0),
        plot.tag              = element_text(face = "bold", size = 16,
                                             family = "Nimbus Sans"),
        plot.tag.position     = c(0, 1),
        plot.margin           = margin(5, 0, 5, 8))

# --- Panel B: ver-2 treatability bar (shared y with A — strip labels) ------
p_B <- readRDS(asd_panel("Figure_3B_treatability_bar_ver2_ggplot.rds", V2_RDS)) +
  labs(tag = NULL) +   # shares panel A with the heatmap to its left
  theme(axis.text.y           = element_blank(),
        axis.ticks.y          = element_blank(),
        axis.line.y           = element_blank(),
        legend.position       = "top",
        legend.direction      = "horizontal",
        legend.justification  = "left",
        legend.title          = element_text(size = 10),
        legend.text           = element_text(size = 9),
        legend.box.margin     = margin(0, 0, 2, 0),
        panel.grid.major.y    = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor      = element_blank(),
        plot.tag              = element_text(face = "bold", size = 16,
                                             family = "Nimbus Sans"),
        plot.tag.position     = c(0, 1),
        plot.margin           = margin(5, 8, 5, 0))

top_row <- (p_A | p_B) + plot_layout(widths = c(1.1, 1.4))

# --- Panel C: CHD8 lollipop (ver-2) ---------------------------------------
p_C <- readRDS(asd_panel("Figure_3C_CHD8_lollipop_ver2_ggplot.rds", V2_RDS)) +
  labs(tag = "B") +
  theme(plot.tag          = element_text(face = "bold", size = 16,
                                         family = "Nimbus Sans"),
        plot.tag.position = c(0, 1),
        plot.margin       = margin(10, 8, 5, 8))

# --- Compose --------------------------------------------------------------
composite <- (top_row / p_C) + plot_layout(heights = c(2, 1))

W <- 16; H <- 18
out_pdf <- file.path(OUT_DIR, "Figure_4.pdf")
out_png <- file.path(OUT_DIR, "Figure_4.png")
out_rds <- file.path(OUT_DIR, "Figure_4.rds")
ggsave(out_pdf, composite, width = W, height = H,
       device = cairo_pdf, bg = "transparent")
ggsave(out_png, composite, width = W, height = H,
       dpi = 150, bg = "white")
saveRDS(composite, out_rds)
cat("wrote:\n  ", out_pdf, "\n  ", out_png, "\n  ", out_rds, "\n", sep = "")
