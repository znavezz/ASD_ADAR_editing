root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Figure 1B (ver-2) — ClinVar × SFARI donut grid with the project Spectral palette.
#
# Reads the cached plotdata RDS from ver-1 (no DB hit), swaps the two slice
# colours, and re-renders.
#
# 2026-07-08 iteration:
#   - keep only the "white" variant (other = white; no yellow)
#   - inner donut-hole outline linewidth = 0.5 (was 0.3)
#   - centre label = fraction "G>A / Total" (e.g. "138/606")
#   - small donuts (radius < 0.20) → label pushed outside via geom_label_repel
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(scatterpie)
  library(ggforce)
  library(ggrepel)
})

# ggrepel resolves label collisions with random jitter, so an unseeded run places
# domain labels differently every time. Nothing in the data moves -- only where a
# label lands and how its leader line is drawn -- but it is the reason these
# figures were not byte-reproducible. Seeded here so they are.
set.seed(42)

VER1     <- asd_source_v1()
VER2     <- asd_source()
PALETTE  <- asd_code("_palette.R")
HELPERS  <- asd_code("ver1", "_helpers.R")   # theme_publication, save_panel
IN_RDS   <- asd_panel("Figure_1C_clinvar_donut_plotdata.rds",
                      asd_source_v1("Intermediate_tables/Figure_1"))
OUT_INT <- asd_intermediate_dir()
OUT_MAIN <- asd_panels_dir()

source(HELPERS)
source(PALETTE)

bundle          <- readRDS(IN_RDS)
pie_plot_data   <- bundle$pie_plot_data
category_levels <- bundle$category_levels

# --- fraction label, all sitting a little above the donut ------------------
pie_plot_data <- pie_plot_data %>%
  mutate(
    label_frac = sprintf("%s/%s",
                         formatC(`G>A`,  format = "d", big.mark = ","),
                         formatC(Total, format = "d", big.mark = ",")),
    label_y    = y_pos + radius + 0.10
  )

# Slice column order controls which side of the pie the G>A wedge lands.
# "other" first → other draws clockwise from 12 o'clock into the right half,
# leaving G>A on the LEFT half (the agreed layout).
slice_cols   <- c("other", "G>A")
b_fill_white <- c("G>A" = unname(b_fill["G>A"]), "other" = "white")

p <- ggplot() +
  geom_scatterpie(
    data = pie_plot_data,
    aes(x = x_pos, y = y_pos, r = radius),
    cols = slice_cols,
    color = "black", linewidth = 0.7, alpha = 0.95
  ) +
  # Fraction labels — all above their own pie, uniform style
  geom_text(
    data = pie_plot_data,
    aes(x = x_pos, y = label_y, label = label_frac),
    fontface = "bold", size = 3.6, family = "Nimbus Sans",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = b_fill_white,
                    breaks = c("G>A", "other"),
                    labels = c("G>A", "other"),
                    name   = NULL) +
  scale_x_continuous(breaks = 1:3, labels = category_levels,
                     expand = expansion(mult = 0.25)) +
  scale_y_continuous(breaks = 1:3, labels = c("1", "2", "3"),
                     expand = expansion(mult = 0.30)) +
  coord_fixed() +
  labs(x = "ClinVar classification",
       y = "SFARI gene score",
       tag = "B") +
  theme_publication(base_size = 12, x_angle = 0, y_bold = TRUE) +
  theme(
    legend.position    = "top",
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major   = element_line(color = "gray95"),
    panel.grid.minor   = element_blank(),
    plot.tag           = element_text(face = "bold", size = 16,
                                      family = "Nimbus Sans"),
    plot.tag.position  = c(0, 1)
  )

save_panel(p, "Figure_1B_clinvar_donut_ver2_white",
           pdf_dir = OUT_MAIN, rds_dir = OUT_INT,
           width = 7, height = 6)
ggsave(file.path(OUT_MAIN, "Figure_1B_clinvar_donut_ver2_white.png"),
       p, width = 7, height = 6, dpi = 200, bg = "white")

# Keep the composite-consumed filename in sync
saveRDS(p, file.path(OUT_INT, "Figure_1B_clinvar_donut_ver2_ggplot.rds"))

cat("Panel B ver-2 (white, fractions above donut, 0.5 ring) written.\n")
print(pie_plot_data[, c("sfari_gene_score","Category","label_frac","radius","label_y")])
