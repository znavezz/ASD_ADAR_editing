root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Figure 1D + 1E (ver-2) — off-target and bystander bin gradients with
# the project Spectral palette.
#
# Reads the cached plotdata + variants RDS from ver-1 (no DB hit), rebuilds
# both gradient panels using `bin_gradient_spectral` (0=teal → 5+=red), saves
# PDF + PNG + ggplot RDS into ver-2 tree.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

VER1     <- asd_source_v1()
VER2     <- asd_source()
PALETTE  <- asd_code("_palette.R")
HELPERS  <- asd_code("ver1", "_helpers.R")

VAR_RDS       <- asd_panel("Figure_1D_offtargets_variants.rds",
                           asd_source_v1("Intermediate_tables/Figure_1"))
PLOTDATA_RDS  <- asd_panel("Figure_1D_offtargets_plotdata.rds",
                           asd_source_v1("Intermediate_tables/Figure_1"))

OUT_INT <- asd_intermediate_dir()
OUT_MAIN <- asd_panels_dir()

source(HELPERS)
source(PALETTE)

# -------------------------------------------------------------------------
# Reload cached plotdata bundle (totals_dna, totals_bys will be recomputed
# for bystander to guarantee identical binning to the ver-1 pipeline)
# -------------------------------------------------------------------------
bundle        <- readRDS(PLOTDATA_RDS)
totals_dna    <- bundle$totals_dna
bin_levels    <- bundle$bin_levels

variants_long <- readRDS(VAR_RDS)

# Recompute bystander totals (ver-1's panel_bystander chunk used a
# clean_subset filtered on !is.na(dna_hits_85))
bin_count <- function(x, top = 5L) {
  ifelse(is.na(x), NA_character_,
         ifelse(x >= top, paste0(top, "+"), as.character(x)))
}

clean_subset <- variants_long %>% filter(!is.na(dna_hits_85))
clean_subset$bystander_bin <- factor(
  bin_count(clean_subset$dna_bystanders_cadd25),
  levels = bin_levels
)

totals_bys <- clean_subset %>%
  filter(!is.na(bystander_bin)) %>%
  count(bystander_bin, name = "Total") %>%
  arrange(bystander_bin)

# -------------------------------------------------------------------------
# Gradient panel builder — matches ver-1's make_gradient_panel exactly
# except for the fill palette (bin_gradient → bin_gradient_spectral).
# -------------------------------------------------------------------------
make_gradient_panel <- function(totals, x_var, x_lab, tag_label) {
  ggplot(totals, aes(x = .data[[x_var]], y = Total, fill = .data[[x_var]])) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.4) +
    geom_text(aes(label = format(Total, big.mark = ",")),
              vjust = -0.4, fontface = "bold", size = 3.5,
              family = "Nimbus Sans") +
    scale_fill_manual(values = bin_gradient_spectral, guide = "none",
                      drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
    labs(x = x_lab, y = "Number of variants", tag = tag_label) +
    theme_publication(base_size = 12, x_angle = 0, y_bold = TRUE) +
    theme(panel.grid.major.x = element_blank(),
          plot.tag           = element_text(face = "bold", size = 16,
                                            family = "Nimbus Sans"),
          plot.tag.position  = c(0, 1))
}

# --- Panel D (DNA off-targets) -------------------------------------------
p_dna <- make_gradient_panel(
  totals_dna, "dna_bin",
  "Number of DNA off-target hits (85% identity)",
  "D"
)
save_panel(p_dna, "Figure_1D_offtargets_dna_gradient_ver2",
           pdf_dir = OUT_MAIN, rds_dir = OUT_INT,
           width = 8, height = 5)
ggsave(file.path(OUT_MAIN, "Figure_1D_offtargets_dna_gradient_ver2.png"),
       p_dna, width = 8, height = 5, dpi = 200, bg = "white")

# --- Panel E (bystander burden CADD25) -----------------------------------
p_bys <- make_gradient_panel(
  totals_bys, "bystander_bin",
  "Number of bystander edits (CADD-phred ≥ 25)",
  "E"
)
save_panel(p_bys, "Figure_1E_bystander_cadd25_gradient_ver2",
           pdf_dir = OUT_MAIN, rds_dir = OUT_INT,
           width = 8, height = 5)
ggsave(file.path(OUT_MAIN, "Figure_1E_bystander_cadd25_gradient_ver2.png"),
       p_bys, width = 8, height = 5, dpi = 200, bg = "white")

cat("DNA off-target totals:\n"); print(totals_dna)
cat("\nBystander (CADD25) totals:\n"); print(totals_bys)
cat("\nPanels D and E ver-2 written.\n")
