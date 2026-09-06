# Project-wide helpers for the ASD ADAR paper.
# Source from any rmd:  source("../_helpers.R")  or  source(asd_code("ver1", "_helpers.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# --- Color anchors (publication guide file 01) ---
fig_colors <- list(
  gtex        = "#1F77B4",  # bulk-reference (publication blue)
  single_cell = "#D4AF37",  # primary signal (gold) — also used for "G>A treatable"
  alu         = "#000000"   # Alu / inner overlay / black
)

# Mismatch-panel palette: A>G hallmark in gold, all others in soft green.
# For ASD we adapt this so G>A (the ADAR-treatable A>I-on-the-other-strand class)
# gets the gold and everything else gets the soft fill.
mismatch_fill <- c(
  "A>C" = "#AAD297", "A>G" = "#AAD297", "A>T" = "#AAD297",
  "C>A" = "#AAD297", "C>G" = "#AAD297", "C>T" = "#AAD297",
  "G>A" = "#D4AF37",  # ADAR-treatable mismatch — gold
  "G>C" = "#AAD297", "G>T" = "#AAD297",
  "T>A" = "#AAD297", "T>C" = "#AAD297", "T>G" = "#AAD297"
)

# Bin-gradient palette for off-target / bystander panels (Panel D / E).
# Six steps: green → yellow → red, conveying "safe → dangerous" as the bin
# count grows.  Used when bars are coloured by the *bin* (no consequence
# stacking), with the legend hidden (the x-axis labels are the legend).
bin_gradient <- c(
  "0"  = "#2C9E25",
  "1"  = "#9EDC3B",
  "2"  = "#FCF708",
  "3"  = "#F6A609",
  "4"  = "#DD0005",
  "5+" = "#95303E"
)

# Alternative "severity decreasing" palette: dark red → pink across bins.
# Bin 0 (the priority bucket — most editable variants) gets the strongest
# red; the colour fades to pink as bin index increases, de-emphasising
# variants with many off-targets.
bin_red_pink <- c(
  "0"  = "#7A0019",   # very dark red
  "1"  = "#A8132E",   # dark red
  "2"  = "#D62246",   # red
  "3"  = "#E55B7B",   # rose
  "4"  = "#F08CA1",   # pink
  "5+" = "#FBCBD0"    # light pink
)

# Stacked-by-consequence palette (panel B inner colors).
# Palette agreed with the collaborators for the ASD paper — kept as named R colors
# so it's easy to swap names directly during review with her.
consequence_fill <- c(
  "Missense"   = "orange",
  "Nonsense"   = "red2",
  "Splice"     = "cornflowerblue",
  "Start Lost" = "blue",
  "Synonymous" = "forestgreen",
  "Frameshift" = "grey60"     # not used in current panel B (no SNV frameshifts) — kept for completeness
)

# --- Publication theme (publication guide file 01) ---
theme_publication <- function(base_size = 12,
                              base_family = "Nimbus Sans",
                              x_angle = 0,
                              y_bold = TRUE) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text              = ggplot2::element_text(color = "black", family = base_family),
      axis.text         = ggplot2::element_text(color = "black"),
      axis.title        = ggplot2::element_text(color = "black"),
      axis.title.y      = ggplot2::element_text(
        color = "black", face = if (y_bold) "bold" else "plain"),
      axis.line         = ggplot2::element_line(color = "black"),
      axis.ticks        = ggplot2::element_line(color = "black"),
      axis.text.x       = ggplot2::element_text(
        color = "black", angle = x_angle,
        hjust = if (x_angle == 0) 0.5 else 1,
        vjust = if (x_angle == 0) 0.5 else 1),
      legend.text       = ggplot2::element_text(color = "black"),
      legend.title      = ggplot2::element_text(color = "black"),
      plot.title        = ggplot2::element_text(color = "black"),
      plot.subtitle     = ggplot2::element_text(color = "black"),
      strip.text        = ggplot2::element_text(color = "black"),
      strip.background  = ggplot2::element_blank(),
      plot.background   = ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank()
    )
}

# --- Dual save: cairo_pdf + ggplot RDS (publication guide file 01) ---
# RDS lets a composer rebuild composites without re-running upstream R.
save_panel <- function(p, panel_id, pdf_dir, rds_dir,
                       width = 8, height = 6) {
  dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

  pdf_file <- file.path(pdf_dir, paste0(panel_id, ".pdf"))
  rds_file <- file.path(rds_dir, paste0(panel_id, "_ggplot.rds"))

  ggplot2::ggsave(pdf_file, p,
                  width = width, height = height,
                  device = cairo_pdf,
                  bg = "transparent")
  saveRDS(p, rds_file)
  message("Saved panel: ", pdf_file)
  message("           : ", rds_file)
  invisible(list(pdf = pdf_file, rds = rds_file))
}

# --- SIFT boundary --------------------------------------------------------
# Whether a score of exactly 0.05 counts as tolerated is a convention, not a
# measurement. Ensembl reports SIFT to two decimal places, so a stored 0.05 is
# any true value in [0.045, 0.055), and this run's own raw VEP output contains
# both SIFT=deleterious(0.05) and SIFT=tolerated(0.05) -- VEP does not apply one
# reading consistently either. 750 variants sit exactly on the boundary.
#
# Defined once here so the comparison and the printed threshold cannot drift
# apart. The default reproduces the published figures exactly; set
# SIFT_TOLERATED_STRICT=TRUE to treat 0.05 as deleterious instead, which is the
# reading the Methods sentence describes.
SIFT_CUT <- 0.05

sift_strict_default <- function() {
  toupper(Sys.getenv("SIFT_TOLERATED_STRICT", "FALSE")) %in% c("TRUE", "1", "YES")
}

#' TRUE where `x` counts as tolerated under the active convention.
sift_tolerated <- function(x, strict = sift_strict_default()) {
  if (isTRUE(strict)) x > SIFT_CUT else x >= SIFT_CUT
}

#' Axis/bin labels that must agree with `sift_tolerated()`, since panels print
#' the threshold rather than merely using it.
sift_labels <- function(strict = sift_strict_default()) {
  if (isTRUE(strict)) list(tolerated = "SIFT>0.05",  deleterious = "SIFT\u22640.05")
  else                list(tolerated = "SIFT\u22650.05", deleterious = "SIFT<0.05")
}
