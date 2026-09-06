root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# All 20 SFARI-1 gene lollipops as INDIVIDUAL panels (not the supp grid).
# Same visual style as the Fig 3C CHD8 / ANK2 panels — full-size fonts,
# 14 × 6 in canvas, one PDF + PNG + RDS per gene.
#
# Reuses the UniProt cache written by 05_supp_all_lollipops_ver2.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(scales)
  library(stringr)
})

# ggrepel resolves label collisions with random jitter, so an unseeded run places
# domain labels differently every time. Nothing in the data moves -- only where a
# label lands and how its leader line is drawn -- but it is the reason these
# figures were not byte-reproducible. Seeded here so they are.
set.seed(42)

VER1     <- asd_source_v1()
VER2     <- asd_source()
INT_F3   <- asd_intermediate_dir()
INT_MO   <- asd_intermediate_dir()
OUT_DIR <- asd_out_dir("Figures", "Supplementary", "lollipops_individual")
UP_CACHE <- file.path(INT_F3, "uniprot_cache_ver2.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(asd_code("ver1", "_helpers.R"))
source(asd_code("_palette.R"))

# --- Data + classification -------------------------------------------------
df_cds       <- readRDS(asd_db_snapshot("df_cds.rds"))
ids_nonsense <- readRDS(file.path(INT_F3, "nonsense_rescue_894_ids.rds"))
ids_siftopt  <- as.integer(readRDS(file.path(INT_MO,
                                             "MissenseOptimization_raw.rds"))$id)
top_genes    <- readRDS(file.path(VER1, "Intermediate_tables/Figure_3/top_genes.rds"))$gene_symbol

stopifnot(file.exists(UP_CACHE))
up_cache <- readRDS(UP_CACHE)

# --- Aesthetic scales (identical to Panel C) -------------------------------
treat_colors <- c(
  "Direct Repair"     = unname(fig_palette["teal_green"]),
  "Nonsense Rescue"   = unname(fig_palette["light_green"]),
  "SIFT Optimization" = unname(fig_palette["yellow"]),
  "Untreatable"       = "grey85"
)
treat_shapes <- c(
  "Direct Repair"     = 22,   # square
  "Nonsense Rescue"   = 24,   # triangle
  "SIFT Optimization" = 23,   # diamond
  "Untreatable"       = 21    # circle
)

# --- Plot builder (CHD8/ANK2 sizing) ---------------------------------------
build_lollipop <- function(gene_name, up_entry) {
  if (is.null(up_entry)) return(NULL)
  p_length <- up_entry$p_length
  if (is.na(p_length) || p_length <= 0) return(NULL)

  df <- df_cds %>%
    filter(gene_symbol == gene_name,
           !is.na(protein_position), !is.na(CADD_phred)) %>%
    mutate(
      protein_position = as.numeric(protein_position),
      cons_first = sapply(str_split(New_Consequence, ","), `[`, 1),
      consequence_simple = case_when(
        cons_first == "Missense"   ~ "Missense",
        cons_first == "Splice"     ~ "Splice",
        cons_first == "Nonsense"   ~ "StopGained",
        TRUE                       ~ NA_character_
      ),
      is_G_to_A = (REF == "G" & ALT == "A"),
      treat_class = case_when(
        is_G_to_A               ~ "Direct Repair",
        id %in% ids_nonsense    ~ "Nonsense Rescue",
        id %in% ids_siftopt     ~ "SIFT Optimization",
        TRUE                    ~ "Untreatable"
      )
    ) %>%
    filter(!is.na(consequence_simple))

  if (nrow(df) == 0) return(NULL)

  df$treat_class <- factor(df$treat_class,
                           levels = c("Direct Repair", "Nonsense Rescue",
                                      "SIFT Optimization", "Untreatable"))

  domains <- up_entry$domain_df %>%
    mutate(clean_name = trimws(gsub(" [0-9]+$", "", domain_name)),
           display_name = case_when(
             grepl("disorder", clean_name, ignore.case = TRUE) ~ "Disordered Region",
             TRUE ~ clean_name),
           xmin = start, xmax = end,
           strip_col = "grey25",
           ymin = -5, ymax = 0) %>%
    arrange(xmin)
  domains$strip_col[domains$display_name == "Disordered Region"] <- "grey85"

  label_domains <- domains %>%
    filter(display_name != "Disordered Region") %>%
    group_by(display_name) %>%
    summarise(x_mid = median((xmin + xmax) / 2), .groups = "drop") %>%
    mutate(domain_label = str_wrap(display_name, width = 15),
           y_anchor     = -5)

  max_y <- max(df$CADD_phred, 40, na.rm = TRUE) + 5

  # Immediately before the repel call, not at the top of the file: any code
  # between the two advances the RNG, so a distant seed does not pin the
  # label layout. This is why an earlier seeding attempt fixed Figure S1 but
  # left Figures 1 and 3 varying run to run.
  set.seed(42)
  ggplot() +
    annotate("rect", xmin = 0, xmax = p_length,
             ymin = -5, ymax = 0,
             fill = "gray97", colour = "gray70", linewidth = 0.3) +
    geom_rect(data = domains,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = domains$strip_col,
              colour = "black", alpha = 0.6, linewidth = 0.3) +
    geom_text_repel(seed = 42, data = label_domains,
                    aes(x = x_mid, y = y_anchor, label = domain_label),
                    size = 3.2, fontface = "bold", colour = "black",
                    nudge_y = -18, direction = "both", vjust = 1,
                    segment.size = 0.5, segment.colour = "gray40",
                    min.segment.length = 0, max.overlaps = 80,
                    force = 4, force_pull = 0.3,
                    box.padding = 0.35, point.padding = 0.2) +
    geom_segment(data = df,
                 aes(x = protein_position, xend = protein_position,
                     y = 0, yend = CADD_phred),
                 colour = "gray75", linewidth = 0.5, alpha = 0.75) +
    geom_point(data = df,
               aes(x = protein_position, y = CADD_phred,
                   fill = treat_class, shape = treat_class),
               colour = "black", size = 3.2, stroke = 0.4, alpha = 0.95) +
    geom_hline(yintercept = 25, linetype = "dashed",
               colour = "red2", alpha = 0.6, linewidth = 0.4) +
    scale_fill_manual(values = treat_colors, name = "Treatment class") +
    scale_shape_manual(values = treat_shapes, name = "Treatment class") +
    scale_x_continuous(limits = c(0, p_length),
                       expand = c(0.01, 0), labels = scales::comma) +
    coord_cartesian(ylim = c(-40, max_y), clip = "off") +
    labs(title = gene_name,
         subtitle = paste0("Amino acid position (total: ", p_length, " AA)"),
         x = NULL, y = "CADD PHRED score") +
    theme_publication(base_size = 12, y_bold = TRUE) +
    theme(plot.title       = element_text(face = "bold", size = 15,
                                          family = "Nimbus Sans"),
          plot.subtitle    = element_text(colour = "gray40", size = 10),
          legend.position  = "top",
          legend.box       = "horizontal",
          legend.title     = element_text(size = 10),
          legend.text      = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "gray95"),
          plot.margin      = margin(8, 8, 8, 8))
}

# --- Loop + save each gene -------------------------------------------------
summary_rows <- list()
for (g in top_genes) {
  cat("[", g, "] ", sep = "")
  p <- build_lollipop(g, up_cache[[g]])
  if (is.null(p)) { cat("skipped (no data)\n"); next }

  panel_id <- paste0(g, "_lollipop_ver2")
  pdf_path <- file.path(OUT_DIR, paste0(panel_id, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(panel_id, ".png"))
  rds_path <- file.path(OUT_DIR, paste0(panel_id, ".rds"))
  ggsave(pdf_path, p, width = 14, height = 6, device = cairo_pdf,
         bg = "transparent")
  ggsave(png_path, p, width = 14, height = 6, dpi = 150, bg = "white")
  saveRDS(p, rds_path)
  cat("ok\n")
  summary_rows[[g]] <- data.frame(
    gene = g,
    n_variants = sum(df_cds$gene_symbol == g &
                       !is.na(df_cds$protein_position) &
                       !is.na(df_cds$CADD_phred)),
    p_length = up_cache[[g]]$p_length,
    stringsAsFactors = FALSE
  )
}

summary_tbl <- bind_rows(summary_rows)
write.csv(summary_tbl, file.path(OUT_DIR, "_summary.csv"), row.names = FALSE)
cat("\nWrote", nrow(summary_tbl), "individual lollipops to:\n  ", OUT_DIR, "\n", sep = "")
cat("Summary CSV:", file.path(OUT_DIR, "_summary.csv"), "\n")
