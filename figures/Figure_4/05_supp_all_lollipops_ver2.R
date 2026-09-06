root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Supplementary figure — ver-2 domain lollipops for all 20 top-SFARI-1 genes.
#
# Uses the same ver-2 design as Panel C:
#   SHAPE + COLOUR paired by treatment class:
#     Direct Repair     = filled square,   #66C2A5 teal
#     Nonsense Rescue   = filled triangle, #ABDDA4 light green
#     SIFT Optimization = filled diamond,  #FEE08B yellow
#     Untreatable       = filled circle,   grey55
#   HEIGHT = CADD PHRED, dashed 25 threshold
#   DOMAINS = uniform grey25, alpha 0.6 (Disordered = grey85)
#
# Layout: 3 cols × 7 rows grid, black border round each panel, one shared
# legend across the top.  Same shape as ver-1 supp figure.
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(tidyr)
  library(ggplot2); library(ggrepel); library(patchwork)
  library(stringr); library(scales); library(cowplot)
})

# ggrepel resolves label collisions with random jitter, so an unseeded run places
# domain labels differently every time. Nothing in the data moves -- only where a
# label lands and how its leader line is drawn -- but it is the reason these
# figures were not byte-reproducible. Seeded here so they are.
set.seed(42)

VER1     <- asd_source_v1()
VER2     <- asd_source()
INT_DIR <- asd_intermediate_dir()
INT_MO   <- asd_intermediate_dir()
OUT_MAIN <- asd_panels_dir()
OUT_COMP <- asd_figure_dir(3)
UP_CACHE <- file.path(INT_DIR, "uniprot_cache_ver2.rds")
dir.create(INT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_COMP, showWarnings = FALSE, recursive = TRUE)

source(asd_code("ver1", "_helpers.R"))
source(asd_code("_palette.R"))

# --- Data + classification (once) ------------------------------------------
df_cds        <- readRDS(asd_db_snapshot("df_cds.rds"))
ids_nonsense  <- readRDS(asd_panel("nonsense_rescue_894_ids.rds",
                                   asd_source("Intermediate_tables/Figure_3")))
ids_siftopt   <- as.integer(readRDS(file.path(INT_MO,
                                              "MissenseOptimization_raw.rds"))$id)
top_genes_tbl <- readRDS(asd_panel("top_genes.rds",
                                  asd_source_v1("Intermediate_tables/Figure_3")))
top_genes     <- top_genes_tbl$gene_symbol
cat("Genes to plot:", length(top_genes), "\n")

# --- UniProt fetch cache ---------------------------------------------------
fetch_uniprot <- function(gene_symbol) {
  url <- modify_url("https://rest.uniprot.org/uniprotkb/search",
                    query = list(
                      query = paste0("gene_exact:", gene_symbol,
                                     " AND organism_id:9606 AND reviewed:true"),
                      format = "json"))
  resp <- GET(url, timeout(30))
  parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = FALSE)
  features <- parsed$results$features[[1]]
  p_length <- as.integer(parsed$results$sequence$length[1])
  domain_df <- features %>%
    filter(type %in% c("Domain", "Region", "Zinc finger", "Motif", "Repeat")) %>%
    mutate(domain_name = ifelse(!is.na(description) & description != "",
                                description, type),
           start = as.numeric(location$start$value),
           end   = as.numeric(location$end$value)) %>%
    select(domain_name, start, end)
  domain_df$domain_name <- gsub(" domain", "", domain_df$domain_name, ignore.case = TRUE)
  list(domain_df = domain_df, p_length = p_length)
}

up_cache <- if (file.exists(UP_CACHE)) readRDS(UP_CACHE) else list()
for (g in top_genes) {
  if (is.null(up_cache[[g]])) {
    cat("  fetching UniProt:", g, "\n")
    up_cache[[g]] <- tryCatch(fetch_uniprot(g),
                              error = function(e) {
                                cat("    error:", conditionMessage(e), "\n")
                                NULL
                              })
    Sys.sleep(0.5)   # be nice to UniProt
  }
}
saveRDS(up_cache, UP_CACHE)

# --- Aesthetic scales ------------------------------------------------------
treat_colors <- c(
  "Direct Repair"     = unname(fig_palette["teal_green"]),
  "Nonsense Rescue"   = unname(fig_palette["light_green"]),
  "SIFT Optimization" = unname(fig_palette["yellow"]),
  "Untreatable"       = "grey85"   # light grey to match Panel B
)
treat_shapes <- c(
  "Direct Repair"     = 22,   # square
  "Nonsense Rescue"   = 24,   # triangle
  "SIFT Optimization" = 23,   # diamond
  "Untreatable"       = 21    # circle
)

# --- Plot builder ----------------------------------------------------------
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
                    size = 2.6, fontface = "bold", colour = "black",
                    nudge_y = -18, direction = "both", vjust = 1,
                    segment.size = 0.4, segment.colour = "gray40",
                    min.segment.length = 0, max.overlaps = 80,
                    force = 4, force_pull = 0.3,
                    box.padding = 0.35, point.padding = 0.2) +
    geom_segment(data = df,
                 aes(x = protein_position, xend = protein_position,
                     y = 0, yend = CADD_phred),
                 colour = "gray75", linewidth = 0.4, alpha = 0.75) +
    geom_point(data = df,
               aes(x = protein_position, y = CADD_phred,
                   fill = treat_class, shape = treat_class),
               colour = "black", size = 2.6, stroke = 0.3, alpha = 0.95) +
    geom_hline(yintercept = 25, linetype = "dashed",
               colour = "red2", alpha = 0.6, linewidth = 0.4) +
    scale_fill_manual(values = treat_colors, name = "Treatment class") +
    scale_shape_manual(values = treat_shapes, name = "Treatment class") +
    scale_x_continuous(limits = c(0, p_length),
                       expand = c(0.01, 0), labels = scales::comma) +
    coord_cartesian(ylim = c(-40, max_y), clip = "off") +
    labs(title = gene_name, x = NULL, y = "CADD PHRED score") +
    theme_publication(base_size = 11, y_bold = TRUE) +
    theme(plot.title       = element_text(face = "bold", size = 13,
                                          family = "Nimbus Sans"),
          legend.position  = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "gray95"),
          plot.margin      = margin(8, 8, 8, 8),
          plot.background  = element_rect(colour = "black", fill = NA,
                                          linewidth = 0.7))
}

# --- Build one panel per gene, save RDS, keep list --------------------------
panels <- list()
for (g in top_genes) {
  cat("[", g, "] building panel ... ", sep = "")
  p <- build_lollipop(g, up_cache[[g]])
  if (is.null(p)) { cat("skipped (no data)\n"); next }
  saveRDS(p, file.path(INT_DIR, paste0("supp_lollipop_", g, "_ggplot.rds")))
  panels[[g]] <- p
  cat("ok\n")
}
cat("\nBuilt", length(panels), "panels.\n")

# --- Build shared legend: shape / colour chips only ------------------------
build_legend_panel <- function() {
  y   <- 0.5
  pts <- data.frame(
    x     = c(0.16, 0.36, 0.56, 0.78),   # centered, evenly spaced
    shape = c(22,    24,    23,    21),
    fill  = c(unname(fig_palette["teal_green"]),
              unname(fig_palette["light_green"]),
              unname(fig_palette["yellow"]),
              "grey85"),
    label = c("Direct Repair",
              "Nonsense Rescue",
              "SIFT Optimization",
              "Untreatable"),
    stringsAsFactors = FALSE
  )
  ggplot() +
    geom_point(data = pts,
               aes(x = x, y = y),
               shape = pts$shape, fill = pts$fill, colour = "black",
               size = 7, stroke = 0.5) +
    geom_text(data = pts,
              aes(x = x + 0.015, y = y, label = label),
              hjust = 0, vjust = 0.5, family = "Nimbus Sans",
              fontface = "bold", size = 4.4) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    theme_void() +
    theme(plot.margin = margin(6, 6, 6, 6))
}

shared_legend <- build_legend_panel()
header <- patchwork::wrap_elements(full = shared_legend)

# --- Compose 3 × 7 grid -----------------------------------------------------
grid <- patchwork::wrap_plots(panels, ncol = 3, nrow = 7)
composite <- header / grid + patchwork::plot_layout(heights = c(0.07, 1))

W <- 30; H <- 36
out_pdf <- file.path(OUT_COMP, "Figure_S1.pdf")
out_png <- file.path(OUT_COMP, "Figure_S1.png")
out_rds <- file.path(OUT_COMP, "Figure_S1.rds")
ggsave(out_pdf, composite, width = W, height = H,
       device = cairo_pdf, bg = "transparent")
ggsave(out_png, composite, width = W, height = H,
       dpi = 120, bg = "white")
saveRDS(composite, out_rds)
cat("wrote:\n  ", out_pdf, "\n  ", out_png, "\n  ", out_rds, "\n", sep = "")
