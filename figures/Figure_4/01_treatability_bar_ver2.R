root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Ver-2 Fig 3 Panel B — Treatability bar (top-20 SFARI-1 genes)
#
# CHANGES from ver-1 (figures/ver1/Figure_3/01_panel_B_treatability.rmd):
#   1. FOUR categories instead of three:
#        - Direct Repair            (unchanged, 2,421 total)
#        - Nonsense Rescue              (was 1,092; now 894 with nmd_escaping=false)
#        - SIFT Optimization  [NEW]     (158 non-G>A missense w/ neighbor-edit rescue)
#        - Untreatable                  (remainder)
#      Classification precedence: G>A > Nonsense > SIFT-Opt > Untreatable.
#   2. Palette: the project Spectral palette.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(tibble)
  library(httr); library(jsonlite)
})

VER1     <- asd_source_v1()
VER2     <- asd_source()
source(asd_code("ver1", "_helpers.R"))
source(asd_code("_palette.R"))

IN_INT        <- file.path(VER2, "Intermediate_tables", "Figure_3")  # recovered: read only
INT_DIR       <- asd_intermediate_dir()                              # ours: written to
INT_MISSENSE  <- asd_intermediate_dir()
OUT_MAIN <- asd_panels_dir()
dir.create(INT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)

# --- Palette for 4 categories ---------------------------------------------
treat_colors <- c(
  "Direct Repair" = unname(fig_palette["teal_green"]),  # #66C2A5
  "Nonsense Rescue"   = unname(fig_palette["light_green"]), # #ABDDA4
  "SIFT Optimization" = unname(fig_palette["yellow"]),      # #FEE08B
  "Untreatable"       = "grey90"
)

# --- Pull 894 Nonsense-Rescue IDs (NMD non-escape) -------------------------
# Local path, not asd_panel(): this is a write target guarded by file.exists(),
# so resolving it to the recovered copy would make the rebuild unreachable.
nonsense_ids_rds <- file.path(asd_intermediate_dir(), "nonsense_rescue_894_ids.rds")

if (file.exists(nonsense_ids_rds)) {
  ids_nonsense <- readRDS(nonsense_ids_rds)
  cat("reused cached Nonsense-Rescue IDs:", length(ids_nonsense), "\n")
} else {
  env_lines <- readLines(asd_env_file())
  secret    <- sub("^HASURA_ADMIN_SECRET=\"?([^\"]+)\"?", "\\1",
                   grep("^HASURA_ADMIN_SECRET=", env_lines, value = TRUE))
  HASURA    <- "http://localhost:8789/v1/graphql"

  q <- '{
    variants(
      where: {
        class: {_eq: "SNV"},
        nmd_escaping_variant: {_eq: false},
        _not: { _or: [
          {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}},
          {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
        ]},
        variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}},
        gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}},
        variants_consequences: {consequence: {name: {_eq: "StopGained"}}}
      }
      limit: 100000
    ) { id }
  }'
  resp   <- POST(HASURA, add_headers(`x-hasura-admin-secret` = secret),
                 content_type_json(),
                 body = toJSON(list(query = q), auto_unbox = TRUE),
                 timeout(600))
  parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
  if (!is.null(parsed$errors)) stop("gql: ", toJSON(parsed$errors, auto_unbox = TRUE))
  ids_nonsense <- as.integer(parsed$data$variants$id)
  saveRDS(ids_nonsense, nonsense_ids_rds)
  cat("pulled Nonsense-Rescue IDs (NMD non-escape):", length(ids_nonsense), "\n")
}

# --- Reuse cached 158 SIFT-Optimization IDs --------------------------------
mo_raw <- readRDS(file.path(INT_MISSENSE, "MissenseOptimization_raw.rds"))
ids_siftopt <- as.integer(mo_raw$id)
cat("SIFT-Optimization IDs (from cache):", length(ids_siftopt), "\n")

# --- Load df_cds + top-20 SFARI-1 gene list --------------------------------
df_cds        <- readRDS(asd_db_snapshot("df_cds.rds"))
top_genes_tbl <- readRDS(asd_panel("top_genes.rds",
                                  asd_source_v1("Intermediate_tables/Figure_3")))
top_genes     <- top_genes_tbl$gene_symbol

# --- Classify each df_cds row into 4 categories ----------------------------
# Precedence: G>A > Nonsense Rescue > SIFT Optimization > Untreatable
df_cds <- df_cds %>%
  mutate(
    is_G_to_A_gs = (REF == "G" & ALT == "A"),                     # gene-sense G>A
    status = case_when(
      is_G_to_A_gs                     ~ "Direct Repair",
      id %in% ids_nonsense             ~ "Nonsense Rescue",
      id %in% ids_siftopt              ~ "SIFT Optimization",
      TRUE                             ~ "Untreatable"
    )
  )

cat("\nCategory totals across all 9,962 CDS+splice variants:\n")
print(table(df_cds$status))

# --- Aggregate per top-20 gene ----------------------------------------------
gene_data <- df_cds %>%
  filter(gene_symbol %in% top_genes) %>%
  group_by(gene_symbol) %>%
  mutate(total_variants = n()) %>%
  ungroup() %>%
  count(gene_symbol, status, total_variants, .drop = FALSE) %>%
  mutate(status = factor(status,
                         levels = c("Untreatable",
                                    "SIFT Optimization",
                                    "Nonsense Rescue",
                                    "Direct Repair")))

cat("\nAggregate top-20 gene table (first rows):\n"); print(head(gene_data, 12))

saveRDS(gene_data, file.path(INT_DIR, "Figure_3_treatability_bar_plotdata_ver2.rds"))
write.csv(gene_data,
          file.path(INT_DIR, "Source_data_Figure_3_treatability_bar_ver2.csv"),
          row.names = FALSE)

# --- Render — bars only; one number above each bar = total treatable ------
BAR_WIDTH <- 0.35

gene_data$gene_f <- with(gene_data,
                         reorder(gene_symbol, total_variants))

# Per-gene: total treatable = G>A + Nonsense Rescue + SIFT Optimization.
# Text sits just above the bar at the x where the coloured region ends
# (i.e. the boundary between the last treatable segment and Untreatable).
label_data <- gene_data %>%
  mutate(status_str = as.character(status)) %>%
  group_by(gene_symbol, gene_f) %>%
  summarise(
    treatable = sum(n[status_str != "Untreatable"]),
    total     = sum(n),
    .groups   = "drop"
  )

p <- ggplot(gene_data,
            aes(x = n, y = gene_f, fill = status)) +
  geom_col(width = BAR_WIDTH, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = treat_colors,
                    breaks = c("Direct Repair",
                               "Nonsense Rescue",
                               "SIFT Optimization",
                               "Untreatable")) +
  geom_text(data = label_data,
            aes(x = treatable, y = gene_f, label = treatable),
            inherit.aes = FALSE,
            position = position_nudge(y = 0.36),
            colour = "black", fontface = "bold", size = 3.4,
            family = "Nimbus Sans", vjust = 0) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)), labels = comma) +
  labs(x = "Number of variants", y = NULL,
       fill = NULL, tag = "B") +
  theme_publication(base_size = 12, x_angle = 0, y_bold = TRUE) +
  theme(
    axis.text.y       = element_text(face = "bold", size = 11, family = "Nimbus Sans"),
    legend.position   = "top",
    legend.justification = "left",
    plot.tag          = element_text(face = "bold", size = 16, family = "Nimbus Sans"),
    plot.tag.position = c(0, 1)
  )

save_panel(p, "Figure_3B_treatability_bar_ver2",
           pdf_dir = OUT_MAIN, rds_dir = INT_DIR,
           width = 10, height = 12)
ggsave(file.path(OUT_MAIN, "Figure_3B_treatability_bar_ver2.png"),
       p, width = 10, height = 12, dpi = 150, bg = "white")

cat("\nPanel B ver-2 written.\n")
