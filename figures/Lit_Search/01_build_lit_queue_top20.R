root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Build the ver-2 literature-search queue for TOP-20 SFARI-1 treatable variants.
#
# Universe (per the collaborators, 2026-07-08):
#   1. Top-20 SFARI-1 genes only (same as Fig 3 gene list)
#   2. TREATABLE variants only: Direct Repair + Nonsense Rescue + SIFT-Opt
#      (Untreatable excluded)
#
# Also cross-check overlap with ver-1's completed lit-search outputs so we
# know how much genuinely needs fresh audit vs. can be lifted verbatim.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(stringr) })

VER1        <- asd_source_v1()
VER2        <- asd_source()
INT_F3      <- asd_intermediate_dir()
INT_MO      <- asd_intermediate_dir()
INT_LIT_V2  <- asd_intermediate_dir()
LIT_V1_DIR  <- file.path(VER1, "Results_Plots_Tables", "Lit_Search", "Correctable")
dir.create(INT_LIT_V2, showWarnings = FALSE, recursive = TRUE)

# --- Data + treatment-class classification ---------------------------------
df_cds        <- readRDS(asd_db_snapshot("df_cds.rds"))
ids_nonsense  <- readRDS(file.path(INT_F3, "nonsense_rescue_894_ids.rds"))
ids_siftopt   <- as.integer(readRDS(file.path(INT_MO, "MissenseOptimization_raw.rds"))$id)
top_genes_tbl <- readRDS(file.path(VER1, "Intermediate_tables/Figure_3/top_genes.rds"))
top_genes     <- top_genes_tbl$gene_symbol

# Same classification precedence as Panel B / Panel C
df_treat <- df_cds %>%
  mutate(
    is_G_to_A = (REF == "G" & ALT == "A"),
    treat_class = case_when(
      is_G_to_A               ~ "Direct Repair",
      id %in% ids_nonsense    ~ "Nonsense Rescue",
      id %in% ids_siftopt     ~ "SIFT Optimization",
      TRUE                    ~ "Untreatable"
    )
  )

# --- Queue: top-20 genes + treatable only ----------------------------------
queue <- df_treat %>%
  filter(gene_symbol %in% top_genes,
         treat_class != "Untreatable") %>%
  mutate(
    coord_key = paste(chr, start, REF, ALT, sep = ":"),   # ver-1 style key
    hgvs_p    = paste0("p.",
                       ifelse(is.na(ref_aa_name), "?", ref_aa_name),
                       ifelse(is.na(protein_position), "?", protein_position),
                       ifelse(is.na(alt_aa_name), "?", alt_aa_name))
  ) %>%
  select(id, gene_symbol, chr, start, REF, ALT,
         protein_position, ref_aa_name, alt_aa_name, hgvs_p,
         Consequence, New_Consequence,
         CADD_phred, SIFT_score, SIFT_pred,
         sfari_gene_score, treat_class, coord_key)

cat("--- Queue summary ---\n")
cat("Top-20 gene total universe (from df_cds):",
    sum(df_cds$gene_symbol %in% top_genes), "\n")
cat("Top-20 treatable queue rows:", nrow(queue), "\n\n")

cat("Per-gene breakdown (treatable only):\n")
per_gene <- queue %>%
  count(gene_symbol, treat_class) %>%
  tidyr::pivot_wider(names_from = treat_class, values_from = n, values_fill = 0L) %>%
  mutate(Total = rowSums(across(-gene_symbol))) %>%
  arrange(desc(Total))
print(per_gene, n = 20)

cat("\nBy treatment class (across top-20):\n")
print(queue %>% count(treat_class))

# --- Overlap with ver-1 completed .txt files -------------------------------
# Ver-1 lit-search outputs: per-variant .txt with coord_key filename pattern.
# Check what's already done and what's fresh.
ver1_txt_files <- list.files(LIT_V1_DIR, pattern = "\\.txt$", full.names = FALSE)
cat("\n--- Overlap check with ver-1 lit-search outputs ---\n")
cat("Ver-1 .txt files at", LIT_V1_DIR, ":", length(ver1_txt_files), "\n")

# Filename usually contains the coord_key or a similar id — sample a few
if (length(ver1_txt_files) > 0) {
  cat("first 5 ver-1 filenames:\n"); print(head(ver1_txt_files, 5))
}

# Try to build a coord_key -> ver-1 filename map from the ver-1 triage table
ver1_triage <- file.path(LIT_V1_DIR, "top20_triage_summary_table.csv")
if (file.exists(ver1_triage)) {
  v1 <- read.csv(ver1_triage, stringsAsFactors = FALSE)
  cat("\nver-1 triage table:", nrow(v1), "rows, columns:", paste(names(v1), collapse = ", "), "\n")
} else {
  cat("\nver-1 triage table not found at", ver1_triage, "\n")
}

# --- Save queue -----------------------------------------------------------
saveRDS(queue, file.path(INT_LIT_V2, "queue_top20_treatable.rds"))
write.csv(queue, file.path(INT_LIT_V2, "queue_top20_treatable.csv"),
          row.names = FALSE)
cat("\nwrote:\n  ", file.path(INT_LIT_V2, "queue_top20_treatable.rds"),
    "\n  ", file.path(INT_LIT_V2, "queue_top20_treatable.csv"), "\n", sep = "")
