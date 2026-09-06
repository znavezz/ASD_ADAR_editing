root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Compare ver-2 top-20 treatable queue (631) against ver-1's completed
# triage table (439 rows). Report which variants are shared vs new so we know
# what needs fresh lit-search.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })

VER1       <- asd_source_v1()
VER2       <- asd_source()
INT_LIT_V2 <- asd_intermediate_dir()

queue_v2   <- readRDS(file.path(INT_LIT_V2, "queue_top20_treatable.rds"))
v1_triage  <- read.csv(file.path(VER1,
                                 "Results_Plots_Tables/Lit_Search/Correctable",
                                 "top20_triage_summary_table.csv"),
                       stringsAsFactors = FALSE)

# Build a coord_key for the ver-1 triage table — it uses `chr_pos` in
# "chr:pos" form.  Combine with ref/alt from ref_codon/alt_codon? No —
# simpler: use chr_pos + ref_codon + alt_codon to disambiguate.
# But we have ref_aa/alt_aa/protein_change; the safest cross-key is chr:pos + REF + ALT.
# Ver-1 triage doesn't carry REF/ALT columns, so join on chr_pos alone
# (one variant per position within a gene is typical for missense/nonsense).
v1_triage <- v1_triage %>%
  mutate(coord_pos = chr_pos)   # already in "chr:pos" form

# Ver-2 queue coord_pos — prepend "chr" so it matches ver-1's format
queue_v2 <- queue_v2 %>%
  mutate(coord_pos = paste0("chr", chr, ":", start))

cat("--- Cross-tabulation ---\n")
cat("Ver-1 triage rows:", nrow(v1_triage), "\n")
cat("Ver-2 top-20 treatable rows:", nrow(queue_v2), "\n")

# Overlap on coord_pos
v1_pos     <- v1_triage$coord_pos
v2_pos     <- queue_v2$coord_pos

in_both    <- intersect(v1_pos, v2_pos)
v2_only    <- setdiff(v2_pos, v1_pos)
v1_only    <- setdiff(v1_pos, v2_pos)

cat("\nBy coord_pos ('chr:pos') only:\n")
cat("  In BOTH (already lit-searched in ver-1):", length(in_both), "\n")
cat("  ONLY in ver-2 queue (NEW, need fresh lit-search):", length(v2_only), "\n")
cat("  ONLY in ver-1 (dropped in ver-2 → NMD-escapers etc.):", length(v1_only), "\n")

# --- Per-gene NEW breakdown ------------------------------------------------
new_variants <- queue_v2 %>% filter(coord_pos %in% v2_only)
cat("\n--- NEW variants to lit-search, per gene × class ---\n")
per_gene_new <- new_variants %>%
  count(gene_symbol, treat_class) %>%
  tidyr::pivot_wider(names_from = treat_class, values_from = n, values_fill = 0L) %>%
  mutate(Total = rowSums(across(-gene_symbol))) %>%
  arrange(desc(Total))
print(per_gene_new, n = 20)

cat("\n--- NEW variants by class ---\n")
print(new_variants %>% count(treat_class))

# Save
saveRDS(new_variants, file.path(INT_LIT_V2, "queue_top20_NEW_only.rds"))
write.csv(new_variants, file.path(INT_LIT_V2, "queue_top20_NEW_only.csv"),
          row.names = FALSE)

# Dropped variants (ver-1 audited but no longer treatable in ver-2)
dropped <- v1_triage %>% filter(coord_pos %in% v1_only)
cat("\n--- DROPPED variants (ver-1 audited, ver-2 no longer treatable) ---\n")
if (nrow(dropped) > 0) {
  print(dropped %>% count(correctable_class, sort = TRUE))
} else {
  cat("(none)\n")
}
saveRDS(dropped, file.path(INT_LIT_V2, "dropped_from_ver1.rds"))

cat("\nwrote:\n  ",
    file.path(INT_LIT_V2, "queue_top20_NEW_only.rds"),
    "\n  ",
    file.path(INT_LIT_V2, "dropped_from_ver1.rds"), "\n", sep = "")
