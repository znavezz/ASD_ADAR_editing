root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Top-20 SFARI-1 gene treatability breakdown — counts + % per category.
# Same format as the CHD8 sentence in the previous manuscript.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

VER2    <- asd_source()
INT_DIR <- asd_intermediate_dir()
NOTES   <- file.path(VER2, "Notes")
dir.create(NOTES, showWarnings = FALSE, recursive = TRUE)

gene_data <- readRDS(file.path(INT_DIR, "Figure_3_treatability_bar_plotdata_ver2.rds"))

wide <- gene_data %>%
  select(gene_symbol, total = total_variants, status, n) %>%
  mutate(status = as.character(status)) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  mutate(across(any_of(c("G>A Direct Repair", "Nonsense Rescue",
                         "SIFT Optimization", "Untreatable")),
                as.integer)) %>%
  transmute(
    Gene   = gene_symbol,
    Total  = total,
    `G>A`  = `G>A Direct Repair`,
    `G>A_pct`   = sprintf("%.1f%%", 100 * `G>A Direct Repair` / total),
    Nonsense    = `Nonsense Rescue`,
    Nonsense_pct= sprintf("%.1f%%", 100 * `Nonsense Rescue`   / total),
    SIFT        = `SIFT Optimization`,
    SIFT_pct    = sprintf("%.1f%%", 100 * `SIFT Optimization` / total),
    Untreatable = `Untreatable`,
    Untreatable_pct = sprintf("%.1f%%", 100 * Untreatable      / total),
    Any_Treatable   = `G>A Direct Repair` + `Nonsense Rescue` + `SIFT Optimization`,
    Any_Treatable_pct = sprintf("%.1f%%",
                                100 * (`G>A Direct Repair` + `Nonsense Rescue`
                                       + `SIFT Optimization`) / total)
  ) %>%
  arrange(desc(Total))

print(wide)

# CSV
csv_path <- file.path(INT_DIR, "Top20_treatability_breakdown.csv")
write.csv(wide, csv_path, row.names = FALSE)
cat("\nwrote:", csv_path, "\n")

# Markdown table
md_path <- file.path(NOTES, "Top20_treatability_breakdown.md")
md <- c(
  "# Top-20 SFARI-1 treatability breakdown",
  paste0("Generated: ", Sys.time()),
  "",
  "| Gene | Total | G>A Direct Repair | Nonsense Rescue | SIFT Optimization | Untreatable | Any Treatable |",
  "|------|------:|-----------------:|----------------:|------------------:|-----------:|--------------:|"
)
for (i in seq_len(nrow(wide))) {
  r <- wide[i, ]
  md <- c(md, sprintf(
    "| **%s** | %d | %d (%s) | %d (%s) | %d (%s) | %d (%s) | %d (%s) |",
    r$Gene, r$Total,
    r$`G>A`,       r$`G>A_pct`,
    r$Nonsense,    r$Nonsense_pct,
    r$SIFT,        r$SIFT_pct,
    r$Untreatable, r$Untreatable_pct,
    r$Any_Treatable, r$Any_Treatable_pct
  ))
}
# Aggregate row
tot_any_all <- sum(wide$Any_Treatable)
tot_all     <- sum(wide$Total)
md <- c(md, "",
        sprintf("**Overall (top-20)**: %d / %d = %.1f%% amenable to some ADAR-mediated approach.",
                tot_any_all, tot_all, 100 * tot_any_all / tot_all))
writeLines(md, md_path)
cat("wrote:", md_path, "\n")
