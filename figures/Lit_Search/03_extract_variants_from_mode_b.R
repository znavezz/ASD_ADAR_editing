root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Cross-reference the 196-variant ver-2 NEW queue against ver-1's Mode B
# (gene-level PubMed-extract) audits.
#
# For each of the 196 new variants, scan every mode_b_audit .txt for that
# gene and check whether the variant name (HGVS protein, coord, rsID) appears
# in the paper's `OTHER_GENE_VARIANTS_IN_PAPER` field or elsewhere in the text.
# Report coverage + a per-variant PMID list so we know what still needs fresh
# lit-search vs. what already has publication anchors.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble)
})

VER1        <- asd_source_v1()
VER2        <- asd_source()
INT_LIT_V2  <- asd_intermediate_dir()
MODE_B_DIR  <- file.path(VER1, "Results_Plots_Tables/Lit_Search/Correctable/mode_b_audits")

queue_new <- readRDS(file.path(INT_LIT_V2, "queue_top20_NEW_only.rds"))
cat("NEW queue:", nrow(queue_new), "variants\n")

# --- Build variant search patterns per row ---------------------------------
# The queue uses FULL amino acid names ("Arginine", "Serine"). Convert to
# single-letter + three-letter forms and build every citation pattern we've
# ever seen in ver-1 audits.
full_to_one <- c(Alanine="A", Arginine="R", Asparagine="N", Aspartate="D",
                 "Aspartic acid"="D", Cysteine="C", Glutamate="E",
                 "Glutamic acid"="E", Glutamine="Q", Glycine="G",
                 Histidine="H", Isoleucine="I", Leucine="L", Lysine="K",
                 Methionine="M", Phenylalanine="F", Proline="P",
                 Serine="S", Threonine="T", Tryptophan="W", Tyrosine="Y",
                 Valine="V", Stop="*", Ter="*")
one_to_three <- c(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys",
                  E="Glu", Q="Gln", G="Gly", H="His", I="Ile",
                  L="Leu", K="Lys", M="Met", F="Phe", P="Pro",
                  S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val",
                  "*"="Ter")

to_one <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  if (x %in% names(full_to_one)) return(unname(full_to_one[x]))
  if (nchar(x) == 3 && x %in% names(one_to_three) == FALSE) {
    # It might already be a 3-letter code
    m <- match(x, one_to_three)
    if (!is.na(m)) return(names(one_to_three)[m])
  }
  if (nchar(x) == 1) return(x)
  NA_character_
}

build_patterns <- function(row) {
  ref1 <- to_one(row$ref_aa_name)
  # Detect nonsense variant — alt_aa_name is NA in data for StopGained
  is_stop <- row$treat_class == "Nonsense Rescue" ||
             (!is.na(row$alt_aa_name) && row$alt_aa_name %in% c("Stop", "Ter", "*"))
  alt1 <- if (is_stop) "*" else to_one(row$alt_aa_name)
  pos  <- row$protein_position

  patterns <- character(0)
  if (!is.na(pos) && !is.na(ref1)) {
    ref3 <- unname(one_to_three[ref1])
    if (is_stop) {
      # Nonsense patterns
      patterns <- c(patterns,
                    paste0("p\\.", ref1, pos, "\\*"),
                    paste0(ref1, pos, "\\*"),
                    paste0("p\\.", ref1, pos, "X"),
                    paste0(ref1, pos, "X"),
                    paste0("p\\.", ref3, pos, "Ter"),
                    paste0(ref3, pos, "Ter"),
                    paste0(ref3, pos, "\\*"))
    } else if (!is.na(alt1)) {
      alt3 <- unname(one_to_three[alt1])
      # Missense patterns
      patterns <- c(patterns,
                    paste0("p\\.", ref1, pos, alt1),
                    paste0(ref1, pos, alt1),
                    if (!is.na(alt3))
                      c(paste0("p\\.", ref3, pos, alt3),
                        paste0(ref3, pos, alt3)))
    }
  }
  # Coord + rsID fallbacks
  patterns <- c(patterns,
                paste0("chr", row$chr, ":", row$start),
                paste0("chr", row$chr, "_", row$start))
  patterns <- patterns[!is.na(patterns) & patterns != ""]
  unique(patterns)
}

# --- Scan all mode_b_audit files -------------------------------------------
audit_files <- list.files(MODE_B_DIR, pattern = "\\.txt$", full.names = TRUE)
cat("Mode B audits:", length(audit_files), "files\n")

# Parse (gene, PMID) from filename: <GENE>__<PMID>.txt
audits <- tibble(
  path = audit_files,
  fname = basename(audit_files),
  gene  = str_extract(basename(audit_files), "^[^_]+"),
  pmid  = str_extract(basename(audit_files), "(?<=__)[0-9]+")
) %>%
  mutate(text = sapply(path, function(p) paste(readLines(p, warn = FALSE),
                                               collapse = "\n")))

# --- For each queue variant, find matching audits -------------------------
hits <- list()
for (i in seq_len(nrow(queue_new))) {
  row <- queue_new[i, ]
  gene_audits <- audits %>% filter(gene == row$gene_symbol)
  if (nrow(gene_audits) == 0) next
  patterns <- build_patterns(row)
  if (length(patterns) == 0) next

  for (j in seq_len(nrow(gene_audits))) {
    hit <- any(sapply(patterns, function(p) grepl(p, gene_audits$text[j],
                                                    ignore.case = FALSE)))
    if (hit) {
      hits[[length(hits) + 1]] <- tibble(
        variant_id  = row$id,
        gene        = row$gene_symbol,
        hgvs_p      = row$hgvs_p,
        chr_pos     = paste0("chr", row$chr, ":", row$start),
        treat_class = row$treat_class,
        pmid        = gene_audits$pmid[j],
        audit_file  = gene_audits$fname[j],
        matched_pattern = paste(patterns[sapply(patterns, function(p)
          grepl(p, gene_audits$text[j], ignore.case = FALSE))], collapse = "|")
      )
    }
  }
}
hit_tbl <- if (length(hits) > 0) bind_rows(hits) else tibble()

cat("\n--- Cross-ref results ---\n")
cat("(variant, PMID) hit pairs:", nrow(hit_tbl), "\n")
if (nrow(hit_tbl) > 0) {
  cat("Distinct variants with >=1 hit:",
      length(unique(hit_tbl$variant_id)), "\n")
  cat("Distinct PMIDs:", length(unique(hit_tbl$pmid)), "\n\n")

  cat("Variants with anchors (by class):\n")
  covered_ids <- unique(hit_tbl$variant_id)
  print(queue_new %>%
          mutate(covered = id %in% covered_ids) %>%
          count(treat_class, covered))

  cat("\nVariants WITHOUT any anchor (need fresh lit-search):\n")
  uncovered <- queue_new %>% filter(!id %in% covered_ids)
  print(uncovered %>% count(treat_class))
}

# --- Save outputs ----------------------------------------------------------
saveRDS(hit_tbl, file.path(INT_LIT_V2, "mode_b_hits_new196.rds"))
write.csv(hit_tbl, file.path(INT_LIT_V2, "mode_b_hits_new196.csv"),
          row.names = FALSE)

uncovered_variants <- if (nrow(hit_tbl) > 0) {
  queue_new %>% filter(!id %in% unique(hit_tbl$variant_id))
} else queue_new
saveRDS(uncovered_variants, file.path(INT_LIT_V2, "uncovered_after_mode_b.rds"))
write.csv(uncovered_variants,
          file.path(INT_LIT_V2, "uncovered_after_mode_b.csv"),
          row.names = FALSE)

cat("\nwrote:\n  ",
    file.path(INT_LIT_V2, "mode_b_hits_new196.rds"), "\n  ",
    file.path(INT_LIT_V2, "uncovered_after_mode_b.rds"), "\n", sep = "")
