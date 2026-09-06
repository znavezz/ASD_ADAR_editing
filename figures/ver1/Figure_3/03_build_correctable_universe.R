root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Build the "Correctable" universe — the full ADAR-treatable variant set for
# downstream lit-search and treatment-proposal analysis.
#
#   CDS group (Direct Repair) : 2,421 gene-sense G>A SFARI SNVs in
#                                {Missense, Nonsense, Splice} with a Fix guide
#                                — same scope as Figure 1D/1E
#   Rescue group              :   81 NonG2A StopGained SFARI SNVs passing
#                                CADD ≥ 25 AND rescue_sift ≥ 0.05
#                                — pass quadrant of Figure 2
#
# Total Correctable = 2,502 (no overlap by definition: 94 G>A nonsense are in
# the CDS group as Direct Repair; the 81 Rescue are non-G>A stop-gained).
#
# Output schema matches the existing 644-variant lit-search CSV so the
# triage/lit-search pipeline can plug into the larger universe with no code
# changes — only a new input filename.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(httr); library(jsonlite); library(tidyr)
})

# Writes into the repository's own output tree; the recovered ver-1 tree under
# Output/source/ is an input and is never written to.
OUT_DIR <- asd_supp_tables_dir("Correctable")

env_lines <- readLines(asd_env_file(), warn = FALSE)
secret <- sub("^HASURA_ADMIN_SECRET=\"?([^\"]+)\"?", "\\1",
              grep("^HASURA_ADMIN_SECRET=", env_lines, value = TRUE))
stopifnot(nchar(secret) > 0)
HASURA <- Sys.getenv("ASD_HASURA_URL", "http://localhost:8789/v1/graphql")
gql <- function(q) {
  r <- POST(HASURA, add_headers(`x-hasura-admin-secret` = secret),
            content_type_json(),
            body = toJSON(list(query = q), auto_unbox = TRUE), timeout(120))
  p <- fromJSON(content(r, "text", encoding = "UTF-8"), flatten = TRUE)
  if (!is.null(p$errors)) stop("GraphQL: ", toJSON(p$errors, auto_unbox = TRUE))
  p$data
}

# --- 1. Pull CDS-group variant IDs (gene-sense G>A, MS/Spl/NS, Fix guide) ----
q_cds <- '{ variants(where: {class: {_eq: "SNV"},
  _or: [
    {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}},
    {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
  ],
  variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}},
  gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}},
  variants_consequences: {consequence: {name: {_in: ["StopGained","Missense","Splice"]}}},
  variants_guides: {edit_type: {_eq: "Fix"}}
}, limit: 100000) { id } }'
cds_ids <- as.integer(gql(q_cds)$variants$id)
cat("CDS group (gene-sense G>A):", length(cds_ids), "\n")

# --- 2. Pull NonG2A_StopGained variant IDs, then filter to pass-quadrant ----
q_rescue <- '{ variants(where: {class: {_eq: "SNV"},
  _not: { _or: [
    {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}},
    {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
  ]},
  variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}},
  gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}},
  variants_consequences: {consequence: {name: {_eq: "StopGained"}}}
}, limit: 100000) { id } }'
rescue_universe_ids <- as.integer(gql(q_rescue)$variants$id)
cat("NonG2A StopGained universe:", length(rescue_universe_ids), "\n")

# Per-variant canonical CADD + rescue_sift (chunked)
chunks <- split(rescue_universe_ids, ceiling(seq_along(rescue_universe_ids) / 500))
rescue_quality <- bind_rows(lapply(chunks, function(chunk) {
  q <- sprintf('{ variants(where: {id: {_in: [%s]}}) {
    id
    variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "CADD_phred"}}}) {
      quantitative_score { value }
    }
    variants_guides(where: {edit_type: {_eq: "Rescue"}, guide: {type: {_eq: "pre_mRNA"}}}) {
      rescue_sift
    }
  } }', paste(chunk, collapse = ","))
  d <- gql(q)$variants
  data.frame(
    variant_id  = d$id,
    cadd_phred  = vapply(d$variants_quantitative_scores, function(x) {
      if (!is.null(x) && nrow(x) > 0) as.numeric(x$quantitative_score.value[1]) else NA_real_
    }, numeric(1)),
    rescue_sift = vapply(d$variants_guides, function(x) {
      if (!is.null(x) && nrow(x) > 0) as.numeric(x$rescue_sift[1]) else NA_real_
    }, numeric(1)))
}))
rescue_ids <- rescue_quality %>%
  filter(!is.na(cadd_phred), !is.na(rescue_sift),
         cadd_phred >= 25, rescue_sift >= 0.05) %>%
  pull(variant_id)
cat("Rescue group (pass quadrant):", length(rescue_ids), "\n")

# --- 3. Load df_cds (rich annotations) and join ------------------------------
df_cds <- readRDS(asd_db_snapshot("df_cds.rds"))
cat("df_cds rows:", nrow(df_cds), "\n")

cds_df <- df_cds %>%
  filter(id %in% cds_ids) %>%
  mutate(correctable_class = "Direct Repair")
rescue_df <- df_cds %>%
  filter(id %in% rescue_ids) %>%
  mutate(correctable_class = "Rescue")
cat("CDS rows matched in df_cds:   ", nrow(cds_df), "\n")
cat("Rescue rows matched in df_cds:", nrow(rescue_df), "\n")

# Sanity: should be 2421 + 81 = 2502 (any rows missing from df_cds means
# they didn't pass the 00_prep_data CDS+Splice filter — should be zero)
correctable <- bind_rows(cds_df, rescue_df)
cat("Total Correctable:", nrow(correctable),
    if (nrow(correctable) == 2502) "  ✓" else "  ⚠ expected 2502", "\n\n")

# --- 4. Build the lit-search schema ------------------------------------------
# Mirror the 644-CSV columns; add correctable_class + the Rescue-only scores
correctable <- correctable %>%
  mutate(
    treatment_status = ifelse(correctable_class == "Direct Repair",
                              "G>A (Direct Repair)",
                              "Nonsense Rescue"),
    chr_pos = paste0("chr", chr, ":", start),
    REF_genesense = ifelse(gene_strand == "+", REF.VEP,
                           chartr("ACGT","TGCA", REF.VEP)),
    ALT_genesense = ifelse(gene_strand == "+", ALT.VEP,
                           chartr("ACGT","TGCA", ALT.VEP)),
    ref_codon = paste0(ref_codon_nt1, ref_codon_nt2, ref_codon_nt3),
    alt_codon = paste0(alt_codon_nt1, alt_codon_nt2, alt_codon_nt3),
    Consequence_simple = case_when(
      str_detect(New_Consequence, "stop_gained|Nonsense") ~ "Nonsense",
      str_detect(New_Consequence, "missense|Missense") ~ "Missense",
      str_detect(New_Consequence, "splice|Splice") ~ "Splice",
      str_detect(New_Consequence, "synonymous|Synonymous") ~ "Synonymous",
      str_detect(New_Consequence, "start_lost|Start_lost") ~ "Start Lost",
      TRUE ~ New_Consequence),
    protein_change = ifelse(is.na(ref_aa) | is.na(alt_aa) |
                              ref_aa == "" | alt_aa == "",
                            "",
                            paste0("p.", ref_aa,
                                   protein_position,
                                   ifelse(alt_aa == "*", "Stop", alt_aa))),
    Consequence_VEP = New_Consequence,
    paper_keys = ifelse(is.na(paper_keys), "", paper_keys),
    paper_count = ifelse(paper_keys == "", 0L,
                         lengths(strsplit(paper_keys, ";")))
  ) %>%
  transmute(
    gene_symbol, treatment_status, correctable_class,
    chr_pos, chr, start, end,
    gene_strand, REF_genesense, ALT_genesense,
    REF_VEP = REF.VEP, ALT_VEP = ALT.VEP,
    cds_position, protein_position, exon,
    ref_codon, alt_codon, ref_aa, ref_aa_name,
    alt_aa, alt_aa_name, alt_codon_stop,
    protein_change, Consequence_VEP, Consequence_simple,
    rs_id, patient_count, paper_count, paper_keys,
    CLIN_SIG, SIFT_score, SIFT_pred,
    polyphen_score, polyphen_pred,
    REVEL, am_pathogenicity, am_class,
    CADD_phred, PhastCons100, pLI, LoFtool, BLOSUM62,
    sfari_gene_score, sfari_syndromic, sfari_reports,
    inheritance, nmd_escaping_variant,
    variant_id = id
  )

cat("=== Final breakdown ===\n")
print(correctable %>% count(correctable_class, Consequence_simple))

cat("\n=== By treatment_status ===\n")
print(correctable %>% count(treatment_status))

cat("\n=== By SFARI score ===\n")
print(correctable %>% count(sfari_gene_score))

# --- 5. Save ----------------------------------------------------------------
csv_path <- file.path(OUT_DIR, "Variants_correctable_2502.csv")
rds_path <- file.path(OUT_DIR, "Variants_correctable_2502.rds")
write.csv(correctable, csv_path, row.names = FALSE)
saveRDS(correctable, rds_path)
cat("\nWrote:\n  ", csv_path, "\n  ", rds_path, "\n", sep = "")
cat("Total rows:", nrow(correctable),
    "  Total cols:", ncol(correctable), "\n")
