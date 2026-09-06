root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Figure 3/4 combined (ver-2) — test the MissenseOptimization + RescueVariants
# queries against the live DB (port 8789).
#
# Reports row counts, gene / consequence breakdowns, SIFT distributions,
# best_edit / rescue_sift coverage.  Saves raw pulls as RDS for downstream
# panel scripts and writes a Markdown summary.
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(tibble)
})

VER1     <- asd_source_v1()
VER2     <- asd_source()
INT_DIR <- asd_intermediate_dir()
NOTES    <- file.path(VER2, "Notes")
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(NOTES,   showWarnings = FALSE, recursive = TRUE)

HASURA_URL <- "http://localhost:8789/v1/graphql"
# 2026-05-11: DB rotated to port 8789 with a new secret.  ASD_paper/.env.local
# holds the current one; check.env still has the older secret and 403s.
ENV_FILE     <- asd_env_file()
ENV_FALLBACK <- asd_env_file()

# --- load HASURA_ADMIN_SECRET without echoing it ----------------------------
.load_hasura_secret <- function(env_file) {
  if (Sys.getenv("HASURA_ADMIN_SECRET") != "") return(invisible(NULL))
  if (!file.exists(env_file)) return(invisible(NULL))
  for (l in readLines(env_file, warn = FALSE)) {
    l <- trimws(l)
    if (l == "" || startsWith(l, "#")) next
    m <- regmatches(l, regexec("^HASURA_ADMIN_SECRET=\"?([^\"]*)\"?$", l))[[1]]
    if (length(m) == 2) { Sys.setenv(HASURA_ADMIN_SECRET = m[2]); break }
  }
}
.load_hasura_secret(ENV_FILE)
if (Sys.getenv("HASURA_ADMIN_SECRET") == "") .load_hasura_secret(ENV_FALLBACK)
stopifnot(nchar(Sys.getenv("HASURA_ADMIN_SECRET")) > 0)

gql_read <- function(query, max_retries = 5) {
  q <- trimws(query)
  if (!grepl("^(query\\b|\\{)", q))
    stop("gql_read: query must start with 'query' or '{'.")
  ql <- tolower(q)
  for (tok in c("mutation ", "mutation{", "insert_", "update_", "delete_")) {
    if (grepl(tok, ql, fixed = TRUE))
      stop("gql_read: forbidden token '", tok, "' detected.")
  }
  body <- toJSON(list(query = query), auto_unbox = TRUE)
  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch(
      POST(HASURA_URL,
           add_headers(`x-hasura-admin-secret` = Sys.getenv("HASURA_ADMIN_SECRET")),
           content_type_json(), body = body, timeout(600)),
      error = function(e) { cat("[gql_read] attempt", attempt, "POST error:",
                                conditionMessage(e), "\n"); NULL }
    )
    if (is.null(resp))                         { Sys.sleep(0.5 * 2^(attempt-1)); next }
    if (status_code(resp) %in% c(502,503,504)) { Sys.sleep(0.5 * 2^(attempt-1)); next }
    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
    if (!is.null(parsed$errors))
      stop("GraphQL error: ", toJSON(parsed$errors, auto_unbox = TRUE))
    return(parsed$data)
  }
  stop("gql_read: failed after ", max_retries, " attempts")
}

# =============================================================================
# QUERY 1 — MissenseOptimization (the new panel — neighbor rescue on missense)
# =============================================================================
q_missense <- '
query MissenseOptimization {
  result: variants(
    where: {
      class: {_eq: "SNV"}
      _not: {_or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
      ]}
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_eq: "Missense"}}}
      variants_quantitative_scores: {quantitative_score: {name: {_eq: "SIFT"}, value: {_lt: 0.05}}}
      neighbor_edits: {is_best: {_eq: true}, improves: {_eq: true}}
    }
    order_by: {id: asc}
  ) {
    id
    coordinate { chr start end }
    ref
    alt
    position_in_codon
    cds_position
    ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
    alt_codon { nt1 nt2 nt3 amino_acid { short_name } }
    variants_consequences { consequence { name } }
    gene {
      ensg symbol strand
      genes_quantitative_scores(where: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}) { quantitative_score { value } }
      genes_qualitative_scores(where: {qualitative_score: {name: {_eq: "SFARI Syndromic"}}}) { qualitative_score { value } }
    }
    variant_sift: variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "SIFT"}}}) { quantitative_score { value } }
    best_edit: neighbor_edits(where: {is_best: {_eq: true}}) {
      edited_positions
      edited_cds_positions
      restores_reference
      improves
      pre_codon  { nt1 nt2 nt3 amino_acid { short_name } }
      post_codon { nt1 nt2 nt3 amino_acid { short_name } }
      neighbor_edit_scores { source value label improvement }
      neighbor_edits_consequences { consequence { name } }
    }
    variants_guides(where: {edit_type: {_eq: "Improve"}, guide: {type: {_eq: "pre_mRNA"}}}) {
      guide { sequence hits_85 hits_90 hits_95 hits_100 }
      bystanders {
        id genomic_position cds_position in_cds
        sift_score sift_prediction cadd_phred cadd_raw
        bystanders_consequences { consequence { name } }
        ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
        alt_codon { nt1 nt2 nt3 amino_acid { short_name } }
      }
    }
  }
}'

# =============================================================================
# QUERY 2 — RescueVariants (nonsense rescue, but strict: hits_85=0)
# =============================================================================
q_rescue <- '
query RescueVariants {
  result: variants(
    where: {
      class: {_eq: "SNV"}
      _not: {_or: [
        {ref: {_eq: "G"}, alt: {_eq: "A"}, gene: {strand: {_eq: "+"}}}
        {ref: {_eq: "C"}, alt: {_eq: "T"}, gene: {strand: {_eq: "-"}}}
      ]}
      variants_features: {feature: {biotype: {name: {_eq: "protein_coding"}}}}
      gene: {genes_quantitative_scores: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}}
      variants_consequences: {consequence: {name: {_eq: "StopGained"}}}
      variants_guides: {edit_type: {_eq: "Rescue"}, guide: {type: {_eq: "pre_mRNA"}, hits_85: {_eq: 0}}}
    }
    order_by: {id: asc}
  ) {
    id
    coordinate { chr start end }
    ref
    alt
    cds_position
    nmd_escaping_variant
    ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
    alt_codon { nt1 nt2 nt3 amino_acid { short_name } }
    variants_consequences { consequence { name } }
    gene {
      ensg symbol strand
      genes_quantitative_scores(where: {quantitative_score: {name: {_eq: "SFARI Gene Score"}}}) { quantitative_score { value } }
      genes_qualitative_scores(where: {qualitative_score: {name: {_eq: "SFARI Syndromic"}}}) { qualitative_score { value } }
    }
    variants_quantitative_scores(where: {quantitative_score: {name: {_in: ["SIFT", "CADD_phred"]}}}) { quantitative_score { name value } }
    variants_guides(where: {edit_type: {_eq: "Rescue"}, guide: {type: {_eq: "pre_mRNA"}}}) {
      rescue_sift
      guide { sequence hits_85 hits_90 hits_95 hits_100 }
      bystanders {
        id genomic_position cds_position in_cds
        sift_score sift_prediction cadd_phred
        bystanders_consequences { consequence { name } }
        ref_codon { nt1 nt2 nt3 amino_acid { short_name } }
        alt_codon { nt1 nt2 nt3 amino_acid { short_name } }
      }
    }
  }
}'

# =============================================================================
# RUN + save raw dumps
# =============================================================================
cat("\n--- MissenseOptimization ---\n")
t0 <- Sys.time()
m_raw <- gql_read(q_missense)$result
cat("elapsed:", format(round(Sys.time() - t0, 1)), " rows:", nrow(m_raw), "\n")
saveRDS(m_raw, file.path(INT_DIR, "MissenseOptimization_raw.rds"))

cat("\n--- RescueVariants ---\n")
t0 <- Sys.time()
r_raw <- gql_read(q_rescue)$result
cat("elapsed:", format(round(Sys.time() - t0, 1)), " rows:", nrow(r_raw), "\n")
saveRDS(r_raw, file.path(INT_DIR, "RescueVariants_raw.rds"))

# =============================================================================
# Missense: flatten SIFT + neighbor-edit best score
# =============================================================================
extract_variant_sift <- function(vs_df) {
  if (is.null(vs_df) || (is.data.frame(vs_df) && nrow(vs_df) == 0)) return(NA_real_)
  as.numeric(vs_df$quantitative_score.value[1])
}
pick_best_edit_row <- function(be_df) {
  if (is.null(be_df) || nrow(be_df) == 0) return(NULL)
  be_df[1, , drop = FALSE]
}
extract_neighbor_score <- function(be_df, label_re = NULL) {
  be <- pick_best_edit_row(be_df)
  if (is.null(be)) return(NA_real_)
  sc <- be$neighbor_edit_scores[[1]]
  if (is.null(sc) || nrow(sc) == 0) return(NA_real_)
  if (!is.null(label_re)) {
    hit <- grep(label_re, tolower(sc$label))
    if (length(hit) == 0) return(NA_real_)
    return(as.numeric(sc$value[hit[1]]))
  }
  as.numeric(sc$value[1])
}
extract_be_improvement <- function(be_df) {
  be <- pick_best_edit_row(be_df)
  if (is.null(be)) return(NA_real_)
  sc <- be$neighbor_edit_scores[[1]]
  if (is.null(sc) || nrow(sc) == 0) return(NA_real_)
  as.numeric(sc$improvement[1])
}

# Peek at score labels used in best_edit.neighbor_edit_scores so we know
# what we're dealing with before flattening columns.
label_universe <- character(0)
for (i in seq_len(nrow(m_raw))) {
  be <- m_raw$best_edit[[i]]
  if (is.null(be) || nrow(be) == 0) next
  sc <- be$neighbor_edit_scores[[1]]
  if (is.null(sc) || nrow(sc) == 0) next
  label_universe <- c(label_universe, sc$label)
}
label_universe <- sort(unique(label_universe))
cat("\nunique neighbor_edit_score labels seen:\n"); print(label_universe)

# Flatten to a scatter-ready dataframe
missense_df <- tibble(
  variant_id      = m_raw$id,
  gene            = m_raw$gene.symbol,
  strand          = m_raw$gene.strand,
  chr             = m_raw$coordinate.chr,
  pos             = m_raw$coordinate.start,
  ref             = m_raw$ref,
  alt             = m_raw$alt,
  variant_sift    = vapply(m_raw$variant_sift, extract_variant_sift, numeric(1)),
  best_edit_score_first    = vapply(m_raw$best_edit, extract_neighbor_score, numeric(1)),
  best_edit_improvement    = vapply(m_raw$best_edit, extract_be_improvement, numeric(1)),
  best_edit_score_sift     = vapply(m_raw$best_edit, extract_neighbor_score, numeric(1),
                                    label_re = "sift"),
  n_improve_guides         = vapply(m_raw$variants_guides,
                                    function(g) if (is.null(g)) 0L else as.integer(nrow(g)),
                                    integer(1))
)

cat("\nMissense flattened dataframe head:\n"); print(head(missense_df, 6))
cat("\nvariant_sift summary:\n"); print(summary(missense_df$variant_sift))
cat("\nbest_edit_score_sift summary:\n"); print(summary(missense_df$best_edit_score_sift))
cat("\nbest_edit_improvement summary:\n"); print(summary(missense_df$best_edit_improvement))

# Basic pass-quadrant candidates (mirror ver-1's CADD≥25 & rescue_sift≥0.05
# threshold but on SIFT-axis: variant SIFT < 0.05 already required in query;
# candidate quadrant is "best_edit_score_sift ≥ 0.05" — neighbor edit lifts
# the SIFT into tolerated range).
n_pass_sift <- sum(!is.na(missense_df$best_edit_score_sift) &
                     missense_df$best_edit_score_sift >= 0.05, na.rm = TRUE)
cat("\ncandidate pass quadrant (best-edit SIFT ≥ 0.05): ", n_pass_sift, "/",
    nrow(missense_df), " (",
    round(100 * n_pass_sift / nrow(missense_df), 1), "%)\n", sep = "")

saveRDS(missense_df, file.path(INT_DIR, "MissenseOptimization_flat.rds"))
write.csv(missense_df, file.path(INT_DIR, "MissenseOptimization_flat.csv"),
          row.names = FALSE)

# =============================================================================
# Rescue: flatten SIFT / CADD / rescue_sift for comparison against ver-1
# =============================================================================
extract_qscore <- function(vs_df, name) {
  if (is.null(vs_df) || nrow(vs_df) == 0) return(NA_real_)
  hit <- vs_df$quantitative_score.name == name
  if (!any(hit)) return(NA_real_)
  as.numeric(vs_df$quantitative_score.value[hit][1])
}
extract_rescue_sift <- function(vg_df) {
  if (is.null(vg_df) || nrow(vg_df) == 0) return(NA_real_)
  as.numeric(vg_df$rescue_sift[1])
}

rescue_df <- tibble(
  variant_id  = r_raw$id,
  gene        = r_raw$gene.symbol,
  strand      = r_raw$gene.strand,
  chr         = r_raw$coordinate.chr,
  pos         = r_raw$coordinate.start,
  variant_sift= vapply(r_raw$variants_quantitative_scores, extract_qscore,
                       numeric(1), name = "SIFT"),
  cadd_phred  = vapply(r_raw$variants_quantitative_scores, extract_qscore,
                       numeric(1), name = "CADD_phred"),
  rescue_sift = vapply(r_raw$variants_guides, extract_rescue_sift, numeric(1))
)

cat("\nRescue flattened dataframe head:\n"); print(head(rescue_df, 6))
cat("\nvariant_sift summary:\n"); print(summary(rescue_df$variant_sift))
cat("\ncadd_phred summary:\n");    print(summary(rescue_df$cadd_phred))
cat("\nrescue_sift summary:\n");   print(summary(rescue_df$rescue_sift))

n_pass_rescue <- sum(!is.na(rescue_df$rescue_sift) & rescue_df$rescue_sift >= 0.05 &
                     !is.na(rescue_df$cadd_phred)  & rescue_df$cadd_phred  >= 25,
                     na.rm = TRUE)
cat("\ncandidate pass quadrant (CADD≥25 & rescue_sift≥0.05): ", n_pass_rescue, "/",
    nrow(rescue_df), " (",
    round(100 * n_pass_rescue / nrow(rescue_df), 1), "%)\n", sep = "")

saveRDS(rescue_df, file.path(INT_DIR, "RescueVariants_flat.rds"))
write.csv(rescue_df, file.path(INT_DIR, "RescueVariants_flat.csv"),
          row.names = FALSE)

# =============================================================================
# Write a Markdown summary
# =============================================================================
md_lines <- c(
  "# Query test — MissenseOptimization + RescueVariants",
  paste0("Run: ", Sys.time()),
  paste0("Endpoint: `", HASURA_URL, "`"),
  "",
  "## MissenseOptimization",
  paste0("- Rows returned: **", nrow(m_raw), "**"),
  paste0("- Unique genes: ", length(unique(m_raw$gene.symbol))),
  paste0("- Variants with non-NA variant SIFT: ",
         sum(!is.na(missense_df$variant_sift))),
  paste0("- Variants with a best-edit SIFT score: ",
         sum(!is.na(missense_df$best_edit_score_sift))),
  paste0("- Variants where best-edit lifts SIFT ≥ 0.05: ",
         n_pass_sift, " (",
         round(100 * n_pass_sift / nrow(missense_df), 1), "%)"),
  paste0("- Unique `neighbor_edit_scores.label` values: ",
         paste(label_universe, collapse = ", ")),
  "",
  "### Missense SIFT distributions",
  "```",
  paste(capture.output(summary(missense_df$variant_sift)),        collapse = "\n"),
  paste(capture.output(summary(missense_df$best_edit_score_sift)), collapse = "\n"),
  "```",
  "",
  "## RescueVariants",
  paste0("- Rows returned: **", nrow(r_raw), "**"),
  paste0("- Unique genes: ", length(unique(r_raw$gene.symbol))),
  paste0("- Variants with non-NA rescue_sift: ",
         sum(!is.na(rescue_df$rescue_sift))),
  paste0("- Variants with non-NA CADD_phred: ",
         sum(!is.na(rescue_df$cadd_phred))),
  paste0("- Candidate pass (CADD≥25 & rescue_sift≥0.05): ",
         n_pass_rescue, " (",
         round(100 * n_pass_rescue / nrow(rescue_df), 1), "%)"),
  "",
  "### Rescue SIFT / CADD distributions",
  "```",
  paste(capture.output(summary(rescue_df$variant_sift)), collapse = "\n"),
  paste(capture.output(summary(rescue_df$cadd_phred)),   collapse = "\n"),
  paste(capture.output(summary(rescue_df$rescue_sift)),  collapse = "\n"),
  "```",
  "",
  "## Files written",
  paste0("- ", file.path(INT_DIR, "MissenseOptimization_raw.rds")),
  paste0("- ", file.path(INT_DIR, "MissenseOptimization_flat.rds")),
  paste0("- ", file.path(INT_DIR, "MissenseOptimization_flat.csv")),
  paste0("- ", file.path(INT_DIR, "RescueVariants_raw.rds")),
  paste0("- ", file.path(INT_DIR, "RescueVariants_flat.rds")),
  paste0("- ", file.path(INT_DIR, "RescueVariants_flat.csv"))
)
writeLines(md_lines, file.path(NOTES, "query_test_2026_07_08.md"))
cat("\nWrote summary:\n  ", file.path(NOTES, "query_test_2026_07_08.md"), "\n")
