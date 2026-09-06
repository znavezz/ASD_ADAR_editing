root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# MissenseOptimization — FULL variant (no SIFT<0.05, no improves=true filters).
# Everything with (a) a SIFT score and (b) a best neighbor-edit, so we can
# plot all quadrants and see where each variant lands.
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr); library(tibble)
})

VER1     <- asd_source_v1()
VER2     <- asd_source()
INT_DIR <- asd_intermediate_dir()
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)

HASURA_URL <- "http://localhost:8789/v1/graphql"
ENV_FILE   <- asd_env_file()

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
stopifnot(nchar(Sys.getenv("HASURA_ADMIN_SECRET")) > 0)

gql_read <- function(query, max_retries = 5) {
  body <- toJSON(list(query = query), auto_unbox = TRUE)
  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch(
      POST(HASURA_URL,
           add_headers(`x-hasura-admin-secret` = Sys.getenv("HASURA_ADMIN_SECRET")),
           content_type_json(), body = body, timeout(600)),
      error = function(e) { cat("[gql_read]", attempt, "err:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(resp)) { Sys.sleep(0.5 * 2^(attempt-1)); next }
    if (status_code(resp) %in% c(502,503,504)) { Sys.sleep(0.5 * 2^(attempt-1)); next }
    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
    if (!is.null(parsed$errors))
      stop("GraphQL error: ", toJSON(parsed$errors, auto_unbox = TRUE))
    return(parsed$data)
  }
  stop("gql_read: failed after ", max_retries, " attempts")
}

# =============================================================================
# Loosened MissenseOptimization query.
# CHANGES from the original:
#   - variants_quantitative_scores : no more value:{_lt: 0.05} — any SIFT
#   - neighbor_edits               : no more improves:{_eq: true} — any best edit
# =============================================================================
q <- '
query MissenseOptimizationFull {
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
      variants_quantitative_scores: {quantitative_score: {name: {_eq: "SIFT"}}}
      neighbor_edits: {is_best: {_eq: true}}
    }
    order_by: {id: asc}
  ) {
    id
    coordinate { chr start }
    ref
    alt
    gene { symbol strand }
    variant_sift: variants_quantitative_scores(where: {quantitative_score: {name: {_eq: "SIFT"}}}) { quantitative_score { value } }
    best_edit: neighbor_edits(where: {is_best: {_eq: true}}) {
      improves
      restores_reference
      neighbor_edit_scores { source value label improvement }
    }
  }
}'

t0 <- Sys.time()
raw <- gql_read(q)$result
cat("full query elapsed:", format(round(Sys.time() - t0, 1)), "  rows:", nrow(raw), "\n")
saveRDS(raw, file.path(INT_DIR, "MissenseOptimization_full_raw.rds"))

# --- Flatten ----------------------------------------------------------------
extract_sift <- function(vs_df) {
  if (is.null(vs_df) || nrow(vs_df) == 0) return(NA_real_)
  as.numeric(vs_df$quantitative_score.value[1])
}
extract_best_sift <- function(be_df) {
  if (is.null(be_df) || nrow(be_df) == 0) return(NA_real_)
  sc <- be_df$neighbor_edit_scores[[1]]
  if (is.null(sc) || nrow(sc) == 0) return(NA_real_)
  sift <- sc[sc$source == "SIFT", , drop = FALSE]
  if (nrow(sift) == 0) return(NA_real_)
  as.numeric(sift$value[1])
}
extract_best_label <- function(be_df) {
  if (is.null(be_df) || nrow(be_df) == 0) return(NA_character_)
  sc <- be_df$neighbor_edit_scores[[1]]
  if (is.null(sc) || nrow(sc) == 0) return(NA_character_)
  sift <- sc[sc$source == "SIFT", , drop = FALSE]
  if (nrow(sift) == 0) return(NA_character_)
  as.character(sift$label[1])
}
extract_improves <- function(be_df) {
  if (is.null(be_df) || nrow(be_df) == 0) return(NA)
  as.logical(be_df$improves[1])
}

df <- tibble(
  variant_id   = raw$id,
  gene         = raw$gene.symbol,
  strand       = raw$gene.strand,
  chr          = raw$coordinate.chr,
  pos          = raw$coordinate.start,
  variant_sift = vapply(raw$variant_sift, extract_sift,       numeric(1)),
  post_sift    = vapply(raw$best_edit,   extract_best_sift,  numeric(1)),
  post_label   = vapply(raw$best_edit,   extract_best_label, character(1)),
  improves     = vapply(raw$best_edit,   extract_improves,   logical(1))
)

cat("\nrows :", nrow(df),
    "\nwith variant_sift  :", sum(!is.na(df$variant_sift)),
    "\nwith post_sift     :", sum(!is.na(df$post_sift)),
    "\nwith BOTH scores   :", sum(!is.na(df$variant_sift) & !is.na(df$post_sift)),
    "\nimproves=TRUE      :", sum(df$improves, na.rm = TRUE),
    "\nimproves=FALSE     :", sum(!df$improves, na.rm = TRUE), "\n")

# Quadrant breakdown (variant_sift on x, post_sift on y; 0.05 threshold both)
d <- df %>% filter(!is.na(variant_sift), !is.na(post_sift))
q_tbl <- d %>%
  mutate(qx = ifelse(variant_sift < 0.05, "lo(del)", "hi(tol)"),
         qy = ifelse(post_sift    >= 0.05, "hi(tol)", "lo(del)")) %>%
  count(qx, qy, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))
cat("\nQuadrant breakdown (all variants with both scores):\n")
print(q_tbl)

saveRDS(df, file.path(INT_DIR, "MissenseOptimization_full_flat.rds"))
write.csv(df, file.path(INT_DIR, "MissenseOptimization_full_flat.csv"),
          row.names = FALSE)
cat("\nwrote:", file.path(INT_DIR, "MissenseOptimization_full_flat.rds"), "\n")
