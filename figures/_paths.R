# =============================================================================
# _paths.R — the one place the figure tree learns where anything lives.
#
# Every figure script used to carry absolute paths into one machine's home
# directory, so the tree ran on exactly one host and leaked that layout into
# the repository. Paths now resolve from two environment variables and fail
# loudly when they cannot, rather than silently writing somewhere unexpected.
#
#   ASD_PAPER_ROOT   the repository root (the directory holding figures/ and Output/)
#   ASD_PAPER_OUT    where results are written        [default: $ASD_PAPER_ROOT/Output]
#   ASD_ENV_FILE     file holding HASURA_ADMIN_SECRET [default: $ASD_PAPER_ROOT/.env.local]
#
# Only the first is required. Source this file as:
#
#   root <- Sys.getenv("ASD_PAPER_ROOT")
#   if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
#   source(file.path(root, "figures", "_paths.R"))
#
# The three-line preamble is repeated rather than factored out because there is
# nowhere to factor it *to*: locating the helper is the problem it solves.
# =============================================================================

ASD_ROOT <- local({
  root <- Sys.getenv("ASD_PAPER_ROOT")
  if (!nzchar(root)) {
    stop("ASD_PAPER_ROOT is unset. Set it to the repository root:\n",
         "  export ASD_PAPER_ROOT=/path/to/repo\n",
         "See figures/README.md.", call. = FALSE)
  }
  normalizePath(root, mustWork = TRUE)
})

ASD_OUT <- local({
  out <- Sys.getenv("ASD_PAPER_OUT")
  normalizePath(if (nzchar(out)) out else file.path(ASD_ROOT, "Output"),
                mustWork = FALSE)
})

#' A path under the repository's figures/ tree.
asd_code <- function(...) file.path(ASD_ROOT, "figures", ...)

#' A path under the output tree. Creates nothing; see `asd_out_dir()`.
asd_out <- function(...) file.path(ASD_OUT, ...)

#' A path under the output tree, with the directory created.
#' Panels write their intermediates and PDFs on first use, so the alternative is
#' a `dir.create` beside every `file.path` call.
asd_out_dir <- function(...) {
  p <- asd_out(...)
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
  p
}

#' The file holding HASURA_ADMIN_SECRET.
#'
#' Never committed: it is a credential. Returns the first candidate that exists,
#' because the repository has used both names -- `.env.example` documents `.env`,
#' while `_run_everything.R` asserts `.env.local`. Falls back to the last
#' candidate when none exists, so the caller reports a path rather than "".
#'
#' Callers must tolerate the file being absent: the panels only need a secret on
#' a regenerating run, and a cached run should not demand one it will not use.
asd_env_file <- function() {
  f <- Sys.getenv("ASD_ENV_FILE")
  if (nzchar(f)) return(normalizePath(f, mustWork = FALSE))
  for (name in c(".env.local", ".env")) {
    p <- file.path(ASD_ROOT, name)
    if (file.exists(p)) return(normalizePath(p))
  }
  file.path(ASD_ROOT, ".env.local")
}

# --- Output layout ----------------------------------------------------------
# One scheme for the whole tree, so a reader can find a figure without knowing
# which script wrote it. Composed figures are numbered as the manuscript numbers
# them; panel intermediates are not, because their filenames carry the v1 panel
# lettering and renaming them would orphan every cached .rds (see figures/README.md).
#
#   Figures/Figure_<n>/     the composed figure, as .pdf/.png/.rds
#   Figures/Supplementary/  supplementary figures
#   Figures/panels/         individual panels, for review
#   Intermediate_tables/    cached .rds and Source_data_*.csv per panel
#   Supp_Tables/            per-figure supplementary tables

#' Directory for a composed manuscript figure, e.g. `asd_figure_dir(2)`.
asd_figure_dir <- function(n) asd_out_dir("Figures")

#' Where individual panel PDFs/PNGs go.
asd_panels_dir <- function() asd_out_dir("Figures", "_build", "panels")

#' Where a panel's cached .rds and source-data CSV go.
asd_intermediate_dir <- function(name = "panels") {
  asd_out_dir("Intermediate_tables", name)
}

#' Where a figure's supplementary tables go.
asd_supp_tables_dir <- function(name) asd_out_dir("Supp_Tables", name)

#' The database snapshot every panel reads.
asd_db_snapshot <- function(...) asd_out("Intermediate_tables", "00_db_snapshot", ...)

#' Use a caller-supplied override when there is one, else the default.
#'
#' The rmds keep their path `params` so `_run_everything.R` can still override
#' them, but the YAML defaults are now empty rather than one machine's layout.
asd_param_path <- function(value, default) {
  if (is.null(value) || !nzchar(value)) default else normalizePath(value, mustWork = FALSE)
}

# --- Recovered source artifacts -------------------------------------------
# The published figures were produced in a separate working tree, whose data
# artifacts are recovered under <out>/source/. `asd_source()` is that tree;
# `asd_source_v1()` is the earlier tree it in turn reads from. Both are inputs
# only -- nothing regenerates into them, so the published composites they carry
# stay available as the comparison baseline.
asd_source <- function(...) asd_out("source", ...)

asd_source_v1 <- function(...) asd_out("source", "ver1", ...)

#' A panel intermediate: the locally regenerated one if it exists, else the
#' recovered source artifact.
#'
#' The migration is incremental -- some panels have a regeneration script in the
#' migrated tree and some are inherited from the earlier one. This resolves to
#' whichever exists and *announces the fallback*, so a run states plainly which
#' panels it rebuilt and which it merely reused. A composite that silently mixes
#' the two is how a figure comes to disagree with the code that claims to draw it.
asd_panel <- function(name, recovered_dir) {
  fresh <- file.path(asd_intermediate_dir(), name)
  if (file.exists(fresh)) return(fresh)
  message("  [inherited] ", name, " -- no local regeneration yet")
  file.path(recovered_dir, name)
}

# --- The manuscript view ---------------------------------------------------
# `Figures/` is what a reader opens: one file per manuscript figure, named as the
# manuscript names it, and nothing else. Panels, intermediates and the recovered
# source artifacts all stay on disk -- reproducibility needs them -- but under
# `Figures/_build/` so they do not clutter the thing being looked at.
#
# The separation is deliberate: a directory that mixes "the figure" with forty
# intermediates is one where nobody can tell which file went into the paper.

#' The clean, manuscript-facing figure directory.
asd_manuscript_dir <- function() asd_out_dir("Figures")

#' Build artifacts: panels, per-panel PDFs, caches. Present, but not in the way.
asd_build_dir <- function(...) asd_out_dir("Figures", "_build", ...)
