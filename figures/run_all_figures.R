# =============================================================================
# run_all_figures.R — rebuild every manuscript figure, end to end.
#
#   Rscript figures/run_all_figures.R
#
# Needs no environment set: ASD_PAPER_ROOT is derived from this file's location.
# A fresh run pulls from the Hasura API; set REGEN=0 to reuse cached panel data
# where a panel supports it.
#
# Order matters. Panels are built before the composites that assemble them, and
# the two panels that feed Figure 3's bar chart depend on the missense set built
# for Figure 2, so Figure 2 runs first.
#
# Every composite reports which panels it rebuilt and which it reused. A run that
# prints "[inherited]" is telling you a panel has no regeneration path and came
# from the recovered artifacts instead -- worth knowing before trusting the
# output, because a composite that silently mixes the two is how a figure comes
# to disagree with the code that claims to draw it.
# =============================================================================

if (!nzchar(Sys.getenv("ASD_PAPER_ROOT"))) {
  a    <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", a[grep("^--file=", a)])
  if (!length(self)) {
    stop("ASD_PAPER_ROOT is unset and this script's path is unknown ",
         "(sourced rather than run with Rscript). Set it explicitly.",
         call. = FALSE)
  }
  Sys.setenv(ASD_PAPER_ROOT = dirname(dirname(normalizePath(self))))
}
source(file.path(Sys.getenv("ASD_PAPER_ROOT"), "figures", "_paths.R"))

REGEN <- toupper(Sys.getenv("REGEN", "TRUE")) %in% c("TRUE", "1", "YES")
SRC   <- file.path(ASD_ROOT, "figures")

message("ASD_PAPER_ROOT : ", ASD_ROOT)
message("output         : ", ASD_OUT)
message("REGEN          : ", REGEN, "\n")

# --- runners ---------------------------------------------------------------

run_r <- function(rel) {
  message("\n=== ", rel)
  system2("Rscript", file.path(SRC, rel), stdout = "", stderr = "")
}

#' Knit a panel rmd without pandoc: chunks execute and write their side effects,
#' which is all a panel needs. Params are the YAML defaults with `regenerate`
#' overridden, so the file stays runnable on its own.
run_rmd <- function(rel) {
  message("\n=== ", rel)
  path <- file.path(SRC, rel)
  fm   <- rmarkdown::yaml_front_matter(path)$params
  td   <- tempdir()
  knitr::opts_chunk$set(fig.path = file.path(td, paste0(basename(rel), "_")))
  env <- new.env()
  assign("params", modifyList(fm, list(regenerate = REGEN)), envir = env)
  old <- setwd(td); on.exit(setwd(old), add = TRUE)
  knitr::knit(path, output = file.path(td, paste0(basename(rel), ".md")),
              envir = env, quiet = TRUE)
}

# --- bootstrap --------------------------------------------------------------
# Pulls the master variant table from the database. Everything below reads it,
# so a fresh clone must run this first -- there is no snapshot in the repository.
run_rmd("ver1/00_prep_data.rmd")

# The missense-optimization set: the input to Figure 3B, and to the SIFT
# Optimization class that Figure 4 and Figure S1 colour by.
run_r("Figure_3/00_test_queries.R")        # writes MissenseOptimization_raw
run_r("Figure_3/00b_test_missense_full.R") # writes MissenseOptimization_full_flat

# --- panels ----------------------------------------------------------------
# Figures 1 and 2 (the split composite): two panels come from the earlier tree,
# which has no later replacement for them; three from the later one.
run_rmd("ver1/Figure_1/00_panel_A_funnel.rmd")          # Fig 2A - filtering funnel
run_rmd("ver1/Figure_1/01_panel_B_mismatch.rmd")        # Fig 2C - G>A consequence bar
run_rmd("ver1/Figure_1/02_panel_C_clinvar_donut.rmd")   # plotdata for Fig 2B
run_rmd("ver1/Figure_1/03_panel_D_offtargets.rmd")      # plotdata for Fig 2D, 2E
run_r("Figures_1_2/panel_B_clinvar_donut_ver2.R")       # Fig 2B - ClinVar x SFARI donut
run_r("Figures_1_2/panels_D_E_offtargets_ver2.R")       # Fig 2D, 2E - off-targets, bystanders

# Figure 3: the missense-optimization scatter.
run_r("Figure_3/03_panel_B_missense_scatter.R")         # Fig 3B

# Figure 4: heatmap from the earlier tree, bars and lollipop from the later one.
run_rmd("ver1/Figure_3/00_panel_brain_heatmap.rmd")     # Fig 4A, left
run_r("Figure_4/01_treatability_bar_ver2.R")            # Fig 4A, right
run_r("Figure_4/04_CHD8_lollipop_ver2.R")               # Fig 4B

# --- composites ------------------------------------------------------------
run_r("Figures_1_2/99_compose_figures_1_2.R")
run_r("Figure_3/99_compose_figure_3.R")
run_r("Figure_4/99_compose_figure_4.R")

# --- supplementary ---------------------------------------------------------
run_r("Figure_4/05_supp_all_lollipops_ver2.R")          # Figure S1

message("\nDone. Figures under: ", asd_out("Figures"))
message("Compare against the published renders with:")
message("  Output/source/Results_Plots_Tables/Composed/")
