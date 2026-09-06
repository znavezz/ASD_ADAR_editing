root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# Ver-2 Fig 3 Panel C — CHD8 domain lollipop.
#
# DESIGN (per the collaborators, 2026-07-08):
#   SHAPE + COLOUR both by treatment class (ver-1 pairing style):
#     Direct Repair  → circle,   #66C2A5 (teal green)
#     Nonsense Rescue    → triangle, #ABDDA4 (light green)
#     SIFT Optimization  → diamond,  #FEE08B (yellow)      [new for ver-2]
#     Untreatable        → square,   grey55
#   HEIGHT = CADD PHRED score
#   DOMAINS = alternating stripes using project palette (blue / orange)
# =============================================================================

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr)
  library(ggplot2); library(ggrepel); library(stringr); library(scales)
})

# ggrepel resolves label collisions with random jitter, so an unseeded run places
# domain labels differently every time. Nothing in the data moves -- only where a
# label lands and how its leader line is drawn -- but it is the reason these
# figures were not byte-reproducible. Seeded here so they are.
set.seed(42)

VER1     <- asd_source_v1()
VER2     <- asd_source()
IN_INT   <- file.path(VER2, "Intermediate_tables", "Figure_3")  # recovered: read only
INT_DIR  <- asd_intermediate_dir()                              # ours: written to
OUT_MAIN <- asd_panels_dir()
INT_MO   <- asd_intermediate_dir()
dir.create(INT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_MAIN, showWarnings = FALSE, recursive = TRUE)

source(asd_code("ver1", "_helpers.R"))
source(asd_code("_palette.R"))

GENE <- Sys.getenv("GENE", unset = "CHD8")

# --- Data + treatment-class classification ---------------------------------
df_cds       <- readRDS(asd_db_snapshot("df_cds.rds"))
ids_nonsense <- readRDS(asd_panel("nonsense_rescue_894_ids.rds", IN_INT))
ids_siftopt  <- as.integer(readRDS(file.path(INT_MO,
                                             "MissenseOptimization_raw.rds"))$id)

df <- df_cds %>%
  filter(gene_symbol == GENE,
         !is.na(protein_position), !is.na(CADD_phred)) %>%
  mutate(
    protein_position = as.numeric(protein_position),
    cons_first = sapply(str_split(New_Consequence, ","), `[`, 1),
    consequence_simple = case_when(
      cons_first == "Missense"   ~ "Missense",
      cons_first == "Splice"     ~ "Splice",
      cons_first == "Nonsense"   ~ "StopGained",
      TRUE                       ~ NA_character_
    ),
    is_G_to_A = (REF == "G" & ALT == "A"),
    treat_class = case_when(
      is_G_to_A               ~ "Direct Repair",
      id %in% ids_nonsense    ~ "Nonsense Rescue",
      id %in% ids_siftopt     ~ "SIFT Optimization",
      TRUE                    ~ "Untreatable"
    )
  ) %>%
  filter(!is.na(consequence_simple))   # drops Synonymous / Start_lost

cat("CHD8 plotting rows:", nrow(df), "\n")
cat("shape (consequence):\n"); print(table(df$consequence_simple))
cat("colour (treat_class):\n"); print(table(df$treat_class))
cat("counts by shape × colour:\n")
print(table(df$consequence_simple, df$treat_class))

# --- UniProt fetch ---------------------------------------------------------
fetch_uniprot <- function(gene_symbol) {
  url <- modify_url("https://rest.uniprot.org/uniprotkb/search",
                    query = list(
                      query = paste0("gene_exact:", gene_symbol,
                                     " AND organism_id:9606 AND reviewed:true"),
                      format = "json"))
  resp <- GET(url, timeout(20))
  parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = FALSE)
  features <- parsed$results$features[[1]]
  p_length <- as.integer(parsed$results$sequence$length[1])
  domain_df <- features %>%
    filter(type %in% c("Domain", "Region", "Zinc finger", "Motif", "Repeat")) %>%
    mutate(domain_name = ifelse(!is.na(description) & description != "",
                                 description, type),
           start = as.numeric(location$start$value),
           end   = as.numeric(location$end$value)) %>%
    select(domain_name, start, end)
  domain_df$domain_name <- gsub(" domain", "", domain_df$domain_name,
                                ignore.case = TRUE)
  list(domain_df = domain_df, p_length = p_length)
}

up <- fetch_uniprot(GENE)
cat("\nCHD8 protein length:", up$p_length, "AA\n")
cat("CHD8 domains fetched:", nrow(up$domain_df), "\n")

domains <- up$domain_df %>%
  mutate(clean_name = trimws(gsub(" [0-9]+$", "", domain_name)),
         display_name = case_when(
           grepl("disorder", clean_name, ignore.case = TRUE) ~ "Disordered Region",
           TRUE ~ clean_name),
         xmin = start, xmax = end) %>%
  arrange(xmin) %>%
  mutate(strip_col = "grey25",   # all domains same dark grey
         ymin = -5,               # y_track_bottom, hardcoded to keep the
         ymax =  0)               # ggplot RDS env-independent

label_domains <- domains %>%
  filter(display_name != "Disordered Region") %>%
  group_by(display_name) %>%
  summarise(x_mid = median((xmin + xmax) / 2), .groups = "drop") %>%
  mutate(domain_label = str_wrap(display_name, width = 15),
         y_anchor     = -5)   # y_track_bottom baked into the data so the
                              # saved ggplot doesn't reference the env

# --- Aesthetic scales ------------------------------------------------------
treat_colors <- c(
  "Direct Repair" = unname(fig_palette["teal_green"]),  # #66C2A5
  "Nonsense Rescue"   = unname(fig_palette["light_green"]), # #ABDDA4
  "SIFT Optimization" = unname(fig_palette["yellow"]),      # #FEE08B
  "Untreatable"       = "grey85"   # light grey to match Panel B fill (grey90)
)
# Shape by TREATMENT CLASS (paired with colour — ver-1 style).
# Note: Direct Repair AND Untreatable both use circle now (per the collaborators, 2026-07-08).
treat_shapes <- c(
  "Direct Repair"     = 22,   # filled square
  "Nonsense Rescue"   = 24,   # filled triangle
  "SIFT Optimization" = 23,   # filled diamond
  "Untreatable"       = 21    # filled circle
)

df$treat_class <- factor(df$treat_class,
                          levels = c("Direct Repair",
                                     "Nonsense Rescue",
                                     "SIFT Optimization",
                                     "Untreatable"))

# --- Plot -----------------------------------------------------------------
y_track_top    <- 0
y_track_bottom <- -5
y_label_space  <- -25
max_y <- max(df$CADD_phred, 40, na.rm = TRUE) + 5

# Immediately before the repel call, not at the top of the file: any code
# between the two advances the RNG, so a distant seed does not pin the
# label layout. This is why an earlier seeding attempt fixed Figure S1 but
# left Figures 1 and 3 varying run to run.
set.seed(42)
p <- ggplot() +
  # domain-track background (annotate → eagerly evaluated so p_length is
  # substituted immediately, keeping the ggplot RDS free of env references)
  annotate("rect", xmin = 0, xmax = up$p_length,
           ymin = y_track_bottom, ymax = y_track_top,
           fill = "gray97", colour = "gray70", linewidth = 0.3) +
  # domain stripes — all grey25, alpha 0.6 (ymin/ymax now live in `domains`)
  geom_rect(data = domains,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin,   ymax = ymax),
            fill = domains$strip_col,
            colour = "black", alpha = 0.6, linewidth = 0.3) +
  # domain labels
  geom_text_repel(seed = 42, data = label_domains,
                  aes(x = x_mid, y = y_anchor, label = domain_label),
                  size = 3.2, fontface = "bold", colour = "black",
                  nudge_y = -8, direction = "both", vjust = 1,
                  segment.size = 0.5, segment.colour = "gray40",
                  min.segment.length = 0, max.overlaps = 50) +
  # lollipop stems
  geom_segment(data = df,
               aes(x = protein_position, xend = protein_position,
                   y = 0, yend = CADD_phred),
               colour = "gray75", linewidth = 0.5, alpha = 0.75) +
  # variant points — shape + fill both encode treatment class
  geom_point(data = df,
             aes(x = protein_position, y = CADD_phred,
                 fill = treat_class, shape = treat_class),
             colour = "black", size = 3.2, stroke = 0.4, alpha = 0.95) +
  # CADD 25 dashed threshold
  geom_hline(yintercept = 25, linetype = "dashed",
             colour = "red2", alpha = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = treat_colors, name = "Treatment class") +
  scale_shape_manual(values = treat_shapes, name = "Treatment class") +
  scale_x_continuous(limits = c(0, up$p_length),
                     expand = c(0.01, 0), labels = scales::comma) +
  coord_cartesian(ylim = c(y_label_space, max_y), clip = "off") +
  labs(title = GENE,
       subtitle = paste0("Amino acid position (total: ", up$p_length, " AA)"),
       x = NULL, y = "CADD PHRED score", tag = "C") +
  theme_publication(base_size = 12, y_bold = TRUE) +
  theme(plot.title       = element_text(face = "bold", size = 15,
                                        family = "Nimbus Sans"),
        plot.subtitle    = element_text(colour = "gray40", size = 10),
        legend.position  = "top",
        legend.box       = "horizontal",
        legend.title     = element_text(size = 10),
        legend.text      = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "gray95"),
        plot.tag         = element_text(face = "bold", size = 16,
                                        family = "Nimbus Sans"),
        plot.tag.position= c(0, 1))

panel_id <- paste0("Figure_3C_", GENE, "_lollipop_ver2")
save_panel(p, panel_id,
           pdf_dir = OUT_MAIN, rds_dir = INT_DIR,
           width = 14, height = 6)
ggsave(file.path(OUT_MAIN, paste0(panel_id, ".png")),
       p, width = 14, height = 6, dpi = 150, bg = "white")

cat("\nPanel C ver-2 (", GENE, " lollipop) written.\n", sep = "")
