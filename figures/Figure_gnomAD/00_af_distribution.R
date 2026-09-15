root <- Sys.getenv("ASD_PAPER_ROOT")
if (!nzchar(root)) stop("ASD_PAPER_ROOT is unset -- see figures/README.md")
source(file.path(root, "figures", "_paths.R"))

# =============================================================================
# gnomAD allele-frequency distribution across the three ADAR editing classes.
#
# Rarity is a MAGNITUDE, so the measured bands use one hue, light to dark.
# "Absent from gnomAD" is held out in neutral grey because it is a CATEGORY,
# not a magnitude: 43.7% of the Direct Repair set has no gnomAD row at all, and
# absence from a population reference is the strongest rarity evidence there is
# (ACMG PM2). Treating it as a low number - or filtering it away - would invert
# the biology for nearly half the data.
#
# Uses gnomade_af (gnomAD exomes r2.1), NOT MAX_AF: MAX_AF is VEP's maximum over
# 1000G/ESP/gnomAD and is not a gnomAD value. gnomADg is NULL for every variant,
# since the GRCh37 cache carries exomes only.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})
source(asd_code("ver1", "_helpers.R"))   # theme_publication

d <- read.csv(asd_out("Intermediate_tables", "gnomad_af_by_class.csv"),
              stringsAsFactors = FALSE)

# Collapse the two highest bands: >=1% is the threshold that carries meaning
# (ACMG BS1), and separating 1-5% from >=5% produced slivers too small to read.
# The >=5% count is called out in the caption instead.
d <- d %>%
  mutate(band = ifelse(band %in% c("0.01 - 0.05", ">= 0.05"), ">= 0.01", band)) %>%
  group_by(editing_class, band, total) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(pct = 100 * n / total)

BANDS <- c("Absent from gnomAD", "< 0.0001", "0.0001 - 0.001", "0.001 - 0.01", ">= 0.01")
af_fill <- c(
  "Absent from gnomAD" = "#9E9E9E",
  "< 0.0001"           = "#BDD7E7",
  "0.0001 - 0.001"     = "#6BAED6",
  "0.001 - 0.01"       = "#3182BD",
  ">= 0.01"            = "#08519C"
)
CLASS_ORDER <- c("Missense Optimization", "Nonsense Rescue", "Direct Repair")

d <- d %>%
  mutate(band = factor(band, levels = BANDS),
         editing_class = factor(editing_class, levels = CLASS_ORDER))

lab <- d %>% filter(pct >= 4)
n_by_class <- d %>% distinct(editing_class, total)

p <- ggplot(d, aes(x = pct, y = editing_class, fill = band)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.9,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = lab, aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            colour = ifelse(lab$band %in% c("0.001 - 0.01", ">= 0.01"), "white", "grey15"),
            size = 3.1, fontface = "bold") +
  scale_fill_manual(values = af_fill, name = "gnomAD exome allele frequency",
                    guide = guide_legend(nrow = 1, title.position = "top", title.hjust = 0)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = function(v) {
    tot <- n_by_class$total[match(v, as.character(n_by_class$editing_class))]
    paste0(v, "\n(n = ", trimws(format(tot, big.mark = ",")), ")")
  }) +
  labs(
    title = "gnomAD allele frequency of variants amenable to ADAR-mediated editing",
    x = "Share of variants in the class", y = NULL
  ) +
  theme_publication(base_size = 12) +
  theme(
    legend.position    = "top",
    legend.title       = element_text(size = 9.5),
    legend.text        = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    axis.text.y        = element_text(face = "bold", lineheight = 1.05),
    plot.title         = element_text(face = "bold", size = 12.5, hjust = 0),
    legend.margin      = margin(b = 4),
    axis.title.x       = element_text(size = 10, colour = "grey25")
  )

OUT <- asd_out_dir("Figures")
ggsave(file.path(OUT, "Figure_gnomAD_AF_distribution.png"), p,
       width = 9.5, height = 4.2, dpi = 300, bg = "white")
# cairo_pdf, as the rest of the suite uses: the default pdf device cannot
# embed the bold face this theme asks for and errors with "invalid font type".
ggsave(file.path(OUT, "Figure_gnomAD_AF_distribution.pdf"), p,
       width = 9.5, height = 4.2, bg = "white", device = cairo_pdf)
write.csv(d %>% arrange(editing_class, band) %>% select(editing_class, band, n, total, pct),
          file.path(OUT, "Figure_gnomAD_AF_distribution_data.csv"), row.names = FALSE)

cat("wrote to", OUT, "\n")
print(d %>% arrange(editing_class, band) %>% as.data.frame())
