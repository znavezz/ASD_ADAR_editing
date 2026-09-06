# =============================================================================
# Spectral 8-color palette — ASD_paper ver-2.
# Requested 2026-07-08 for Figure 1 panels B, D, E.
# =============================================================================
# All eight colours, in Spectral order (red → orange → yellow → green → blue):
fig_palette <- c(
  red         = "#D53E4F",
  orange      = "#F46D43",
  light_orange= "#FDAE61",
  yellow      = "#FEE08B",
  pale_green  = "#E6F598",
  light_green = "#ABDDA4",
  teal_green  = "#66C2A5",
  blue        = "#3288BD"
)

# --- Panel B (ClinVar × SFARI donut, 3×3 grid) -------------------------------
# Two slices per donut: G>A (ADAR-treatable, highlighted) vs other.
# G>A gets the strong red side; "other" gets the soft yellow.
b_fill <- c(
  "G>A"   = unname(fig_palette["red"]),      # #D53E4F
  "other" = unname(fig_palette["yellow"])    # #FEE08B
)

# --- Panels D & E (off-target / bystander bin gradient, 6 bins) --------------
# Bin 0 = safe / desirable → teal-green.
# Bin 5+ = many hits / dangerous → red.
# Uses 6 of the 8 palette colours evenly (skips FDAE61 and 3288BD).
#
# NOTE 2026-07-08 (per the collaborators): the low end reads a bit off — the two green shades
# (0/1) and the pale-yellow (2) feel too close in luminance.  Leaving as-is
# for now; may swap in FDAE61 or reshuffle the bins later.
bin_gradient_spectral <- c(
  "0"  = unname(fig_palette["teal_green"]),  # safest  #66C2A5
  "1"  = unname(fig_palette["light_green"]), #         #ABDDA4
  "2"  = unname(fig_palette["pale_green"]),  #         #E6F598
  "3"  = unname(fig_palette["yellow"]),      #         #FEE08B
  "4"  = unname(fig_palette["orange"]),      #         #F46D43
  "5+" = unname(fig_palette["red"])          # worst   #D53E4F
)
