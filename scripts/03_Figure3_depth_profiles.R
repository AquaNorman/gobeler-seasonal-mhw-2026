## ---------------------------------------------------------------------------
## Figure 3: Organic matter content and porewater nutrient depth profiles
##
## Input:
##   - PW_OM_nutrients_depth_profile.csv
##
## Output:
##   - Figure3.tiff
## ---------------------------------------------------------------------------

library(tidyverse)
library(gridExtra)

## ---- 1. Read depth-profile data -------------------------------------------

Depth_profile_converted <- read.csv("../data/PW_OM_nutrients_depth_profile.csv", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::rename(Season = season, Date = date, Phase = phase, Treatment = treatment, Depth = depth,
                Parameter = parameter, Phase_alpha = phase_alpha) %>%
  mutate(
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    Parameter = factor(Parameter, levels = c("OM_LOI", "NH4_N", "NO3.NO2_N", "PO4_P", "Si")),
    Treatment = factor(Treatment, levels = c("Control", "MHW")),
    value_mmol = if_else(Parameter == "OM_LOI", value, value / 1000)
  )

## ---- 2. Build one depth-profile panel per parameter -----------------------

make_depth_plot <- function(df, param_name, x_limits, x_breaks, x_label = "Value") {
  ggplot(
    data = subset(df, Parameter == param_name),
    aes(x = value_mmol, y = Depth, color = Treatment, shape = Phase)
  ) +
    geom_path(aes(alpha = Phase_alpha, group = interaction(Treatment, Phase)), linewidth = 0.75) +
    geom_point(aes(alpha = Phase_alpha), size = 2) +
    scale_color_manual(values = c("Control" = "steelblue1", "MHW" = "#D95F02")) +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    scale_y_reverse(breaks = c(0, 2.5, 5, 7.5, 10)) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    facet_grid(. ~ Season) +
    labs(x = x_label, y = "Depth [cm]", title = param_name) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.ticks.length = unit(3, "pt")
    )
}

p_om  <- make_depth_plot(Depth_profile_converted, "OM_LOI",     c(9, 18.5), c(10, 13, 16), "OM [% LOI]")
p_nh4 <- make_depth_plot(Depth_profile_converted, "NH4_N",      c(0, 1.05), c(0.0, 0.4, 0.8), expression(NH[4]~"(mmol L"^{-1}*")"))
p_nox <- make_depth_plot(Depth_profile_converted, "NO3.NO2_N",  c(0, 0.02), c(0, 0.007, 0.014), expression(NO[x]~"(mmol L"^{-1}*")"))
p_po4 <- make_depth_plot(Depth_profile_converted, "PO4_P",      c(0, 0.3),  c(0.0, 0.1, 0.2), expression(PO[4]~"(mmol L"^{-1}*")"))
p_si  <- make_depth_plot(Depth_profile_converted, "Si",         c(0, 0.7),  c(0.0, 0.3, 0.6), expression(Si~"(mmol L"^{-1}*")"))

## ---- 3. Stack into Figure 3 (legend kept only on the bottom panel) -------

p_om  <- p_om  + theme(legend.position = "none")
p_nh4 <- p_nh4 + theme(legend.position = "none")
p_nox <- p_nox + theme(legend.position = "none")
p_po4 <- p_po4 + theme(legend.position = "none")
# p_si keeps its legend

Figure3 <- grid.arrange(p_om, p_nh4, p_nox, p_po4, p_si, ncol = 1)

ggsave("../outputs/figures/Figure3.tiff", plot = Figure3,
       width = 6, height = 12, dpi = 300, compression = "lzw")
