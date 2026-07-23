## ---------------------------------------------------------------------------
## Figure 4: Seasonal oxygen and nutrient flux rates (Control vs. MHW)
##
## Computes incubation flux rates (Oxygen, NH4, NOx, PO4, Si) from raw
## Start/End incubation nutrient concentrations, converts to mmol m^-2 d^-1,
## and plots the mean +/- SD rate per sampling occasion, faceted by
## parameter (rows) and season (columns).
##
## Input:
##   - Seasonal_MHW_nutrients.csv
##
## Output:
##   - Figure4.tiff
## ---------------------------------------------------------------------------

library(tidyverse)
library(lubridate)

## ---- 1. Read raw incubation nutrient/oxygen data ---------------------------

Nutrients <- read.csv("../data/Seasonal_MHW_nutrients.csv", sep = ",", dec = ".",
                       header = TRUE, stringsAsFactors = TRUE) %>%
  dplyr::rename(Season = season, Date = date, DateTime = datetime, Sampling = sampling, Phase = phase,
                Treatment = treatment, Replicate = replicate, Timepoint = timepoint,
                Oxygen_mg.L = oxygen, NH4_N = nh4_n, NO3.NO2_N = no3_no2_n, PO4_P = po4_p, Si = si,
                Water_volume.L. = water_volume)

filtered_Nutrients <- Nutrients %>%
  mutate(
    Date = as.Date(as.character(Date)),
    DateTime = dmy_hm(as.character(DateTime))
  )

## ---- 2. Pivot Start/End timepoints wide and compute incubation duration --

calculate_duration <- function(df) {
  df %>%
    filter(Timepoint %in% c("Start", "End")) %>%
    pivot_wider(
      id_cols = c(Season, Date, Sampling, Phase, Treatment, Replicate, Water_volume.L.),
      names_from = Timepoint,
      values_from = c(DateTime, Oxygen_mg.L, NH4_N, NO3.NO2_N, PO4_P, Si),
      names_glue = "{Timepoint}_{.value}"
    ) %>%
    mutate(across(ends_with("DateTime"), as.POSIXct, origin = "1970-01-01")) %>%
    mutate(Duration = as.numeric(difftime(End_DateTime, Start_DateTime, units = "hours")))
}

duration_data <- calculate_duration(filtered_Nutrients)

## ---- 3. Compute flux rate per parameter ------------------------------------
## Rate = (End - Start) / Duration * Water_volume * 24h * area-conversion factor.
## Conversion factor scales from the experimental incubation footprint (GEMAX
## cores, ~9 cm diameter) to per-m^2 rates.

calculate_rate <- function(df, variable_name, conversion_factor) {
  start_col <- paste0("Start_", variable_name)
  end_col <- paste0("End_", variable_name)

  df %>%
    mutate(
      Rate = ((get(end_col) - get(start_col)) / Duration) * Water_volume.L. * 24 * conversion_factor,
      Parameter = variable_name
    ) %>%
    select(Season, Sampling, Phase, Treatment, Replicate, Rate, Parameter)
}

Rates <- bind_rows(
  calculate_rate(duration_data, "NH4_N", 157.2),
  calculate_rate(duration_data, "Oxygen_mg.L", 157.2),
  calculate_rate(duration_data, "NO3.NO2_N", 157.2),
  calculate_rate(duration_data, "PO4_P", 157.2),
  calculate_rate(duration_data, "Si", 157.2)
)

## ---- 4. Convert to mmol m^-2 d^-1 ------------------------------------------

Rates_converted <- Rates %>%
  mutate(Rate_mmol = case_when(
    Parameter == "Oxygen_mg.L" ~ Rate / 31.999,  # mg -> mmol (O2 molar mass)
    TRUE ~ Rate / 1000                            # ug -> mmol for NH4/PO4/Si/NOx
  )) %>%
  mutate(Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn")))

## ---- 5. Mean +/- SD per Season x Phase x Treatment x Parameter -----------

Rates_converted_av <- Rates_converted %>%
  group_by(Season, Phase, Treatment, Parameter) %>%
  summarise(
    Rate_av = mean(Rate_mmol, na.rm = TRUE),
    Rate_sd = sd(Rate_mmol, na.rm = TRUE),
    Rate_n = sum(!is.na(Rate_mmol)),
    Rate_se = Rate_sd / sqrt(Rate_n),
    .groups = "drop"
  ) %>%
  mutate(Parameter = ordered(Parameter, levels = c("Oxygen_mg.L", "NH4_N", "PO4_P", "NO3.NO2_N", "Si")))

## ---- 6. Figure 4 -----------------------------------------------------------

pd <- position_dodge(0.75)

Figure4 <- ggplot(Rates_converted_av, aes(x = Phase, y = Rate_av, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = Rate_av - Rate_sd, ymax = Rate_av + Rate_sd, color = Treatment),
                width = 0.75, linewidth = 1.5, position = pd) +
  scale_color_manual(values = c("Control" = "steelblue1", "MHW" = "#D95F02")) +
  facet_grid(Parameter ~ Season, scales = "free_y") +
  xlab("Sampling") +
  ylab(bquote(Rate ~ "[" * "mmol " ~ m^-2 ~ d^-1 * "]")) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1.25, "cm"),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 26, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

Figure4

ggsave("../outputs/figures/Figure4.tiff", plot = Figure4,
       width = 8, height = 12, dpi = 300, compression = "lzw")
