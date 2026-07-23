## ---------------------------------------------------------------------------
## Figure 1: Seasonal Marine Heatwaves (Reference Period 1931-2020)
##
## Input:
##   - Daily_average_Temperature_MHW.csv       experimental heatwave treatment
##   - Daily_average_Temperature_controls.csv  experimental control treatment
##   - YSI_Storfjarden_MHW.csv                 MONICOAST field monitoring data
##
## Output:
##   - Figure1.tiff
## ---------------------------------------------------------------------------

library(tidyverse)
library(heatwaveR)
library(ggpubr)

## ---- 1. MHW experimental temperature: detect heatwave events -------------

mhw_daily <- read.csv("../data/Daily_average_Temperature_MHW.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename(Date = date, Treatment = treatment, Season = season, Phase = phase,
                MeanTemp = mean_temp, RP31.20_mean = rp31_20_mean, RP31.20_Q10 = rp31_20_q10,
                RP31.20_Q90 = rp31_20_q90) %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  select(Date, doy, Phase, Treatment, MeanTemp, RP31.20_mean, RP31.20_Q90) %>%
  arrange(Date) %>%
  rename(t = Date, temp = MeanTemp, seas = RP31.20_mean, thresh = RP31.20_Q90)

# 2x and 3x thresholds, as defined in heatwaveR category conventions
mhw_daily$thresh_2x <- mhw_daily$thresh + (mhw_daily$thresh - mhw_daily$seas)
mhw_daily$thresh_3x <- mhw_daily$thresh_2x + (mhw_daily$thresh - mhw_daily$seas)

mhw <- detect_event(mhw_daily)

## ---- 2. Rename climatology output to figure-ready names ------------------

mhw2 <- mhw$climatology %>%
  rename(Date = t, Temperature = temp, Climatology = seas,
         Threshold = thresh, Threshold_2x = thresh_2x, Threshold_3x = thresh_3x)

## ---- 3. Control incubation temperature ------------------------------------

Daily_Control_Days <- read.csv("../data/Daily_average_Temperature_controls.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename(Date = date, Treatment = treatment, MeanTemp = mean_temp) %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

## ---- 4. Field monitoring data (MONICOAST / Storfjärden) ------------------

Stor_MHW_raw <- read.csv("../data/YSI_Storfjarden_MHW.csv", stringsAsFactors = FALSE) %>%
  dplyr::rename(Date = date, Temperature = temperature, Climatological_Mean = climatological_mean,
                Threshold = threshold, Threshold_2x = threshold_2x, Threshold_3x = threshold_3x) %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  rename(Climatology = Climatological_Mean)

## ---- 5. Build Figure 1 -----------------------------------------------------

mhw2_plot <- mhw2 %>% mutate(Source = "Seasonal Experimental data")
raw_plot  <- Stor_MHW_raw %>% mutate(Source = "MONICOAST/Storfjärden field data")

Figure1 <- ggplot() +
  # each data source gets its own geom_flame() call - combining sources into
  # one data frame before plotting produces incorrect ribbons between sources
  geom_flame(data = mhw2_plot, aes(x = Date, y = Temperature, y2 = Threshold, fill = "Moderate"), alpha = 0.7) +
  geom_flame(data = mhw2_plot, aes(x = Date, y = Temperature, y2 = Threshold_2x, fill = "Strong"), alpha = 0.7) +
  geom_flame(data = raw_plot, aes(x = Date, y = Temperature, y2 = Threshold, fill = "Moderate"), alpha = 0.3) +
  geom_flame(data = raw_plot, aes(x = Date, y = Temperature, y2 = Threshold_2x, fill = "Strong"), alpha = 0.3) +
  geom_line(data = raw_plot, aes(x = Date, y = Temperature, colour = "Long-term observations"), linewidth = 1, alpha = 0.3) +
  geom_line(data = mhw2_plot, aes(x = Date, y = Temperature, colour = "Experimental period"), linewidth = 1) +
  geom_line(data = Daily_Control_Days, aes(x = Date, y = MeanTemp, colour = "Control incubation temperature"), linewidth = 1.2) +
  geom_line(data = mhw2, aes(x = Date, y = Climatology, colour = "Climatological mean"), linewidth = 0.8) +
  geom_line(data = mhw2, aes(x = Date, y = Threshold, linetype = "MHW threshold"), colour = "grey50", linewidth = 0.8) +
  geom_line(data = mhw2, aes(x = Date, y = Threshold_2x, linetype = "Strong MHW threshold"), colour = "grey50", linewidth = 0.8) +
  scale_colour_manual(
    name = "Temperature series",
    values = c(
      "Long-term observations" = "grey50",
      "Experimental period" = "black",
      "Control incubation temperature" = "blue",
      "Climatological mean" = "steelblue1"
    )
  ) +
  scale_linetype_manual(
    name = "Thresholds",
    values = c("MHW threshold" = "dashed", "Strong MHW threshold" = "dotted")
  ) +
  scale_fill_manual(
    name = "MHW intensity",
    values = c("Moderate" = "#ffc866", "Strong" = "#ff6900")
  ) +
  labs(y = expression(paste("Temperature [", degree, "C]")), x = "Date",
       title = "Seasonal Marine Heatwaves RP1931-2020") +
  font("xlab", size = 26) + font("ylab", size = 26) + font("xy.text", size = 22) +
  font("legend.text", size = 24) + font("legend.title", size = 24) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2))

Figure1

## ---- 6. Save ---------------------------------------------------------------

ggsave("../outputs/figures/Figure1.tiff", plot = Figure1, width = 16, height = 7,
       dpi = 300, compression = "lzw")
