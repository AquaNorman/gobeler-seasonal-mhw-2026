## ---------------------------------------------------------------------------
## Figure 2 & Table 1: Seasonal macrofaunal community abundance and biomass
##
## Reproduces:
##   - Table 1 (seasonal abundance & biomass of key taxa)
##   - Figure 2 (vertical distribution of biomass across sediment depth)
##   - Supplementary Figure 1 (abundance boxplot, all taxa)
##   - Supplementary Figure 2 (biomass boxplot, taxa with detectable wet weight)
##
## PERMANOVA/PERMDISP statistics (Supplementary Tables 1-2) are in a separate
## script that sources this one for the abundance/biomass data.
##
## Input:
##   - Community.csv          per-core abundance & biomass counts
##   - Community_sizes.csv    individual specimen sizes (for Macoma size split)
##
## Output:
##   - Figure2.tiff
##   - Table1_abundance.csv, Table1_biomass.csv
##   - SupplementaryFigure1.tiff, SupplementaryFigure2.tiff
## ---------------------------------------------------------------------------

library(tidyverse)
library(gridExtra)

## ---- 1. Read raw community data -------------------------------------------

Abundance_Biomass <- read.csv("../data/Community.csv", header = TRUE, sep = ",", dec = ".",
                               stringsAsFactors = TRUE, na.strings = c("#N/A", "NA")) %>%
  dplyr::rename(Date = date, Season = season, Parameter = parameter, Phase = phase,
                Treatment = treatment, Replicate = replicate, Sediment_depth = sediment_depth,
                Macoma = macoma, Hydrobidae = hydrobidae, Ostracoda = ostracoda,
                Marenzelleria = marenzelleria, Monoporeira = monoporeira, Chironimadae = chironimadae,
                Oligochaete = oligochaete, Theodoxus = theodoxus, Sinelobus = sinelobus, Harmothoe = harmothoe)

Abundance <- Abundance_Biomass %>%
  filter(Parameter == "Community_Abund") %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

Biomass <- Abundance_Biomass %>%
  filter(Parameter == "Community_Mass")

## ---- 2. Split Macoma abundance into adult/juvenile using individual sizes -

Sizes <- read.csv("../data/Community_sizes.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors = FALSE) %>%
  dplyr::rename(Sample_Date = sample_date, Season = season,
                Phase = phase, Treatment = treatment, ID = id, Sediment_depth = sediment_depth,
                Species = species, Size_av = size_av, Size_median = size_median, Amount = amount) %>%
  dplyr::rename_with(~ gsub("^size", "Size", .x), starts_with("size"))

filtered_Sizes <- Sizes %>%
  select(-c(Sample_Date, Size_av, Amount, Phase, Size_median))

Sizes_long <- reshape2::melt(filtered_Sizes,
                              id.vars = c("Season", "Treatment", "ID", "Species", "Sediment_depth"),
                              variable.name = "Parameter", value.name = "Size")
Sizes_long$Parameter <- "Size"

Sizes_long_adjusted <- Sizes_long %>%
  mutate(Species = case_when(
    Species == "Macoma" & Size > 1  ~ "Macoma_ad",
    Species == "Macoma" & Size <= 1 ~ "Macoma_juv",
    TRUE ~ Species
  ))

# individuals per core -> individuals per m^2 (core area correction factor: 157.2)
Size_to_Abundance <- Sizes_long_adjusted %>%
  group_by(Season, Treatment, ID, Species, Sediment_depth) %>%
  summarise(Abundance = sum(!is.na(Size)) * 157.2, .groups = "drop")

Size_to_Abundance_Macoma <- Size_to_Abundance %>%
  filter(Species %in% c("Macoma_ad", "Macoma_juv")) %>%
  dplyr::rename(Replicate = ID)

## ---- 3. Long-format abundance (Macoma replaced by size-split values) ------

Abundance_long <- reshape2::melt(Abundance,
                                  id.vars = c("Date", "Parameter", "Season", "Treatment", "Replicate", "Phase", "Sediment_depth"),
                                  variable.name = "Species", value.name = "Abundance") %>%
  mutate(Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  select(-c(Date, Parameter, Phase)) %>%
  filter(Species != "Macoma") %>%
  bind_rows(Size_to_Abundance_Macoma) %>%
  mutate(Sediment_depth = as.factor(Sediment_depth))

## ---- 4. Long-format biomass -----------------------------------------------

Biomass_long <- reshape2::melt(Biomass,
                                id.vars = c("Date", "Parameter", "Season", "Treatment", "Replicate", "Phase", "Sediment_depth"),
                                variable.name = "Species", value.name = "Biomass") %>%
  mutate(Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn")))

## ---- 5. Published species names -------------------------------------------
## Maps data-entry codes to the taxon names used in the manuscript, including
## two data-entry corrections: "Chironimadae" -> Chironomidae and
## "Harmothoe" -> Bylgides sarsi.

abundance_species_names <- c(
  "Macoma_ad"     = "Macoma balthica (large)",
  "Macoma_juv"    = "Macoma balthica (juvenile)",
  "Hydrobidae"    = "Hydrobiidae",
  "Ostracoda"     = "Ostracoda",
  "Oligochaete"   = "Oligochaeta",
  "Marenzelleria" = "Marenzelleria spp.",
  "Monoporeira"   = "Monoporeia affinis",
  "Chironimadae"  = "Chironomidae",
  "Theodoxus"     = "Theodoxus fluviatilis",
  "Sinelobus"     = "Sinelobus sp.",
  "Harmothoe"     = "Bylgides sarsi"
)

biomass_species_names <- c(
  "Macoma"        = "Macoma balthica",
  "Hydrobidae"    = "Hydrobiidae",
  "Ostracoda"     = "Ostracoda",
  "Oligochaete"   = "Oligochaeta",
  "Marenzelleria" = "Marenzelleria spp.",
  "Monoporeira"   = "Monoporeia affinis",
  "Chironimadae"  = "Chironomidae",
  "Theodoxus"     = "Theodoxus fluviatilis",
  "Sinelobus"     = "Sinelobus sp.",
  "Harmothoe"     = "Bylgides sarsi"
)

## ---- 6. Table 1: seasonal abundance ---------------------------------------

Abundance_long_sum <- Abundance_long %>%
  group_by(Season, Treatment, Replicate, Species) %>%
  summarise(Sum_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

Abundance_av_per_season <- Abundance_long_sum %>%
  group_by(Season, Species) %>%
  summarise(Mean_Abundance = mean(Sum_abundance, na.rm = TRUE),
            sd = sd(Sum_abundance, na.rm = TRUE),
            n = sum(!is.na(Sum_abundance)),
            se = sd / sqrt(n),
            .groups = "drop") %>%
  mutate(Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
         Species = recode(as.character(Species), !!!abundance_species_names))

## ---- 7. Table 1: seasonal biomass (taxa with detectable wet weight) -------

Biomass_long_na <- Biomass_long %>% drop_na()

Biomass_long_sum <- Biomass_long_na %>%
  group_by(Season, Treatment, Replicate, Species) %>%
  summarise(Sum_Biomass = sum(Biomass, na.rm = TRUE), .groups = "drop")

Biomass_av_per_season <- Biomass_long_sum %>%
  group_by(Season, Species) %>%
  summarise(Mean_Biomass = mean(Sum_Biomass, na.rm = TRUE),
            sd = sd(Sum_Biomass, na.rm = TRUE),
            n = sum(!is.na(Sum_Biomass)),
            se = sd / sqrt(n),
            .groups = "drop") %>%
  filter(Species %in% c("Macoma", "Marenzelleria", "Monoporeira", "Chironimadae")) %>%
  mutate(Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
         Species = recode(as.character(Species), !!!biomass_species_names))

Abundance_av_per_season
Biomass_av_per_season

write.csv(Abundance_av_per_season, "../outputs/Table1_abundance.csv", row.names = FALSE)
write.csv(Biomass_av_per_season, "../outputs/Table1_biomass.csv", row.names = FALSE)

## ---- 8. Supplementary Figure 1: abundance boxplot, all taxa --------------

abundance_species_order <- c(
  "Hydrobidae", "Ostracoda", "Marenzelleria", "Monoporeira", "Chironimadae",
  "Oligochaete", "Theodoxus", "Sinelobus", "Harmothoe", "Macoma_ad", "Macoma_juv"
)

Abundance_long_named <- Abundance_long %>%
  mutate(Species = factor(as.character(Species), levels = abundance_species_order),
         Species = recode(Species, !!!abundance_species_names),
         Common_species = if_else(startsWith(as.character(Species), "Macoma balthica"), "Macoma", "Others"))

SupplementaryFigure1 <- ggplot(Abundance_long_named, aes(x = Treatment, y = Abundance)) +
  geom_boxplot(aes(fill = Species), linewidth = 0.25, varwidth = 4) +
  facet_grid(Common_species ~ Season, scales = "free") +
  ylab(bquote(Abundance ~ (Ind ~ m^-2))) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("Seasonal Abundance") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "right",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "transparent", colour = "transparent")
  )

SupplementaryFigure1

ggsave("../outputs/figures/SupplementaryFigure1.tiff", plot = SupplementaryFigure1,
       width = 12, height = 8, dpi = 300, compression = "lzw")

## ---- 9. Supplementary Figure 2: biomass boxplot, selected taxa -----------

biomass_species_order <- c("Macoma", "Marenzelleria", "Monoporeira", "Chironimadae")

Biomass_long_sum_named <- Biomass_long_sum %>%
  filter(Species %in% biomass_species_order) %>%
  mutate(Species = factor(as.character(Species), levels = biomass_species_order),
         Species = recode(Species, !!!biomass_species_names),
         Common_species = if_else(as.character(Species) == "Macoma balthica", "Macoma", "Others"))

SupplementaryFigure2 <- ggplot(Biomass_long_sum_named, aes(x = Treatment, y = Sum_Biomass)) +
  geom_boxplot(aes(fill = Species), linewidth = 0.25, varwidth = 4) +
  facet_grid(Common_species ~ Season, scales = "free") +
  ylab(bquote(Biomass ~ (g ~ m^-2))) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("Seasonal Biomass") +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "right",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 24),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "transparent", colour = "transparent")
  )

SupplementaryFigure2

ggsave("../outputs/figures/SupplementaryFigure2.tiff", plot = SupplementaryFigure2,
       width = 12, height = 8, dpi = 300, compression = "lzw")

## ---- 10. Figure 2: vertical distribution of biomass across sediment depth -

# Depth column: numeric midpoint (cm) for each sediment-depth slice
depth_lookup <- c("0" = 1.5, "3" = 6.5, "10" = 12)

Biomass_long_depth <- Biomass_long %>%
  mutate(Depth = depth_lookup[as.character(Sediment_depth)])

Depth_Biomass_av <- Biomass_long_depth %>%
  filter(Species %in% c("Macoma", "Marenzelleria", "Monoporeira", "Chironimadae")) %>%
  group_by(Season, Treatment, Species, Sediment_depth, Depth) %>%
  summarise(Mean_Biomass = mean(Biomass, na.rm = TRUE),
            sd = sd(Biomass, na.rm = TRUE),
            n = sum(!is.na(Biomass)),
            se = sd / sqrt(n),
            .groups = "drop")

Sum_Biomass <- Biomass_long_depth %>%
  group_by(Season, Treatment, Replicate, Sediment_depth, Depth) %>%
  summarise(community_biomass = if (all(is.na(Biomass))) NA_real_ else sum(Biomass, na.rm = TRUE),
            .groups = "drop")

Depth_Biomass_total <- Sum_Biomass %>%
  group_by(Season, Treatment, Sediment_depth, Depth) %>%
  summarise(Mean_Biomass = mean(community_biomass, na.rm = TRUE),
            sd = sd(community_biomass, na.rm = TRUE),
            n = sum(!is.na(community_biomass)),
            se = sd / sqrt(n),
            .groups = "drop") %>%
  mutate(Species = "Community Biomass")

Depth_Biomass_av <- bind_rows(Depth_Biomass_av, Depth_Biomass_total) %>%
  mutate(
    Season = ordered(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    # Figure 2 uses short display names (not the full Table 1 binomial names)
    Species = recode(Species, "Chironimadae" = "Chironomidae", "Monoporeira" = "Monoporeia"),
    Species = factor(Species, levels = c("Chironomidae", "Macoma", "Marenzelleria", "Monoporeia", "Community Biomass")),
    # offset Control/MHW points so they don't overlap at the same depth
    Depth = case_when(
      Treatment == "Control" ~ Depth + 0.5,
      Treatment == "MHW"     ~ Depth - 0.5,
      TRUE ~ Depth
    )
  )

pd <- position_dodge(2)

make_biomass_plot <- function(df, species_name, x_limits, x_breaks) {
  ggplot(data = subset(df, Species == species_name),
         aes(y = Depth, x = Mean_Biomass, color = Treatment)) +
    geom_point(position = pd, size = 3) +
    geom_errorbar(aes(xmin = Mean_Biomass - sd, xmax = Mean_Biomass + sd, color = Treatment),
                  width = 0.5, linewidth = 0.65) +
    scale_color_manual(values = c("Control" = "steelblue1", "MHW" = "#D95F02")) +
    facet_grid(~Season, scales = "free") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    labs(x = paste0(species_name, " [g/m", "\u00b2", "]"), y = "Sediment depth [cm]") +
    theme_minimal() +
    scale_y_reverse(breaks = c(0, 2.5, 5, 7.5, 10)) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      strip.text = element_text(size = 16),
      plot.title = element_text(size = 22, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.ticks.length = unit(3, "pt")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "gray")
}

p_chiro <- make_biomass_plot(Depth_Biomass_av, "Chironomidae",       c(-1, 7),    c(0, 3, 6))
p_mono  <- make_biomass_plot(Depth_Biomass_av, "Monoporeia",         c(-1, 7),    c(0, 3, 6))
p_maren <- make_biomass_plot(Depth_Biomass_av, "Marenzelleria",      c(-1, 7),    c(0, 3, 6))
p_maco  <- make_biomass_plot(Depth_Biomass_av, "Macoma",             c(-75, 350), c(0, 150, 300))
p_com   <- make_biomass_plot(Depth_Biomass_av, "Community Biomass",  c(-75, 350), c(0, 150, 300))

Figure2 <- grid.arrange(p_chiro, p_mono, p_maren, p_maco, p_com, ncol = 1)

ggsave("../outputs/figures/Figure2.tiff", plot = Figure2,
       width = 6, height = 12, dpi = 300, compression = "lzw")
