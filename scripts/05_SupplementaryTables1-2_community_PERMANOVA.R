## ---------------------------------------------------------------------------
## Supplementary Table 1 (community part) & Supplementary Table 2
## PERMANOVA / PERMDISP on macrofaunal community abundance and biomass
##
## Supplementary Table 1: overall PERMANOVA (Sediment depth * Season * Treatment).
##   This script covers the community half only; the biogeochemical-flux half
##   of Supplementary Table 1, and Supplementary Table 3 / Supplementary
##   Figure 3, are in 06_SupplementaryTable1flux_Table3_PERMANOVA.R and
##   07_SupplementaryFigure3_flux_distances.R.
## Supplementary Table 2: per-season PERMANOVA (Sediment depth * Treatment).
##
## Sources 02_Figure2_Table1_community.R for the abundance/biomass data. Run
## this script from the scripts/ directory (or open the repository as an
## RStudio project) so the relative path below resolves correctly.
## ---------------------------------------------------------------------------

library(tidyverse)
library(vegan)

source("02_Figure2_Table1_community.R")

## ---- 1. Restrict to species present in >33% of cores ----------------------

Abundance_summarized <- Abundance %>%
  group_by(Season, Treatment, Replicate) %>%
  summarise(across(c(Macoma, Hydrobidae, Ostracoda, Oligochaete, Marenzelleria,
                      Monoporeira, Chironimadae, Theodoxus, Sinelobus, Harmothoe),
                    ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

number_of_cores <- nrow(Abundance_summarized)
species_presence <- colSums(Abundance_summarized[, -(1:3)] > 0)
selected_species <- names(species_presence[species_presence > (number_of_cores * 0.33)])
selected_species
# Expected: Macoma, Hydrobidae, Ostracoda, Marenzelleria, Monoporeira, Chironimadae
# (Theodoxus, Sinelobus, Harmothoe/Bylgides sarsi, Oligochaeta excluded - too rare)

## ---- 2. Wide-format abundance & biomass for community-level PERMANOVA ----

Abundance_wide_prep <- Abundance_long %>%
  filter(!Species %in% c("Theodoxus", "Sinelobus", "Harmothoe", "Oligochaete"))

Abundance_wide <- pivot_wider(Abundance_wide_prep, names_from = Species, values_from = Abundance)
Abundance_wide[is.na(Abundance_wide)] <- 0

Abundance_sub <- Abundance_wide[, c("Macoma_ad", "Macoma_juv", "Hydrobidae", "Ostracoda",
                                     "Marenzelleria", "Monoporeira", "Chironimadae")]
Abundance_dist <- vegdist(Abundance_sub, method = "bray")

Biomass_wide <- Biomass %>%
  select(-c(Theodoxus, Sinelobus, Harmothoe, Oligochaete))
Biomass_wide$Dummy <- 0.001
Biomass_wide[is.na(Biomass_wide)] <- 0

Biomass_sub <- Biomass_wide[, c("Macoma", "Hydrobidae", "Ostracoda", "Marenzelleria",
                                 "Monoporeira", "Chironimadae", "Dummy")]
Biomass_dist <- vegdist(Biomass_sub, method = "bray")

## ---- 3. Supplementary Table 1 (community part): overall PERMANOVA --------

permanova_Abundance <- adonis2(Abundance_dist ~ Sediment_depth * Season * Treatment,
                                data = Abundance_wide, permutations = 999, by = "terms")
print(permanova_Abundance)

permanova_Biomass <- adonis2(Biomass_dist ~ Sediment_depth * Season * Treatment,
                              data = Biomass_wide, permutations = 999, by = "terms")
print(permanova_Biomass)

## ---- 4. Supplementary Table 1 (community part): PERMDISP -----------------

Abundance_wide <- Abundance_wide %>%
  mutate(
    DepthSeason = as.factor(paste0(Sediment_depth, Season)),
    DepthTreatment = as.factor(paste0(Sediment_depth, Treatment)),
    SeasonTreatment = as.factor(paste0(Season, Treatment)),
    DepthSeasonTreatment = as.factor(paste0(Sediment_depth, Season, Treatment))
  )

Biomass_wide <- Biomass_wide %>%
  mutate(
    DepthSeason = as.factor(paste0(Sediment_depth, Season)),
    DepthTreatment = as.factor(paste0(Sediment_depth, Treatment)),
    SeasonTreatment = as.factor(paste0(Season, Treatment)),
    DepthSeasonTreatment = as.factor(paste0(Sediment_depth, Season, Treatment))
  )

run_permdisp <- function(dist_obj, grouping, label) {
  pd <- betadisper(dist_obj, grouping)
  pt <- permutest(pd, pairwise = TRUE)
  cat("\n---", label, "---\n")
  print(pt)
  pt
}

permdisp_results_Abundance <- list(
  Depth               = run_permdisp(Abundance_dist, Abundance_wide$Sediment_depth, "Abundance ~ Sediment_depth"),
  Season              = run_permdisp(Abundance_dist, Abundance_wide$Season, "Abundance ~ Season"),
  Treatment           = run_permdisp(Abundance_dist, Abundance_wide$Treatment, "Abundance ~ Treatment"),
  DepthSeason         = run_permdisp(Abundance_dist, Abundance_wide$DepthSeason, "Abundance ~ Depth:Season"),
  DepthTreatment      = run_permdisp(Abundance_dist, Abundance_wide$DepthTreatment, "Abundance ~ Depth:Treatment"),
  SeasonTreatment     = run_permdisp(Abundance_dist, Abundance_wide$SeasonTreatment, "Abundance ~ Season:Treatment"),
  DepthSeasonTreatment = run_permdisp(Abundance_dist, Abundance_wide$DepthSeasonTreatment, "Abundance ~ Depth:Season:Treatment")
)

permdisp_results_Biomass <- list(
  Depth               = run_permdisp(Biomass_dist, Biomass_wide$Sediment_depth, "Biomass ~ Sediment_depth"),
  Season              = run_permdisp(Biomass_dist, Biomass_wide$Season, "Biomass ~ Season"),
  Treatment           = run_permdisp(Biomass_dist, Biomass_wide$Treatment, "Biomass ~ Treatment"),
  DepthSeason         = run_permdisp(Biomass_dist, Biomass_wide$DepthSeason, "Biomass ~ Depth:Season"),
  DepthTreatment      = run_permdisp(Biomass_dist, Biomass_wide$DepthTreatment, "Biomass ~ Depth:Treatment"),
  SeasonTreatment     = run_permdisp(Biomass_dist, Biomass_wide$SeasonTreatment, "Biomass ~ Season:Treatment"),
  DepthSeasonTreatment = run_permdisp(Biomass_dist, Biomass_wide$DepthSeasonTreatment, "Biomass ~ Depth:Season:Treatment")
)

## ---- 5. Supplementary Table 2: per-season PERMANOVA + PERMDISP -----------

run_season_permanova <- function(wide_df, season_name, species_cols) {
  df <- wide_df %>% filter(Season == season_name)
  sub <- df[, species_cols]
  dist_obj <- vegdist(sub, method = "bray")
  cat("\n===", season_name, "===\n")
  pn <- adonis2(dist_obj ~ Sediment_depth * Treatment, data = df, permutations = 999, by = "terms")
  print(pn)
  pd <- betadisper(dist_obj, df$Sediment_depth)
  pt <- permutest(pd, pairwise = TRUE)
  print(pt)
  list(permanova = pn, permdisp = pt)
}

abundance_species_cols <- c("Macoma_ad", "Macoma_juv", "Hydrobidae", "Ostracoda",
                             "Marenzelleria", "Monoporeira", "Chironimadae")
biomass_species_cols <- c("Macoma", "Hydrobidae", "Ostracoda", "Marenzelleria",
                           "Monoporeira", "Chironimadae", "Dummy")

SuppTable2_Abundance <- map(
  c("Winter", "Spring", "Summer", "Autumn"),
  ~ run_season_permanova(Abundance_wide, .x, abundance_species_cols)
) %>% set_names(c("Winter", "Spring", "Summer", "Autumn"))

SuppTable2_Biomass <- map(
  c("Winter", "Spring", "Summer", "Autumn"),
  ~ run_season_permanova(Biomass_wide, .x, biomass_species_cols)
) %>% set_names(c("Winter", "Spring", "Summer", "Autumn"))
