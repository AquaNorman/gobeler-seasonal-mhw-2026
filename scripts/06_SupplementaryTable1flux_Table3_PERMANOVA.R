## ---------------------------------------------------------------------------
## Supplementary Table 1 (biogeochemical fluxes part) & Supplementary Table 3
## PERMANOVA / PERMDISP on oxygen and nutrient flux rates
##
## Supplementary Table 1: overall PERMANOVA (Treatment * Season), each of
##   Oxygen and combined Nutrients (NH4, NOx, PO4, Si) tested separately.
## Supplementary Table 3: per-season PERMANOVA (Treatment * Sampling), same
##   two response groups. "Phase" in the code corresponds to "Sampling" in
##   the manuscript.
##
## Input:
##   - Fluxrates.csv
## ---------------------------------------------------------------------------

library(tidyverse)
library(vegan)

## ---- 1. Read and reshape flux rate data -----------------------------------

Fluxrates_raw <- read.csv("../data/Fluxrates.csv", sep = ",", dec = ".", header = TRUE,
                           stringsAsFactors = TRUE) %>%
  dplyr::rename(Season = season, Phase = phase, Treatment = treatment, Replicate = replicate,
                Rate = rate_raw, Parameter = parameter, Rate_mmol = rate_mmol)

Fluxrates <- Fluxrates_raw %>%
  select(Season, Phase, Treatment, Replicate, Parameter, Rate_mmol) %>%
  filter(Phase != "Acclimation") %>%
  pivot_wider(names_from = Parameter, values_from = Rate_mmol) %>%
  dplyr::rename(O2_rate_mmol = Oxygen_mg.L, NH4_rate_mmol = NH4_N, NO_rate_mmol = NO3.NO2_N,
                Si_rate_mmol = Si, PO4_rate_mmol = PO4_P)

## ---- 2. Supplementary Table 1: overall PERMANOVA (Treatment * Season) ----

norm_data <- Fluxrates %>%
  select(Season, Phase, Treatment, Replicate) %>%
  mutate(
    O2_rate_mmol = as.numeric(scale(Fluxrates$O2_rate_mmol)),
    NH4_rate_mmol = as.numeric(scale(Fluxrates$NH4_rate_mmol)),
    NO_rate_mmol = as.numeric(scale(Fluxrates$NO_rate_mmol)),
    Si_rate_mmol = as.numeric(scale(Fluxrates$Si_rate_mmol)),
    PO4_rate_mmol = as.numeric(scale(Fluxrates$PO4_rate_mmol)),
    ChamberID = paste(Treatment, Replicate, sep = "_")
  )

dist_oxy  <- vegdist(norm_data[, "O2_rate_mmol", drop = FALSE], method = "euclidean")
dist_nuts <- vegdist(norm_data[, c("NH4_rate_mmol", "PO4_rate_mmol", "Si_rate_mmol", "NO_rate_mmol")], method = "euclidean")

SuppTable1_permanova_oxy  <- adonis2(dist_oxy  ~ Treatment * Season, data = norm_data, strata = norm_data$ChamberID, by = "terms")
SuppTable1_permanova_nuts <- adonis2(dist_nuts ~ Treatment * Season, data = norm_data, strata = norm_data$ChamberID, by = "terms")

print(SuppTable1_permanova_oxy)
print(SuppTable1_permanova_nuts)

## ---- 3. Supplementary Table 1: PERMDISP -----------------------------------

norm_data$TreatmentSeason <- paste(norm_data$Treatment, norm_data$Season, sep = "_")

run_permdisp <- function(dist_obj, grouping, label) {
  pd <- betadisper(dist_obj, grouping)
  pt <- permutest(pd, pairwise = TRUE)
  cat("\n---", label, "---\n")
  print(pt)
  pt
}

SuppTable1_permdisp_oxy <- list(
  Treatment       = run_permdisp(dist_oxy, norm_data$Treatment, "Oxygen ~ Treatment"),
  Season          = run_permdisp(dist_oxy, norm_data$Season, "Oxygen ~ Season"),
  TreatmentSeason = run_permdisp(dist_oxy, norm_data$TreatmentSeason, "Oxygen ~ Treatment:Season")
)

SuppTable1_permdisp_nuts <- list(
  Treatment       = run_permdisp(dist_nuts, norm_data$Treatment, "Nutrients ~ Treatment"),
  Season          = run_permdisp(dist_nuts, norm_data$Season, "Nutrients ~ Season"),
  TreatmentSeason = run_permdisp(dist_nuts, norm_data$TreatmentSeason, "Nutrients ~ Treatment:Season")
)

## ---- 4. Supplementary Table 3: per-season PERMANOVA (Treatment * Sampling) ----

normalize_by_season <- function(data) {
  data %>%
    group_by(Season) %>%
    mutate(
      O2_rate_mmol = as.numeric(scale(O2_rate_mmol)),
      NH4_rate_mmol = as.numeric(scale(NH4_rate_mmol)),
      Si_rate_mmol = as.numeric(scale(Si_rate_mmol)),
      PO4_rate_mmol = as.numeric(scale(PO4_rate_mmol)),
      NO_rate_mmol = as.numeric(scale(NO_rate_mmol))
    ) %>%
    ungroup()
}

norm_data_season <- Fluxrates %>%
  select(Season, Phase, Treatment, Replicate, O2_rate_mmol, NH4_rate_mmol, Si_rate_mmol, PO4_rate_mmol, NO_rate_mmol) %>%
  normalize_by_season() %>%
  mutate(TreatmentPhase = paste0(Treatment, Phase))

run_season_flux_permanova <- function(season_name) {
  df <- norm_data_season %>% filter(Season == season_name)

  d_oxy  <- vegdist(df[, "O2_rate_mmol", drop = FALSE], method = "euclidean")
  d_nuts <- vegdist(df[, c("NH4_rate_mmol", "PO4_rate_mmol", "Si_rate_mmol", "NO_rate_mmol")], method = "euclidean")

  cat("\n===", season_name, "- Oxygen ===\n")
  pn_oxy <- adonis2(d_oxy ~ Treatment * Phase, data = df, strata = df$Replicate, by = "terms")
  print(pn_oxy)

  cat("\n===", season_name, "- Nutrients ===\n")
  pn_nuts <- adonis2(d_nuts ~ Treatment * Phase, data = df, strata = df$Replicate, by = "terms")
  print(pn_nuts)

  df$TreatmentPhase <- paste0(df$Treatment, df$Phase)
  pd_oxy_treatment       <- run_permdisp(d_oxy, df$Treatment, paste(season_name, "- Oxygen ~ Treatment"))
  pd_oxy_phase           <- run_permdisp(d_oxy, df$Phase, paste(season_name, "- Oxygen ~ Sampling"))
  pd_oxy_treatmentphase  <- run_permdisp(d_oxy, df$TreatmentPhase, paste(season_name, "- Oxygen ~ Treatment:Sampling"))
  pd_nuts_treatment      <- run_permdisp(d_nuts, df$Treatment, paste(season_name, "- Nutrients ~ Treatment"))
  pd_nuts_phase          <- run_permdisp(d_nuts, df$Phase, paste(season_name, "- Nutrients ~ Sampling"))
  pd_nuts_treatmentphase <- run_permdisp(d_nuts, df$TreatmentPhase, paste(season_name, "- Nutrients ~ Treatment:Sampling"))

  list(permanova_oxy = pn_oxy, permanova_nuts = pn_nuts,
       permdisp_oxy_treatment = pd_oxy_treatment, permdisp_oxy_phase = pd_oxy_phase,
       permdisp_oxy_treatmentphase = pd_oxy_treatmentphase,
       permdisp_nuts_treatment = pd_nuts_treatment, permdisp_nuts_phase = pd_nuts_phase,
       permdisp_nuts_treatmentphase = pd_nuts_treatmentphase)
}

SuppTable3 <- map(c("Winter", "Spring", "Summer", "Autumn"), run_season_flux_permanova) %>%
  set_names(c("Winter", "Spring", "Summer", "Autumn"))
