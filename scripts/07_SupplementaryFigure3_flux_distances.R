## ---------------------------------------------------------------------------
## Supplementary Figure 3: Euclidean dissimilarity (mean +/- SD) between
## Control and MHW of combined nutrient fluxes and oxygen fluxes, per season
## and sampling occasion (Acclimation, MHW_1, MHW_2, MHW_3).
##
## Method: per Season x Sampling combination, the distance between Control
## and MHW is the Euclidean distance between their (per-season standardized)
## rate vectors.
##   - Equal replicate counts: distance = minimum Euclidean distance found
##     across every possible ordering (permutation) of the MHW replicates
##     against the Control replicates. Deterministic.
##   - Unequal replicate counts: the larger group is randomly down-sampled to
##     match the smaller group 100 times; each draw's distance is itself a
##     permutation-minimum (as above), and the 100 draws are averaged
##     (mean +/- SD). Depends on the random seed set below.
##   - Oxygen: single response variable, so only the unequal-replicate
##     Season x Sampling combinations have an SD (from bootstrapping).
##   - Nutrients (NH4, NOx, PO4, Si): each permutation's distance is the mean
##     of the 4 individual nutrient distances; the SD shown is the spread
##     across those 4 nutrients at the best (minimum-mean) permutation, not a
##     bootstrap SD - so nutrient points all have error bars, unlike oxygen.
##
## Input:
##   - Fluxrates.csv
##
## Output:
##   - SupplementaryFigure3.tiff
## ---------------------------------------------------------------------------

library(tidyverse)
library(gtools)  # for permutations()

set.seed(123)
N_BOOTSTRAP <- 100

## ---- 1. Read flux data, keep Acclimation, standardize within season ------

Fluxrates_raw <- read.csv("../data/Fluxrates.csv", sep = ",", dec = ".", header = TRUE,
                           stringsAsFactors = TRUE) %>%
  dplyr::rename(Season = season, Phase = phase, Treatment = treatment, Replicate = replicate,
                Rate = rate_raw, Parameter = parameter, Rate_mmol = rate_mmol)

Fluxrates_acc <- Fluxrates_raw %>%
  select(Season, Phase, Treatment, Replicate, Parameter, Rate_mmol) %>%
  pivot_wider(names_from = Parameter, values_from = Rate_mmol) %>%
  dplyr::rename(O2_rate_mmol = Oxygen_mg.L, NH4_rate_mmol = NH4_N, NO_rate_mmol = NO3.NO2_N,
                Si_rate_mmol = Si, PO4_rate_mmol = PO4_P)

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

norm_data_season_acc <- Fluxrates_acc %>%
  select(Season, Phase, Treatment, Replicate, O2_rate_mmol, NH4_rate_mmol, Si_rate_mmol, PO4_rate_mmol, NO_rate_mmol) %>%
  normalize_by_season() %>%
  dplyr::rename(Sampling = Phase)  # published table/figure uses "Sampling", not "Phase"

## ---- 2. Distance functions -------------------------------------------------

# Euclidean distance between two equal-length rate vectors (one per parameter)
euclid_dist_vec <- function(ctrl_vec, mhw_vec) sqrt(sum((ctrl_vec - mhw_vec)^2))

# Distance for one permutation ordering, averaged across parameter columns
# (for oxygen there is only one column, so this is just that one distance)
permutation_distance <- function(ctrl_df, mhw_df_ordered) {
  col_dists <- map_dbl(seq_len(ncol(ctrl_df)), function(j) euclid_dist_vec(ctrl_df[[j]], mhw_df_ordered[[j]]))
  c(mean_distance = mean(col_dists), sd_across_params = sd(col_dists))
}

# Exhaustive permutation minimum: try every ordering of mhw_df's rows against
# ctrl_df, return the ordering with the smallest mean distance
min_permutation_distance <- function(ctrl_df, mhw_df) {
  n <- nrow(mhw_df)
  perms <- gtools::permutations(n = n, r = n)
  results <- t(apply(perms, 1, function(idx) permutation_distance(ctrl_df, mhw_df[idx, , drop = FALSE])))
  best <- which.min(results[, "mean_distance"])
  results[best, ]
}

# Exhaustive permutation-minimum if replicate counts match; otherwise
# bootstrap (100x down-sample to matching size, permutation-minimum each
# draw, then average)
compute_group_distance <- function(ctrl_df, mhw_df, n_bootstrap = N_BOOTSTRAP) {
  n_ctrl <- nrow(ctrl_df); n_mhw <- nrow(mhw_df)
  if (n_ctrl == n_mhw) {
    best <- min_permutation_distance(ctrl_df, mhw_df)
    tibble(Euclidean_Distance = best[["mean_distance"]],
           Distance_SD = best[["sd_across_params"]],
           Bootstrapped_Mean = NA_real_, Bootstrapped_SD = NA_real_)
  } else {
    n_min <- min(n_ctrl, n_mhw)
    boot_dists <- numeric(n_bootstrap)
    boot_param_sds <- numeric(n_bootstrap)
    for (b in seq_len(n_bootstrap)) {
      ctrl_sub <- ctrl_df[sample(seq_len(n_ctrl), n_min), , drop = FALSE]
      mhw_sub  <- mhw_df[sample(seq_len(n_mhw), n_min), , drop = FALSE]
      best <- min_permutation_distance(ctrl_sub, mhw_sub)
      boot_dists[b] <- best[["mean_distance"]]
      boot_param_sds[b] <- best[["sd_across_params"]]
    }
    tibble(Euclidean_Distance = mean(boot_dists),
           Distance_SD = mean(boot_param_sds),
           Bootstrapped_Mean = mean(boot_dists), Bootstrapped_SD = sd(boot_dists))
  }
}

## ---- 3. Compute distances for every Season x Sampling combination --------

compute_all <- function(data, param_cols, label) {
  combos <- data %>% distinct(Season, Sampling)
  results <- pmap_dfr(combos, function(Season, Sampling) {
    df <- data %>% filter(.data$Season == !!Season, .data$Sampling == !!Sampling)
    ctrl_df <- df %>% filter(Treatment == "Control") %>% select(all_of(param_cols))
    mhw_df  <- df %>% filter(Treatment == "MHW") %>% select(all_of(param_cols))
    bind_cols(Season = Season, Sampling = Sampling, compute_group_distance(ctrl_df, mhw_df))
  })
  results$Flux <- label
  results
}

oxygen_distances    <- compute_all(norm_data_season_acc, "O2_rate_mmol", "Oxygen flux")
nutrient_distances  <- compute_all(norm_data_season_acc, c("NH4_rate_mmol", "PO4_rate_mmol", "Si_rate_mmol", "NO_rate_mmol"), "Nutrient flux")

# Oxygen has only one response variable, so there's no cross-parameter spread
oxygen_distances$Distance_SD <- 0

all_distances <- bind_rows(oxygen_distances, nutrient_distances) %>%
  mutate(
    # combined error bar: bootstrap spread (if any) + cross-parameter spread (0 for oxygen)
    SD = if_else(!is.na(Bootstrapped_SD), Bootstrapped_SD + Distance_SD, Distance_SD),
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    Sampling = factor(Sampling, levels = c("Acclimation", "MHW_1", "MHW_2", "MHW_3"),
                       labels = c("Acclimation", "MHW 1", "MHW 2", "MHW 3")),
    Flux = factor(Flux, levels = c("Nutrient flux", "Oxygen flux"))
  )

## ---- 4. Supplementary Figure 3 --------------------------------------------

SupplementaryFigure3 <- ggplot(all_distances, aes(x = Sampling, y = Euclidean_Distance, color = Season, group = Season)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_line(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = Euclidean_Distance - SD, ymax = Euclidean_Distance + SD),
                width = 0.2, na.rm = TRUE, position = position_dodge(0.5)) +
  geom_vline(xintercept = seq_along(levels(all_distances$Sampling)) + 0.5,
             color = "gray", linetype = "dashed") +
  facet_grid(. ~ Flux) +
  labs(y = "Euclidean Distance", x = "Sampling") +
  scale_color_manual(values = c("Winter" = "lightblue", "Spring" = "lightgreen", "Summer" = "gold", "Autumn" = "tomato")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "transparent"),
    axis.text = element_text(size = 16),
    legend.position = "bottom"
  )

SupplementaryFigure3

ggsave("../outputs/figures/SupplementaryFigure3.tiff", plot = SupplementaryFigure3,
       width = 10, height = 6, dpi = 300, compression = "lzw")
