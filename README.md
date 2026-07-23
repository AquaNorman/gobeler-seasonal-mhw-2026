# Season-specific effects of Marine Heatwaves on macrofaunal community composition and ecosystem processes

Code accompanying:

> Göbeler N, Norkko A, Norkko J (2026) Season-specific effects of Marine Heatwaves
> on macrofaunal community composition and ecosystem processes. *Ecosystems*.
> [DOI to be added on publication]

## Summary

Marine heatwaves increasingly affect coastal ecosystems, but most research has
focused on summer extremes. This study exposed intact sediment cores from
Storfjärden, Baltic Sea, to control and marine-heatwave temperatures across all
four seasons of 2021, measuring macrofaunal community structure, vertical
distribution, sediment organic matter, porewater nutrients, and oxygen/nutrient
flux rates. This repository contains the R scripts used to process the data and
reproduce all figures and statistical results reported in the manuscript.

## Data availability

The data underlying these scripts are openly available in the Bolin Centre
Database: [DOI to be added].

## Repository structure

```
scripts/    R scripts, numbered in the order they should be run
data/       not included here - download from the Bolin Centre Database (see above)
outputs/    figures/ and tables/ are created here when the scripts are run
```

Each script expects the `data/` folder (downloaded separately from Bolin) to sit
alongside `scripts/` and `outputs/` as shown above, and reads/writes using
relative paths. Run scripts from within the `scripts/` folder, or open the
repository as an RStudio Project.

## Scripts

| Script | Produces |
|---|---|
| `01_Figure1_seasonal_MHW.R` | Figure 1 |
| `02_Figure2_Table1_community.R` | Figure 2, Table 1, Supplementary Figures 1-2 |
| `03_Figure3_depth_profiles.R` | Figure 3 |
| `04_Figure4_flux_rates.R` | Figure 4 |
| `05_SupplementaryTables1-2_community_PERMANOVA.R` | Supplementary Table 1 (community part), Supplementary Table 2 |
| `06_SupplementaryTable1flux_Table3_PERMANOVA.R` | Supplementary Table 1 (flux part), Supplementary Table 3 |
| `07_SupplementaryFigure3_flux_distances.R` | Supplementary Figure 3 |

Script `05` sources `02` to reuse the community abundance/biomass data; all
other scripts are self-contained given the data files they require (see the
comment header of each script for its specific inputs).

## Requirements

R (tested under R 4.x) with the following packages:

`tidyverse`, `heatwaveR`, `ggpubr`, `gridExtra`, `vegan`, `gtools`, `lubridate`

## Reproducibility notes

* Supplementary Figure 3 involves a bootstrap step for Season x Sampling
  combinations with unequal replicate numbers (see script header for details).
  A seed is set for reproducibility within a run, but exact values may differ
  slightly between R versions/platforms.
* PERMANOVA/PERMDISP p-values (Supplementary Tables 1-3) are derived from
  random permutation and may vary slightly (typically ±0.01-0.02) between runs;
  test statistics (R², pseudo-F) are deterministic given the data.

## License

Code is released under the MIT License (see `LICENSE`).

## Contact

Norman Göbeler, Tvärminne Zoological Station, University of Helsinki
norman.gobeler@helsinki.fi
