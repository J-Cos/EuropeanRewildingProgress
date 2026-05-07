# European Rewilding Progress

This repository contains the computational pipeline for quantifying long-term vegetation dynamics across European rewilding sites. By extracting and rigorously processing MODIS time-series data, it computes and statistically tests annual vegetation productivity metrics (INDVI, minNDVI, maxNDVI).

## Repository Contents

*   **`MODIS_Demo.ipynb`**: Initial prototyping notebook demonstrating the extraction of MODIS 250m inputs, the application of quality-assurance (QA) masking, and temporal smoothing based on the framework by Garonna et al. (2009).
*   **`NDVI_Trend_Analysis.ipynb`**: The main analysis pipeline, processing a `FeatureCollection` of rewilding site polygons at two resolutions in parallel:
    *   **250 m** (MOD13Q1, 16-day) — QA-filtered and Garonna-smoothed
    *   **1000 m** (MOD13A3, monthly) — raw scaled NDVI (Williams et al. sensitivity test)
    *   At each resolution, the pipeline:
        *   Applies Dynamic World V1 land-cover masking (full-year 2022) to filter out non-vegetated pixels
        *   Computes annual INDVI, minNDVI, maxNDVI per pixel per year
        *   Runs pixel-wise Mann-Kendall significance testing with tie-corrected variance (Sen, 1968) and continuity-corrected Z-scores
        *   Calculates Sen's Slope (robust median rate of change)
        *   Exports a combined multi-band `.tif` per site to Google Drive
    *   **Exported bands** (Float32, identical structure at both resolutions):
        *   **MK summary** (9 bands): `MK_S_INDVI`, `PValue_INDVI`, `SenSlope_INDVI`, `MK_S_minNDVI`, `PValue_minNDVI`, `SenSlope_minNDVI`, `MK_S_maxNDVI`, `PValue_maxNDVI`, `SenSlope_maxNDVI`
        *   **Annual metrics** (3 × N_years bands): `INDVI_YYYY`, `minNDVI_YYYY`, `maxNDVI_YYYY` for each year in the site's time range
    *   **Demo site** (Wild Ennerdale): additionally exports full temporal NDVI stacks (every 16-day/monthly observation) for inspection
    *   **Export folders**: `GEE_MK_250m/` and `GEE_MK_1000m/`
    *   **Time range**: site-specific start year through 2022 (last year of available climate data)
*   **`visualise_mk_results.R`**: R visualisation script using `tidyterra`, `ggplot2`, and `cowplot`. Generates a publication-quality 4×5 multi-panel figure of Sen's Slope (INDVI) across all 20 rewilding sites:
    *   Panels are ordered by rewilding start year (2000–2013) with year ranges in titles
    *   All panels are equal-sized squares with consistent 1 km scale bars (bottom-left)
    *   Shared diverging colour scale (RdYlGn) with symmetric limits based on the 2nd/98th percentile
    *   Site metadata (start years) hardcoded from the GEE `FeatureCollection` properties
*   **`visualise_ndvi_spiral.R`**: R visualisation script generating PNAS-style seasonality spirals and temporal trend plots using `ggplot2` and `cowplot`.
    *   Creates point-based Cartesian spirals with angular jitter for multi-pixel density visualization (resolving polar geometric interpolation artifacts).
    *   Generates a 2-column composite: Spiral plot (A) alongside iNDVI, minNDVI, and maxNDVI temporal trajectories (B-D).
    *   Contains standalone diagnostic logic for tracing individual pixel phenology loops.

## Outputs

*   **`Outputs/GEE_MK_250m/`**: 250 m GeoTIFFs — MK results + annual metrics per site, plus `NDVI_TimeSeries_250m_Wild_Ennerdale.tif`
*   **`Outputs/GEE_MK_1000m/`**: 1000 m GeoTIFFs — MK results + annual metrics per site, plus `NDVI_TimeSeries_1000m_Wild_Ennerdale.tif`
*   **`Outputs/plots/SenSlope_INDVI_multipanel.png`**: Multi-panel figure produced by `visualise_mk_results.R`
*   **`Outputs/plots/NDVI_spiral_composite_Wild_Ennerdale_250m.png`**: Composite visualization (250m resolution, ±8 day jitter).
*   **`Outputs/plots/NDVI_spiral_composite_Wild_Ennerdale_1000m.png`**: Composite visualization (1000m resolution, ±15 day jitter).
*   **`Outputs/plots/NDVI_spiral_Wild_Ennerdale_..._single_pixel.png`**: Single-pixel phenology diagnostic plots.

## Functions & Unit Tests

All substantive GEE pipeline logic is defined in top-level functions (Cell 2). The notebook executes 18 unit tests (Cell 3) covering every function before any data processing begins. The `visualise_ndvi_spiral.R` logic is fully functionalized and tested in `test_visualise_ndvi_spiral.R`.

| # | Function | Test |
|---|----------|------|
| 1 | `mk_sign` | +1, -1, 0 for synthetic pairs |
| 2–4 | `mk_statistic` | S=10 (monotonic up), S=-10 (down), S=0 (flat) |
| 5 | `mk_variance` | No-tie variance matches analytical formula |
| 6 | `mk_z_score` | Continuity-corrected Z matches hand calculation |
| 7 | `mk_p_value` | p < 0.05 for strong monotonic trend |
| 8 | `ee_cdf` | CDF(0) = 0.5 |
| 9–10 | `sens_slope` | Slope ≈ 0.1 (linear), slope ≈ 0 (flat) |
| 11 | `build_dw_mask` | Returns binary single-band image |
| 12 | `scale_only` | Integer 5000 → float 0.5 |
| 13 | `prep_modis` | Good QA passes, cloud-shadow QA masks |
| 14 | `compute_annual_metrics` | INDVI, minNDVI, maxNDVI correctly aggregated |
| 15 | `garonna_smooth` | Returns correct structure, preserves collection size |
| 16 | `stack_mk` | 9 bands with correct names |
| 17 | `stack_annual` | 3 × N_years bands, correct naming pattern |
| 18 | `stack_temporal` | Band count matches collection size, NDVI_YYYYMMDD naming |
| 19 | `polar_to_cartesian` | (R) Computes correct circle geometry with/without uniform jitter |
| 20 | `interpolate_arcs` | (R) Handles intra-year segments and >50 day gaps via NA breaks |
| 21 | `make_grid` | (R) Correctly generates 4 rings and 12 spokes for polar visualization |

## Figure Captions

**Figure 1: NDVI seasonality and long-term change at an example rewilding site: Wild Ennerdale (250m resolution, 2003–2022).** (A) Seasonality spiral depicting the raw, pixel-level NDVI trajectory. Data points represent individual 16-day MODIS (MOD13Q1) composite observations, with angular jitter (±8 days) applied. Radial distance shows NDVI, and point color shows year. (B–D) Annual trajectories for pixel-level integrated NDVI (iNDVI), minimum NDVI, and maximum NDVI. Semi-transparent lines represent individual pixel trajectories; solid black lines show an estimated mean (LOESS) across the site.

**Figure 2: NDVI seasonality and long-term change at an example rewilding site: Wild Ennerdale (1000m resolution, 2003–2022).** (A) Seasonality spiral depicting the raw, pixel-level NDVI trajectory. Data points represent individual monthly MODIS (MOD13A3) composite observations, with an expanded angular jitter (±15 days) applied. Radial distance shows NDVI, and point color shows year. Point and line opacities are increased to account for the reduced spatial density (16x fewer pixels). (B–D) Annual trajectories for pixel-level integrated NDVI (iNDVI), minimum NDVI, and maximum NDVI. Semi-transparent lines represent individual pixel trajectories; solid black lines show an estimated mean (LOESS) across the site. Note the dampening of extreme values compared to the 250m resolution due to spatial aggregation.

## Agentic Programming & Refactoring

This repository was constructed using agentic programming (Claude and Gemini with human supervision and verification). Fully human-authored Javascript GEE analyses were refactored to the Python API by AI, all key code logic is unit tested and human-validated. Accountability for any errors is with the human authors.
