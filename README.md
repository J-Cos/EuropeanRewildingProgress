# European Rewilding Progress

This repository contains the computational pipeline for quantifying long-term vegetation dynamics across European rewilding sites. By extracting and rigorously processing MODIS time-series data, it computes and statistically tests annual vegetation productivity metrics (INDVI, minNDVI, maxNDVI).

## Repository Contents

*   **`MODIS_Demo.ipynb`**: Initial prototyping notebook demonstrating the extraction of MODIS 250m inputs, the application of quality-assurance (QA) masking, and temporal smoothing based on the framework by Garonna et al. (2009).
*   **`NDVI_Trend_Analysis.ipynb`**: The main analysis pipeline, processing a `FeatureCollection` of rewilding site polygons at two resolutions in parallel:
    *   **250 m** (MOD13Q1, 16-day) — QA-filtered and Garonna-smoothed
    *   **1000 m** (MOD13A3, monthly) — raw scaled NDVI (Williams et al. sensitivity test)
    *   At each resolution, the pipeline:
        *   Applies Dynamic World V1 land-cover masking to filter out undesirable pixels (water, crops, built, bare, snow/ice)
        *   Runs pixel-wise Mann-Kendall significance testing with tie-corrected variance (Sen, 1968) and continuity-corrected Z-scores
        *   Calculates Sen's Slope (robust median rate of change)
        *   Exports a 9-band stacked `.tif` per site to Google Drive
    *   **Exported bands** (Float32, identical structure at both resolutions):
        *   `MK_S_INDVI`, `PValue_INDVI`, `SenSlope_INDVI`
        *   `MK_S_minNDVI`, `PValue_minNDVI`, `SenSlope_minNDVI`
        *   `MK_S_maxNDVI`, `PValue_maxNDVI`, `SenSlope_maxNDVI`
    *   **Export folders**: `GEE_MK_250m/` and `GEE_MK_1000m/`
*   **`visualise_mk_results.R`**: A lightweight R plotting utility using base `terra::plot` to render the exported multi-band GeoTIFFs with independent colour scales per statistic.

## Functions & Unit Tests

All substantive logic is defined in top-level functions (Cell 2). The notebook executes 15 unit tests (Cell 3) covering every function before any data processing begins:

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

## Agentic Programming & Refactoring

This repository was constructed using agentic programming (Claude and Gemini with human supervision and verification). Fully human-authored Javascript GEE analyses were refactored to the Python API by AI, all key code logic is unit tested and human-validated. Accountability for any errors is with the human authors.
