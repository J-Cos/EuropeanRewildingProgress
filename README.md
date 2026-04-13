# European Rewilding Progress

This repository contains the computational pipeline for quantifying long-term vegetation dynamics across European rewilding sites. By extracting and rigorously processing MODIS time-series data, it computes and statistically tests annual vegetation productivity metrics (INDVI, minNDVI, maxNDVI).

## Repository Contents

*   **`MODIS_Demo.ipynb`**: Initial prototyping notebook demonstrating the extraction of MODIS 250m inputs, the application of quality-assurance (QA) masking, and temporal smoothing based on the framework by Garonna et al. (2009).
*   **`NDVI_Trend_Analysis.ipynb`**: A robust, scaleable pipeline iterating over a complex `ee.FeatureCollection` of rewilding polygons. This pipeline:
    *   Synthesises annual metrics alongside Dynamic World V1 static masks to filter out undesirable land-cover (e.g. water, crops, built environments).
    *   Applies pixel-wise, monotonic Mann-Kendall significance testing with tie-corrected variance (Sen, 1968) and continuity-corrected Z-scores.
    *   Calculates a robust median rate of change using Sen's Slope.
    *   Exports a dense, 9-band `.tif` stack mapping metric statistics directly to Google Drive. The exported raster contains the following Float32 bands:
        *   `MK_S_INDVI`, `PValue_INDVI`, `SenSlope_INDVI`
        *   `MK_S_minNDVI`, `PValue_minNDVI`, `SenSlope_minNDVI`
        *   `MK_S_maxNDVI`, `PValue_maxNDVI`, `SenSlope_maxNDVI`
*   **`visualise_mk_results.R`**: A lightweight R plotting utility using base `terra::plot` to natively handle the geospatial rendering of the exported Earth Engine multi-band grids, applying independent sequential and diverging aesthetics to the generated S, P-value, and slope matrices.

## Agentic Programming & Refactoring

This repository was constructed using agentic programming (Claude and Gemini with human supervision and verification). Fully human-authored Javascript GEE analyses were refactored to the Python API by AI, all key code logic is unit tested and human-validated. Accountability for any errors is with the human authors.
