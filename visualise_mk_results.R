# ==============================================================================
# Simple visualisor for Earth Engine Mann-Kendall Stacked Rasters
# ==============================================================================

library(tidyverse)
library(terra)
library(tidyterra)

# 1. Provide the path to the downloaded .tif file from Google Drive
# Replace this with your actual downloaded file path
raster_path <- "MannKendall_Stacked_Knepp.tif"

if(!file.exists(raster_path)){
  message("Please provide the correct path to your downloaded .tif file on line 10.")
} else {
  
  # 2. Load the stacked multi-band raster
  r <- rast(raster_path)
  
  # Print the imported raster to see the available 9 bands:
  # MK_S_INDVI, PValue_INDVI, SenSlope_INDVI, etc.
  print(r)
  
  # 3. Choose a metric to visualise (e.g., INDVI)
  # We extract only the 3 relevant bands for an easier plot
  metric <- "INDVI" 
  r_metric <- r[[paste0(c("MK_S_", "PValue_", "SenSlope_"), metric)]]
  
  # Rename them for cleaner facet labels
  names(r_metric) <- c("Mann-Kendall S", "P-Value", "Sen's Slope")
  
  # 4. Visualisation
  # This sets up a grid with independent colour scales for each statistic
  p <- ggplot() +
    geom_spatraster(data = r_metric) +
    facet_wrap(~lyr, scales = "free") +
    scale_fill_whitebox_c(
      palette = "muted", 
      na.value = "transparent"
    ) +
    theme_minimal() +
    labs(
      title = paste("NDVI Mann-Kendall Trend Results:", metric),
      subtitle = "European Rewilding Progress Space",
      fill = "Value"
    ) +
    theme(
      axis.text = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid = element_blank()
    )
  
  # Print the plot
  print(p)
  
  # 5. Save the plot to disk
  ggsave(
    filename = paste0("MK_results_plot_", metric, ".png"), 
    plot = p, 
    width = 12, height = 5, dpi = 300, bg = "white"
  )
  
  message("Plot saved to working directory.")
}
