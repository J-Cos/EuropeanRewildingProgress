# ==============================================================================
# Simple visualisor for Earth Engine Mann-Kendall Stacked Rasters
# ==============================================================================

library(terra)

# 1. Provide the path to the downloaded .tif file
raster_path <- "MannKendall_Stacked_Knepp.tif"

if(!file.exists(raster_path)){
  message("Please provide the correct path to your downloaded .tif file on line 8.")
} else {
  
  # 2. Load the stacked multi-band raster
  r <- rast(raster_path)
  
  # Choose a metric to visualise (e.g., INDVI)
  metric <- "INDVI" 
  
  # Extract individual bands
  r_mk    <- r[[paste0("MK_S_", metric)]]
  r_p     <- r[[paste0("PValue_", metric)]]
  r_slope <- r[[paste0("SenSlope_", metric)]]
  
  # Rename for plot headings
  names(r_mk) <- "Mann-Kendall S"
  names(r_p) <- "P-Value"
  names(r_slope) <- "Sen's Slope"
  
  # Combine them back into a slim stack
  r_plot <- c(r_mk, r_p, r_slope)
  
  # Set up output file
  out_file <- paste0("MK_results_plot_", metric, ".png")
  png(out_file, width = 1600, height = 500, res = 150)
  
  # Plot side-by-side. 'terra::plot' natively handles independent colour scales
  # per layer in the plot stack!
  plot(r_plot, 
       main = names(r_plot),
       mar = c(3, 3, 3, 5),
       axes = FALSE)
       
  dev.off()
  
  message(paste("Plot saved to:", out_file))
}
