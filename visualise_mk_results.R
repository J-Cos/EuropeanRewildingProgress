# ==============================================================================
# Batch plot Sen's Slope (INDVI) across all 250m sites
# ==============================================================================

library(terra)

input_dir  <- "Outputs/GEE_MK_250m"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tifs <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract site name from filename
get_site_name <- function(path) {
  base <- tools::file_path_sans_ext(basename(path))
  sub("^MK_Stacked_250m_", "", base) |> gsub("_", " ", x = _)
}

# ── Individual site plots ────────────────────────────────
for (f in tifs) {
  site <- get_site_name(f)
  r <- rast(f)
  NAflag(r) <- -9999

  slope <- r[["SenSlope_INDVI"]]
  names(slope) <- "Sen's Slope (INDVI)"

  out_file <- file.path(output_dir, paste0("SenSlope_INDVI_", gsub(" ", "_", site), ".png"))
  png(out_file, width = 800, height = 700, res = 150)
  plot(slope,
       main = paste0(site, "\nSen's Slope (INDVI)"),
       col = hcl.colors(50, "RdYlGn"),
       axes = FALSE,
       mar = c(2, 2, 3, 4))
  dev.off()
  message(paste("  Plotted:", site))
}

message(paste("\n✅ All", length(tifs), "site plots saved to", output_dir))
