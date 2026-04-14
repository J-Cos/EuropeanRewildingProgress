# ==============================================================================
# Multi-panel figure: Sen's Slope (INDVI) across all 250m sites (4x5)
# ==============================================================================

library(terra)

input_dir  <- "Outputs/GEE_MK_250m"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

tifs <- sort(list.files(input_dir, pattern = "\\.tif$", full.names = TRUE))

get_site_name <- function(path) {
  base <- tools::file_path_sans_ext(basename(path))
  sub("^MK_Stacked_250m_", "", base) |> gsub("_", " ", x = _)
}

# ── Load all slope layers and compute shared range ───────
slopes <- lapply(tifs, function(f) {
  r <- rast(f)
  NAflag(r) <- -9999
  r[["SenSlope_INDVI"]]
})
names(slopes) <- sapply(tifs, get_site_name)

all_vals <- unlist(lapply(slopes, function(r) values(r, na.rm = TRUE)))
zlim <- c(-1, 1) * max(abs(quantile(all_vals, c(0.02, 0.98))), na.rm = TRUE)
ncols <- 50
pal <- hcl.colors(ncols, "RdYlGn")

# ── Helper: draw a 1 km scale bar (black, bottom-left) ───
add_scalebar <- function(r) {
  e <- ext(r)
  bar_deg <- 1000 / (111320 * cos(mean(c(e$ymin, e$ymax)) * pi / 180))
  
  # Fixed proportional position: 5% in from left, 6% up from bottom
  x0 <- e$xmin + (e$xmax - e$xmin) * 0.05
  y0 <- e$ymin + (e$ymax - e$ymin) * 0.06
  
  lines(c(x0, x0 + bar_deg), c(y0, y0), lwd = 2, col = "black")
  text(x0 + bar_deg / 2, y0 + (e$ymax - e$ymin) * 0.04,
       "1 km", cex = 0.45, col = "black", font = 2)
}

# ── Build the layout: 4 rows x 5 cols + 1 legend column ─
# Matrix: 20 panels (1-20) in 4x5, plus column 6 = legend (21)
mat <- matrix(c(
  1,  2,  3,  4,  5, 21,
  6,  7,  8,  9, 10, 21,
  11, 12, 13, 14, 15, 21,
  16, 17, 18, 19, 20, 21
), nrow = 4, byrow = TRUE)

out_file <- file.path(output_dir, "SenSlope_INDVI_multipanel.png")
png(out_file, width = 3000, height = 2400, res = 200)

layout(mat, widths = c(rep(1, 5), 0.45))
par(oma = c(0, 0, 2, 0))  # outer margin for title

for (i in seq_along(slopes)) {
  par(mar = c(1, 1, 2, 0.5))
  r <- slopes[[i]]
  
  plot(r, col = pal, range = zlim, legend = FALSE, axes = FALSE,
       main = names(slopes)[i], cex.main = 0.8, font.main = 1)
  add_scalebar(r)
}

# ── Legend panel ─────────────────────────────────────────
par(mar = c(4, 1, 4, 4.5))
plot.new()
plot.window(xlim = c(0, 1), ylim = zlim)

# Draw colour ramp
n <- length(pal)
yvals <- seq(zlim[1], zlim[2], length.out = n + 1)
for (j in seq_len(n)) {
  rect(0.15, yvals[j], 0.55, yvals[j + 1], col = pal[j], border = NA)
}
rect(0.15, zlim[1], 0.55, zlim[2], border = "grey30", lwd = 0.5)

axis(4, las = 1, cex.axis = 0.75, line = -1.8)
mtext("Sen's Slope\n(INDVI yr\u207b\u00b9)", side = 4, line = 1.8, cex = 0.65)

dev.off()
message(paste("✅ Multi-panel figure saved to", out_file))
