# ==============================================================================
# Multi-panel figure: Sen's Slope (INDVI) across all 250m sites (4x5)
# Ordered by start_year (earliest first)
# ==============================================================================

library(terra)

input_dir  <- "Outputs/GEE_MK_250m"
output_dir <- "Outputs/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Site metadata (from GEE FeatureCollection) ───────────
meta <- data.frame(
  site_id = c(
    "Alladale", "Border_Meuse", "Central_Apennines", "Chernobyl_Exclusion_Zone",
    "Creag_Meagaidh", "Doberitzer_Heide", "Greater_Coa_Valley", "Naliboki_Forest",
    "Oostvaardersplassen", "Pentezug", "Swiss_National_Park",
    "Knepp", "Wild_Ennerdale", "Tarutino_Steppe",
    "Dundreggan", "Kempen_Broek",
    "Southern_Carpathians", "Velebit_Mountains", "Vanatori_Neamt_Nature_Park",
    "Danube_Delta"
  ),
  start = c(
    2000, 2000, 2000, 2000,
    2000, 2000, 2000, 2000,
    2000, 2000, 2000,
    2001, 2003, 2004,
    2008, 2010,
    2011, 2012, 2012,
    2013
  ),
  end = 2025,
  stringsAsFactors = FALSE
)
meta <- meta[order(meta$start, meta$site_id), ]

# ── Load slope layers in order ───────────────────────────
slopes <- list()
titles <- character()
for (i in seq_len(nrow(meta))) {
  sid <- meta$site_id[i]
  fname <- file.path(input_dir, paste0("MK_Stacked_250m_", sid, ".tif"))
  if (!file.exists(fname)) { message("  Skipping: ", sid); next }
  r <- rast(fname)
  NAflag(r) <- -9999
  slopes[[length(slopes) + 1]] <- r[["SenSlope_INDVI"]]
  label <- gsub("_", " ", sid)
  titles <- c(titles, paste0(label, " (", meta$start[i], "\u2013", meta$end[i], ")"))
}

all_vals <- unlist(lapply(slopes, function(r) values(r, na.rm = TRUE)))
zlim <- c(-1, 1) * max(abs(quantile(all_vals, c(0.02, 0.98))), na.rm = TRUE)
pal <- hcl.colors(50, "RdYlGn")

# ── Helper: 1 km scale bar using plot usr coords ────────
add_scalebar <- function() {
  u <- par("usr")
  lat_mid <- (u[3] + u[4]) / 2
  bar_deg <- 1000 / (111320 * cos(lat_mid * pi / 180))
  x0 <- u[1] + (u[2] - u[1]) * 0.05
  y0 <- u[3] + (u[4] - u[3]) * 0.08
  lines(c(x0, x0 + bar_deg), c(y0, y0), lwd = 2, col = "black")
  text(x0 + bar_deg / 2, y0 + (u[4] - u[3]) * 0.04,
       "1 km", cex = 0.45, col = "black", font = 2)
}

# ── Layout: 4x5 + legend ────────────────────────────────
mat <- matrix(0, nrow = 4, ncol = 6)
for (idx in seq_along(slopes)) {
  mat[((idx - 1) %/% 5) + 1, ((idx - 1) %% 5) + 1] <- idx
}
mat[, 6] <- length(slopes) + 1

out_file <- file.path(output_dir, "SenSlope_INDVI_multipanel.png")
png(out_file, width = 3000, height = 2400, res = 200)
layout(mat, widths = c(rep(1, 5), 0.45))
par(oma = c(0, 0, 2, 0))

for (i in seq_along(slopes)) {
  par(mar = c(1, 1, 2.5, 0.5))
  plot(slopes[[i]], col = pal, range = zlim, legend = FALSE, axes = FALSE,
       main = titles[i], cex.main = 0.65, font.main = 1)
  add_scalebar()
}

# Empty panels
empties <- 4 * 5 - length(slopes)
for (e in seq_len(max(empties, 0))) { par(mar = c(1, 1, 2.5, 0.5)); plot.new() }

# ── Legend ───────────────────────────────────────────────
par(mar = c(4, 1, 4, 4.5))
plot.new()
plot.window(xlim = c(0, 1), ylim = zlim)
n <- length(pal)
yvals <- seq(zlim[1], zlim[2], length.out = n + 1)
for (j in seq_len(n)) rect(0.15, yvals[j], 0.55, yvals[j + 1], col = pal[j], border = NA)
rect(0.15, zlim[1], 0.55, zlim[2], border = "grey30", lwd = 0.5)
axis(4, las = 1, cex.axis = 0.75, line = -1.8)
mtext("Sen's Slope\n(INDVI yr\u207b\u00b9)", side = 4, line = 1.8, cex = 0.65)

dev.off()
message(paste("\u2705 Multi-panel figure saved to", out_file))
