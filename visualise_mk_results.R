# ==============================================================================
# Multi-panel figure: Sen's Slope (INDVI) across all 250m sites (4x5)
# Reads site_metadata.csv for year ranges; orders by start_year
# ==============================================================================

library(terra)

input_dir  <- "Outputs/GEE_MK_250m"
output_dir <- "Outputs/plots"
meta_file  <- "Outputs/site_metadata.csv"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load site metadata and sort by start_year ────────────
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[order(meta$start_year, meta$site_id), ]

# ── Load slope layers in metadata order ──────────────────
slopes <- list()
titles <- character()
for (i in seq_len(nrow(meta))) {
  sid <- meta$site_id[i]
  fname <- file.path(input_dir, paste0("MK_Stacked_250m_", sid, ".tif"))
  if (!file.exists(fname)) {
    message("  Skipping (file not found): ", sid)
    next
  }
  r <- rast(fname)
  NAflag(r) <- -9999
  slopes[[length(slopes) + 1]] <- r[["SenSlope_INDVI"]]
  label <- gsub("_", " ", sid)
  titles <- c(titles, paste0(label, "\n(", meta$start_year[i], "\u2013", meta$end_year[i], ")"))
}

all_vals <- unlist(lapply(slopes, function(r) values(r, na.rm = TRUE)))
zlim <- c(-1, 1) * max(abs(quantile(all_vals, c(0.02, 0.98))), na.rm = TRUE)
ncols <- 50
pal <- hcl.colors(ncols, "RdYlGn")

# ── Helper: 1 km scale bar in absolute plot coords ──────
add_scalebar <- function(r) {
  # Use usr coords (plot region) for consistent placement
  u <- par("usr")  # c(x1, x2, y1, y2)
  lat_mid <- (u[3] + u[4]) / 2
  bar_deg <- 1000 / (111320 * cos(lat_mid * pi / 180))

  # 5% in from left edge, 8% up from bottom edge
  x0 <- u[1] + (u[2] - u[1]) * 0.05
  y0 <- u[3] + (u[4] - u[3]) * 0.08

  lines(c(x0, x0 + bar_deg), c(y0, y0), lwd = 2, col = "black", xpd = FALSE)
  text(x0 + bar_deg / 2, y0 + (u[4] - u[3]) * 0.04,
       "1 km", cex = 0.45, col = "black", font = 2)
}

# ── Layout: 4x5 + legend column ─────────────────────────
n_panels <- length(slopes)
n_cols <- 5
n_rows <- ceiling(n_panels / n_cols)

mat_vec <- seq_len(n_panels)
# Pad with 0 (empty) if fewer than n_rows * n_cols
if (length(mat_vec) < n_rows * n_cols) {
  mat_vec <- c(mat_vec, rep(0, n_rows * n_cols - length(mat_vec)))
}
mat <- matrix(c(rbind(matrix(mat_vec, nrow = n_rows, byrow = TRUE),
                      rep(n_panels + 1, n_rows))),
              nrow = n_rows, byrow = FALSE)
# Simpler: build manually
mat <- matrix(0, nrow = n_rows, ncol = n_cols + 1)
for (idx in seq_len(n_panels)) {
  row <- ((idx - 1) %/% n_cols) + 1
  col <- ((idx - 1) %% n_cols) + 1
  mat[row, col] <- idx
}
mat[, n_cols + 1] <- n_panels + 1  # legend column

out_file <- file.path(output_dir, "SenSlope_INDVI_multipanel.png")
png(out_file, width = 3000, height = 2400, res = 200)

layout(mat, widths = c(rep(1, n_cols), 0.45))
par(oma = c(0, 0, 2, 0))

for (i in seq_along(slopes)) {
  par(mar = c(1, 1, 2.5, 0.5))
  r <- slopes[[i]]

  plot(r, col = pal, range = zlim, legend = FALSE, axes = FALSE,
       main = titles[i], cex.main = 0.7, font.main = 1)
  add_scalebar(r)
}

# Fill empty panels if any
empties <- n_rows * n_cols - n_panels
if (empties > 0) {
  for (e in seq_len(empties)) {
    par(mar = c(1, 1, 2.5, 0.5))
    plot.new()
  }
}

# ── Legend panel ─────────────────────────────────────────
par(mar = c(4, 1, 4, 4.5))
plot.new()
plot.window(xlim = c(0, 1), ylim = zlim)

n <- length(pal)
yvals <- seq(zlim[1], zlim[2], length.out = n + 1)
for (j in seq_len(n)) {
  rect(0.15, yvals[j], 0.55, yvals[j + 1], col = pal[j], border = NA)
}
rect(0.15, zlim[1], 0.55, zlim[2], border = "grey30", lwd = 0.5)

axis(4, las = 1, cex.axis = 0.75, line = -1.8)
mtext("Sen's Slope\n(INDVI yr\u207b\u00b9)", side = 4, line = 1.8, cex = 0.65)

dev.off()
message(paste("\u2705 Multi-panel figure saved to", out_file))
