#######################################
##
## Script name: RMM Figure 4 combined (T95 + predicted vs observed panels)
##
## Purpose of script: Combine four pre-rendered tiff panels into a single 2x2 figure
##                    and save as a high-resolution tiff with panel labels (a–d).
##
## Author: Samuel Tabor
##
#######################################

###############################################
####################STARTUP####################
###############################################

# Import libraries
# Set options
# Define variables

options(stringsAsFactors = FALSE)

library(tiff)
library(grid)

in_paths <- c(
  "../DataProcessed/Plots/5- RMM/T95 of the regimens.tiff",
  "../DataProcessed/Plots/5- RMM/Model 1 T95 predicted vs observed.tiff",
  "../DataProcessed/Plots/5- RMM/Model 2 T95 predicted vs observed.tiff",
  "../DataProcessed/Plots/5- RMM/Model 3 T95 predicted vs observed.tiff"
)

out_path <- "../DataProcessed/Plots/Figure 4 combined.tiff"

###############################################
#################FIGURE ASSEMBLY###############
###############################################

imgs  <- lapply(in_paths, function(x) readTIFF(x, native = TRUE))
grobs <- lapply(imgs, rasterGrob)

tiff(
  filename = out_path,
  width = 6,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.newpage()
pushViewport(viewport(
  layout = grid.layout(
    nrow = 2, ncol = 2,
    widths  = unit(rep(1, 2), "null"),
    heights = unit(rep(1, 2), "null")
  )
))

draw_panel <- function(g, r, c, label) {
  pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
  grid.draw(g)
  grid.text(
    paste0(label, ")"),
    x = unit(2, "mm"),
    y = unit(1, "npc") - unit(2, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 12, col = "black", fontface = "bold")
  )
  popViewport()
}

draw_panel(grobs[[1]], 1, 1, "a")
draw_panel(grobs[[2]], 1, 2, "b")
draw_panel(grobs[[3]], 2, 1, "c")
draw_panel(grobs[[4]], 2, 2, "d")

dev.off()

out_path_eps <- "../DataProcessed/Plots/Figure 4 combined.eps"

cairo_ps(
  filename = out_path_eps,
  width = 6,
  height = 6,
  onefile = FALSE,
  fallback_resolution = 600
)

grid.newpage()
pushViewport(viewport(
  layout = grid.layout(
    nrow = 2, ncol = 2,
    widths  = unit(rep(1, 2), "null"),
    heights = unit(rep(1, 2), "null")
  )
))

draw_panel(grobs[[1]], 1, 1, "a")
draw_panel(grobs[[2]], 1, 2, "b")
draw_panel(grobs[[3]], 2, 1, "c")
draw_panel(grobs[[4]], 2, 2, "d")

dev.off()
