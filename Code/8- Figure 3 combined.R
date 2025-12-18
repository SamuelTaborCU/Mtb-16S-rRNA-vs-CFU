#######################################
##
## Script name: RMM Figure 3 combined (T95 + predicted vs observed panels)
##
## Purpose of script: Combine four pre-rendered JPEG panels into a single 2x2 figure
##                    and save as a high-resolution JPEG with panel labels (a–d).
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

library(jpeg)
library(grid)

in_paths <- c(
  "../DataProcessed/Plots/5- RMM/T95 of the regimens.jpg",
  "../DataProcessed/Plots/5- RMM/Model 1 T95 predicted vs observed.jpg",
  "../DataProcessed/Plots/5- RMM/Model 2  T95 predicted vs observed.jpg",
  "../DataProcessed/Plots/5- RMM/Model 3  T95 predicted vs observed.jpg"
)

out_path <- "../DataProcessed/Plots/Figure 3 combined.jpg"

###############################################
#################FIGURE ASSEMBLY###############
###############################################

imgs  <- lapply(in_paths, readJPEG)
grobs <- lapply(imgs, rasterGrob)

jpeg(out_path, width = 1800, height = 1800, res = 300)

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

