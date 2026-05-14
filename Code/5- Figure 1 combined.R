#######################################
##
## Script name: Figure 1 combined (Untreated + Single drugs + Regimens)
##
## Purpose of script: Combine pre-generated manuscript panels into a single multi-panel Figure 1
##                    and save as high-resolution TIFF and EPS with (a–g) labels.
##
## Author: Samuel Tabor
##
#######################################

###############################################
# Import libraries / Define paths

library(tiff)
library(grid)

out_file_tiff <- "../DataProcessed/Plots/Figure 1 combined.tiff"
out_file_eps  <- "../DataProcessed/Plots/Figure 1 combined.eps"

dir.create(dirname(out_file_tiff), recursive = TRUE, showWarnings = FALSE)

img_paths <- c(
  "../DataProcessed/Plots/1-Untreated mice/Untreated_mice_panel_A_timecourse.tiff",                       # a
  "../DataProcessed/Plots/1-Untreated mice/Untreated_mice_panel_B_scatter.tiff",                         # b
  "../DataProcessed/Plots/2-Single drugs/Single_drugs_panel_B_scatter.tiff",                             # c
  "../DataProcessed/Plots/3- Regimen treated/Regimen_CFU vs_16S_28_panel_d_scatter_without_PreRx.tiff", # d
  "../DataProcessed/Plots/2-Single drugs/Single Drugs CFU vs 16S rRNA log Reduction from PreRx.tiff",   # e
  "../DataProcessed/Plots/2-Single drugs/Single Drugs ranking.tiff",                                     # f
  "../DataProcessed/Plots/3- Regimen treated/Regimen ranking.tiff"                                       # g
)

panel_labels <- letters[1:length(img_paths)]  # a-g

###############################################
# Read images -> grobs

imgs  <- lapply(img_paths, function(x) readTIFF(x, native = TRUE))
grobs <- lapply(imgs, function(x) rasterGrob(x, interpolate = TRUE))

###############################################
# Layout + label helper
# Add a little space between ALL plot rows
# Move a,b,c... slightly higher

ncol_layout <- 2

# 4 plot rows + 3 spacer rows between them
layout_heights <- unit.c(
  unit(1, "null"), unit(3, "mm"),
  unit(1, "null"), unit(3, "mm"),
  unit(1, "null"), unit(3, "mm"),
  unit(1, "null")
)

add_label <- function(label) {
  grid.text(
    paste0(label, ")"),
    x = unit(1.5, "mm"),
    y = unit(1, "npc") + unit(0.0, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 12, col = "black", fontface = "bold")
  )
}

draw_combined_figure <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 7,
        ncol = ncol_layout,
        heights = layout_heights
      )
    )
  )
  
  for (i in seq_along(grobs)) {
    
    # Plot rows are 1, 3, 5, 7 because 2, 4, 6 are spacer rows
    if (i == 5) {
      # Panel e spans both columns
      vp <- viewport(layout.pos.row = 5, layout.pos.col = 1:2)
    } else {
      row <- ifelse(i >= 6, 7, c(1, 1, 3, 3)[i])
      col <- ifelse(i == 6, 1, ifelse(i == 7, 2, i %% ncol_layout))
      col <- ifelse(col == 0, ncol_layout, col)
      vp  <- viewport(layout.pos.row = row, layout.pos.col = col)
    }
    
    pushViewport(vp)
    
    pushViewport(viewport(
      x = 0.5, y = 0.48,
      width = 1, height = 0.96,
      just = c("center", "center")
    ))
    grid.draw(grobs[[i]])
    popViewport()
    
    add_label(panel_labels[i])
    popViewport()
  }
}

###############################################
# Save combined figure as high-resolution TIFF

tiff(
  filename = out_file_tiff,
  width = 6,
  height = 8.3,
  units = "in",
  res = 600,
  compression = "lzw"
)

draw_combined_figure()
dev.off()

###############################################
# Save combined figure as EPS

cairo_ps(
  filename = out_file_eps,
  width = 6,
  height = 8.3,
  onefile = FALSE,
  fallback_resolution = 600
)

draw_combined_figure()
dev.off()

