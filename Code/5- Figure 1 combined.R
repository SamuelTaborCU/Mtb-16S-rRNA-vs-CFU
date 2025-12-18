#######################################
##
## Script name: Figure 1 combined (Untreated + Single drugs + Regimens)
##
## Purpose of script: Combine pre-generated manuscript panels into a single multi-panel Figure 1
##                    and save as a high-resolution JPG with (a–g) labels.
##
## Author: Samuel Tabor
##
#######################################

###############################################
# Import libraries / Define paths

library(jpeg)
library(grid)

out_file <- "../DataProcessed/Plots/Figure 1 combined.jpg"

img_paths <- c(
  "../DataProcessed/Plots/1-Untreated mice/Untreated_mice_panel_A_timecourse.jpg",                 # a
  "../DataProcessed/Plots/1-Untreated mice/Untreated_mice_panel_B_scatter.jpg",                   # b
  "../DataProcessed/Plots/2-Single drugs/Single_drugs_panel_B_scatter.jpg",                       # c
  "../DataProcessed/Plots/3- Regimen treated/Regimen_CFU vs_16S_28_panel_d_scatter_without_PreRx.jpg", # d
  "../DataProcessed/Plots/2-Single drugs/Single Drugs CFU vs 16S rRNA log reduction from PreRx.jpg",   # e (spans two columns)
  "../DataProcessed/Plots/2-Single drugs/Single Drugs ranking.jpg",                               # f
  "../DataProcessed/Plots/3- Regimen treated/Regimen ranking.jpg"                                 # g
)

panel_labels <- letters[1:length(img_paths)]  # a-g

###############################################
# Read images -> grobs

imgs  <- lapply(img_paths, readJPEG)
grobs <- lapply(imgs, rasterGrob)

###############################################
# Layout + label helper

ncol_layout <- 2
nrow_layout <- 4

add_label <- function(label) {
  grid.text(
    paste0(label, ")"),
    x = unit(2, "mm"),
    y = unit(1, "npc") - unit(0, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 12, col = "black", fontface = "bold")
  )
}

###############################################
# Save combined figure 

jpeg(out_file, width = 1800, height = 2400, res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow_layout, ncol_layout)))

for (i in seq_along(grobs)) {
  
  # e panel spans row 3 across both columns
  if (i == 5) {
    vp <- viewport(layout.pos.row = 3, layout.pos.col = 1:2)
  } else {
    row <- ifelse(i >= 6, 4, ceiling(i / ncol_layout))
    col <- ifelse(i == 6, 1, ifelse(i == 7, 2, i %% ncol_layout))
    col <- ifelse(col == 0, ncol_layout, col)
    vp  <- viewport(layout.pos.row = row, layout.pos.col = col)
  }
  
  pushViewport(vp)
  grid.draw(grobs[[i]])
  add_label(panel_labels[i])
  popViewport()
}

dev.off()
