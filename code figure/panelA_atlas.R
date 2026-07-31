sample_names_fig <- c("M70", "I1", "L1", "L15")

col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2-CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

seumerged_sansvoid <- subset(x = seumerged, subset = regions_md != "void")
for (smpl in sample_names_fig){
  df <- seumerged_sansvoid@meta.data %>% 
    filter(orig.ident == smpl)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)
  
  p <-  df %>% 
    ggplot(aes(x = col, y = row, color = regions_md))+
    geom_point(size = 1.9)+
    scale_color_manual(name = "Regions", values = col_regions) +
    scale_y_reverse() +
    coord_fixed() +  # keep aspect ratio equal
    theme_minimal() +
    theme(
      axis.text = element_blank(),         # Remove axis numbers
      axis.ticks = element_blank(),        # Remove axis ticks
      axis.title = element_blank(),        # Remove axis titles
      panel.grid = element_blank(),        # Remove grid lines
      panel.border = element_blank()  # turn off default border
      # legend.position = "none"
    ) +
    #+ labs(color = "Regions")
    + guides(
      color = guide_legend(
        override.aes = list(shape = 15, size = 5)
      )
    ) +  
    theme(legend.box.background = element_rect(color = "black", size = .25)) # Add a box around the legend

  pdf_file <- file.path(fig_path,"figure2_panels",paste0("panelA_", smpl,"_atlas_regions_md.pdf"))
  pdf(pdf_file, height = 5, width = 4)
  print(p) #needs to be explicitely printed
  dev.off()
}

legend_atlas <- get_legend(p)

pdf_file <- file.path(fig_path,"figure2_panels","legend_atlas.pdf")
pdf(pdf_file, height = 5, width = 4, useDingbats = FALSE)
grid::grid.draw(legend_atlas) #needs to be explicitely printed
dev.off()