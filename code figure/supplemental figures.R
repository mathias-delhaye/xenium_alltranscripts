## Code to save dot plots of case images for regions and regions_md

col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300", "grey"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void", "hilus"))

color_ILAE <- setNames(c("#C67BFF","#F8766C","#7BAD00","#00C0C5"),
                       c("P", "ILAE 0","ILAE 1","ILAE 2"))

save_path <- file.path(fig_path,"figure2/panels_mainfigure/panelA")

dir.create(path = save_path, recursive = T)

# name of all samples in megadata in order assigned in megadata
sample_names <- names(megadata@images)
panelA_regions <- c("I1", "L1", "L15", "M70")

# loading files
seumerged <- read_rds(file.path(data_path, "seumerged_log30_v6.rds"))
# removing void (and NAs)
seumerged_sansvoid <- subset(x = seumerged, subset = regions != "void")

for (smpl in panelA_regions){
  pdf(file.path(save_path,paste0("panelA_", smpl,"_atlas_regions_md.pdf")), width = 4, height = 5)
  
  df <- seumerged_sansvoid@meta.data %>% 
    filter(orig.ident == smpl)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)
  
   p <-  df %>% 
    ggplot(aes(x = col, y = row, color = regions_md))+
    geom_point(shape = 15, size = 2.5)+
    scale_color_manual(name = "Regions", values = col_regions) +
    scale_y_reverse() +
    coord_fixed() +  # keep aspect ratio equal
    theme_minimal() +
    theme(
      axis.text = element_blank(),         # Remove axis numbers
      axis.ticks = element_blank(),        # Remove axis ticks
      axis.title = element_blank(),        # Remove axis titles
      panel.grid = element_blank(),        # Remove grid lines
      panel.border = element_blank(),  # turn off default border
      legend.position = "none"
    )
    # + labs(color = "Regions")
    # + guides(
    #   color = guide_legend(
    #     override.aes = list(shape = 15, size = 5)
    #   )
    # ) +  
    #theme(legend.box.background = element_rect(color = "black", size = .25)) # Add a box around the legend
  print(p) #needs to be explicitely printed
  dev.off()
}

# removing void (and NAs)
seumerged_CAs <- subset(x = seumerged, subset = regions %in% c("CA2_CA4","CA1","GABA","SL.SR.SLM"))

# calculate the proportion of each subregion and each case
proportion_CAs <- feature_proportion(seumerged_CAs, feature = "smpl_ILAE", loc = "regions_md") %>% 
  separate(smpl_ILAE, into = c("smpl", "ILAE"),sep = "_") %>%
  mutate(ILAE = case_when(
    ILAE == "0" ~ "ILAE 0",
    ILAE == "1" ~ "ILAE 1",
    ILAE == "2" ~ "ILAE 2",
    TRUE ~ ILAE
  ),
  ILAE = factor(ILAE, c( "P", "ILAE 0", "ILAE 1", "ILAE 2")))

# region we want to plot
boxplot_comparison <- function(df, loc = "CA1", ttl = "CA1", y_limit = c(0,50), s = 1){
  
  # isolate for the specific region
  df_sub <- filter(df, regions_md %in% loc) %>% 
    group_by(smpl, ILAE) %>% 
    summarise(avg_percent_reg = sum(avg_percent_reg))
  
  dunn <- dunnTest(avg_percent_reg ~ ILAE, data = df_sub, method = "bh")
  dunn <- dunn$res %>% 
    mutate(significance = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "**",
      TRUE ~ "ns"
    )) %>% 
    separate(Comparison, into = c("group1", "group2"), sep = " - ")
  dunn_plot <- dunn %>% 
    filter(P.adj<0.05)

  p <- ggplot(df_sub, aes(x = ILAE, y = avg_percent_reg, fill = ILAE))+
    geom_boxplot(lwd = 0.2, colour = "black")+ #can be adjusted
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          #legend.position = "none"
    )+
    scale_fill_manual(values = color_ILAE)+
    geom_jitter(color="black", size=s, alpha=0.9, width = 0.15)+
    ggtitle(ttl) +
    ylab("Proportion")+
    coord_cartesian(ylim = y_limit, clip = "off") # can be changed with the next line
  
  if (nrow(dunn_plot)>0){
    max_df <- max(df_sub$avg_percent_reg)
    position_comparison = seq(from = max_df + 1, to = max_df + nrow(dunn_plot)*1, by = 1)
    p <- p +
      stat_pvalue_manual(dunn_plot, label = "significance", tip.length = 0.01, y.position = position_comparison)
  }
  print(p)
}

#region we want to plot
reg <- c("CA2_CA4","SL.SR.SLM", "GABA")

p <- boxplot_comparison(proportion_CAs, loc = reg, ttl = paste(reg,collapse = "_"), y_limit = c(0,100), s = 0.5)

# save path and creating folder
save_path <- file.path(fig_path,"figure2/supplemental/S2B_hipporegionproportions")
#dir.create(save_path, recursive = TRUE)
pdf(file.path(save_path,paste0("proportions_ILAE_", paste(reg,collapse = "_"), ".pdf")), width = 5, height = 4)
print(p)
dev.off()
