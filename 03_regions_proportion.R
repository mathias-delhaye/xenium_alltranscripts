##################################################
## Script purpose: Run marker genes analysis for the different clusters/regions
## Date: 2026-07-28
## Update date: 2026-07-28
## Author: Mathias Delhaye
##################################################

##################################################
## Section 1: Libraries
##################################################

library(tidyverse)
library(here)
library(future)
library(pals)

##################################################
## Section 2: Environment
##################################################

# Set paths
choose.files()
data_path <- choose.dir()
rds_path <- file.path(data_path, "01_rds_files")
fig_path <- file.path(data_path, "02_plots","05_region_proportions")
dir.create(fig_path,recursive = T, showWarnings = F)

# Load files
seumerged <- read_rds(file.path(rds_path,"07_seumerged_v8.rds"))

# Color scheme for hippocampal regions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))
color_ILAE <- setNames(c("#C67BFF","#F8766C","#7BAD00","#00C0C5","#E69F00"),
                       c("P", "ILAE0","ILAE1","ILAE2","ILAE3"))

##################################################
## Section 3: helper functions 
##################################################

# function calculating the proportions and the spot counts based on the feature selected
feature_proportion <- function(seu,feature = "orig.ident", loc = "region"){
  seu@meta.data %>%
    filter(region != "void") %>% # removing the void regions
    group_by_at(c(feature, loc)) %>%
    tally() %>% 
    ungroup() %>%
    group_by_at(c(feature)) %>% 
    mutate(percent_feature = 100*n/sum(n)) %>% 
    ungroup()
}

# function to plot the proportions stacked for the feature selected
plot_proportions <- function(df_proportion, feature = "orig.ident", loc = "region", y = "percent_feature", cols, y_lab = "percent/region"){
  # Convert column names from character to symbols for tidy evaluation
  feature_sym <- sym(feature)
  loc_sym <- sym(loc)
  y_sym <- sym(y)
  
  df_proportion %>% 
    ggplot(aes(x = !!feature_sym, y = !!y_sym, fill = !!loc_sym))+
    geom_bar(stat = "identity", width = .95) +
    ylab("Proportion") +
    scale_fill_manual(name = loc, values = cols) +
    theme_minimal()+
    theme(axis.text.x = element_text(size = 8,
                                     colour = "black",
                                     angle = 45),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10,
                                     colour = "black"),
          axis.title.y = element_text(size = 12,
                                      colour = "black"),
          panel.grid = element_blank())
}

##################################################
## Section 4: Calculate and plot proportions for each case and ILAE 
##################################################

seumerged_ILAE_reg_proportion <- feature_proportion(seumerged, feature = "ILAE_score", loc = "region_md")

p <- plot_proportions(seumerged_ILAE_reg_proportion, feature = "ILAE_score", loc = "region_md", cols = col_regions)
p
ggsave(
  filename = file.path(fig_path, "02_regionmdprop_perILAE.png"),
  plot = p,
  dpi = 300,
  width = 6,
  height = 6
)

seumerged_case_reg_proportion <- feature_proportion(seumerged, feature = "smpl_ILAE", loc = "region")
p <- plot_proportions(seumerged_case_reg_proportion, feature = "smpl_ILAE", loc = "region", cols = col_regions)
p
ggsave(
  filename = file.path(fig_path, "03_regionprop_percase.png"),
  plot = p,
  dpi = 300,
  width = 9,
  height = 6
)

##################################################
## Section 5: plot proportions CAs regions
##################################################

# removing void (and NAs)
proportion_CAs <- subset(x = seumerged, subset = region_md %in% c("CA2_CA4","CA1","GABA","SL.SR.SLM")) %>% 
  feature_proportion(feature = "smpl_ILAE", loc = "region_md") %>% 
  separate(smpl_ILAE, into = c("smpl", "ILAE"),sep = "_") %>%
  mutate(ILAE = case_when(
    ILAE == "0" ~ "ILAE0",
    ILAE == "1" ~ "ILAE1",
    ILAE == "2" ~ "ILAE2",
    ILAE == "3" ~ "ILAE3",
    TRUE ~ ILAE
  ),
  ILAE = factor(ILAE, c( "P", "ILAE0", "ILAE1", "ILAE2", "ILAE3")))

CA_regions <- c("CA2_CA4","CA1","GABA","SL.SR.SLM")
for (reg in CA_regions){
  # isolate for the specific region
  df_sub <- filter(proportion_CAs, region_md == reg)
  dunn <- dunnTest(percent_feature ~ ILAE, data = df_sub, method = "bh")
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
  p <- df_sub %>% 
    ggplot(aes(x = ILAE, y = percent_feature, fill = ILAE))+
    geom_boxplot(lwd = 0.2, colour = "black",,outlier.shape = NA)+ #can be adjusted
    theme_minimal()+
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black",
                                     size = 9),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_blank(),
          title = element_text(colour = "black",
                               size = 9),
    )+
    scale_fill_manual(values = color_ILAE)+
    geom_jitter(color="black", size=0, alpha=0.9, width = 0.15)+
    ggtitle(paste0(reg)) +
    ylab("Proportion")+
    coord_cartesian(ylim = c(0,60), clip = "off") # can be changed with the next line
  
  if (nrow(dunn_plot)>0){
    max_df <- max(df_sub$percent_feature)
    position_comparison = seq(from = max_df + 1, to = max_df + nrow(dunn_plot)*3, by = 3)
    p <- p +
      stat_pvalue_manual(dunn_plot, label = "significance", tip.length = 0.01, y.position = position_comparison)
  }
  ggsave(
    filename = file.path(fig_path,paste0("05_proportion_",reg,".png")),
    plot = p,
    width = 2,
    height = 4,
    dpi = 300
  )
}

##################################################
## Section 6: additional plot for QC (number of counts per case/region)
##################################################

seumerged@meta.data %>%
  filter(region != "void") %>% # removing the void regions
  ggplot(aes(x = cluster_precast, y = nCount_RNA))+
  geom_violin()

seumerged@meta.data %>%
  filter(region != "void") %>% # removing the void regions
  ggplot(aes(x = region, y = nCount_RNA))+
  geom_violin()

seumerged@meta.data %>% 
  filter(region != "void") %>%
  group_by(cluster_precast) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = cluster_precast, y = avg_ncount))+
  geom_bar(stat = "identity")

seumerged@meta.data %>% 
  filter(region != "void") %>%
  group_by(region) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = region, y = avg_ncount))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

seumerged@meta.data %>% 
  filter(region != "void") %>%
  group_by(smpl_ILAE) %>% 
  summarise(tot_ncount = sum(nCount_RNA)) %>% 
  ggplot(aes(x = smpl_ILAE, y = tot_ncount))+
  geom_bar(stat = "identity")

seumerged@meta.data %>% 
  filter(region != "void") %>%
  group_by(smpl_ILAE) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = smpl_ILAE, y = avg_ncount))+
  geom_bar(stat = "identity")