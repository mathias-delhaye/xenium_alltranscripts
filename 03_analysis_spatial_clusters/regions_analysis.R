## Quality control plotting for seumerged_log30 object

## Loading packages

library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)

## setting paths
path = here()
data_path = file.path(str_remove(path,"xenium_epilepsy_repo"))
fig_path = file.path(data_path,"figures_v6")

## Load files
megadata <- read_rds(file.path(data_path,"MegaData_v6.rds"))
metadata_megadata <- read_rds(file.path(data_path,"metadata_md_v6.rds"))
seulist <- readRDS(file.path(data_path,"seulist_v6.rds"))
seulist_log30 <- readRDS(file.path(data_path,"seulist_log30_v6.rds"))
PRECASTObj_log30 <- readRDS(file.path(data_path,"PRECASTObj_log30_v6.rds"))
PRECASTObj_log30_k23 <- readRDS(file.path(data_path,"PRECASTObj_log30_k23_v6.rds"))
seuint <- readRDS(file.path(data_path,"seuint_log30_k23_v6.rds"))
seumerged <- readRDS(file.path(data_path,"seumerged_log30_v6.rds"))

# sample names and genes list
sample_names <- names(megadata@images)
gene_list <- Features(megadata)

# defining gradient of colors for the regions and subregions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2-CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

#plotting regions proportions

# function calculating the proportions and the spot counts based on the feature selected
feature_proportion <- function(seu,feature = "orig.ident", loc = "regions"){
  seu@meta.data %>%
    filter(regions != "void") %>% # removing the void regions
    group_by_at(c(feature,"orig.ident", loc)) %>%
    tally() %>% 
    mutate(percent_feature = 100*n/sum(n)) %>% 
    ungroup() %>% 
    group_by_at(c(feature, loc)) %>% 
    summarise(avg_percent_reg = mean(percent_feature)) %>%
    ungroup()
}

feature_n <- function(seu,feature = "smpl_ILAE", loc = "regions"){
  seu@meta.data %>%
    filter(regions != "void") %>% # removing the void regions
    group_by_at(c(feature,"smpl_ILAE", loc)) %>%
    tally()
}

plot_proportions <- function(df_proportion, feature = "orig.ident", loc = "regions", y = "avg_percent_reg", cols, y_lab = "percent/regions"){
  # Convert column names from character to symbols for tidy evaluation
  feature_sym <- sym(feature)
  loc_sym <- sym(loc)
  y_sym <- sym(y)
  
  df_proportion %>% 
    ggplot(aes(x = !!feature_sym, y = !!y_sym, fill = !!loc_sym))+
    geom_bar(stat = "identity", width = .95) +
    ylab("Proportion") +
    scale_fill_manual(name = loc, values = cols) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
}

seumerged_ILAE_reg_proportion <- feature_proportion(seumerged, feature = "ILAE_score", loc = "regions")
png(file.path(fig_path,"regions_proportions_ILAE.png"))
p <- plot_proportions(seumerged_ILAE_reg_proportion, feature = "ILAE_score", loc = "regions", cols = col_regions)
print(p)
dev.off()

seumerged_ILAE_reg_md_proportion <- feature_proportion(seumerged, feature = "ILAE_score", loc = "regions_md")
png(file.path(fig_path,"regions_md_proportions_ILAE.png"))
p <- plot_proportions(seumerged_ILAE_reg_md_proportion, feature = "ILAE_score", loc = "regions_md", cols = col_regions)
print(p)
dev.off()

seumerged_case_reg_proportion <- feature_proportion(seumerged, feature = "smpl_ILAE", loc = "regions")
png(file.path(fig_path,"regions_proportions_cases.png"), width = 960, height = 480)
p <- plot_proportions(seumerged_case_reg_proportion, feature = "smpl_ILAE", loc = "regions", cols = col_regions)
print(p)
dev.off()

seumerged_case_reg_md_proportion <- feature_proportion(seumerged, feature = "smpl_ILAE", loc = "regions_md")
pdf(file.path(fig_path, "figure2_panels", "panelB","regions_md_proportions_cases.pdf"), width = 7, height = 6)
p <- plot_proportions(seumerged_case_reg_md_proportion, feature = "smpl_ILAE", loc = "regions_md", cols = col_regions)
print(p)
dev.off()

seumerged_cases_n_reg <- feature_n(seumerged)
png(file.path(fig_path,"regions_proportions_cases_nspots.png"), width = 960, height = 480)
p <- plot_proportions(seumerged_cases_n_reg, feature = "smpl_ILAE", y = "n",cols = col_regions, y_lab = "n_spot")
print(p)
dev.off()

seumerged@meta.data %>%
  filter(regions != "void") %>% # removing the void regions
  ggplot(aes(x = cluster_precast, y = nCount_RNA))+
  geom_violin()

seumerged@meta.data %>%
  filter(regions != "void") %>% # removing the void regions
  ggplot(aes(x = regions, y = nCount_RNA))+
  geom_violin()

seumerged@meta.data %>% 
  filter(regions != "void") %>%
  group_by(cluster_precast) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = cluster_precast, y = avg_ncount))+
  geom_bar(stat = "identity")

seumerged@meta.data %>% 
  filter(regions != "void") %>%
  group_by(regions) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = regions, y = avg_ncount))+
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

seumerged@meta.data %>% 
  filter(regions != "void") %>%
  group_by(smpl_ILAE) %>% 
  summarise(tot_ncount = sum(nCount_RNA)) %>% 
  ggplot(aes(x = smpl_ILAE, y = tot_ncount))+
  geom_bar(stat = "identity")

seumerged@meta.data %>% 
  filter(regions != "void") %>%
  group_by(smpl_ILAE) %>% 
  summarise(avg_ncount = mean(nCount_RNA)) %>% 
  ggplot(aes(x = smpl_ILAE, y = avg_ncount))+
  geom_bar(stat = "identity")

png(file.path(fig_path,"avg_transcriptcount_cluster.png"))
print(p)
dev.off()