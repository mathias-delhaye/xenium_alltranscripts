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
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)

##################################################
## Section 2: Environment
##################################################

# Set paths
choose.files()
data_path <- choose.dir()
rds_path <- file.path(data_path, "01_rds_files")
fig_path <- file.path(data_path, "02_plots")

# Load files
megadata <- read_rds(file.path(rds_path,"01_MegaData_v8.rds"))
seumerged <- read_rds(file.path(rds_path,"07_seumerged_v8.rds"))

# Color scheme for hippocampal regions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

##################################################
## Section 3: Marker genes for each region
##################################################

# Join layers for seurat functions
seumerged_joined <- JoinLayers(seumerged)

# Remove void (and NAs)
seumerged_joined <- subset(x = seumerged_joined, subset = region != "void")

# Find all markers for precast_clusters, regions and subregions
markers_regions <- FindAllMarkers(seumerged_joined, group.by = "region", only.pos = TRUE)

# Keeping the top 10 for each, with log2 fold change >1
top10_regions <- markers_regions %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC>1) %>% 
  slice_head(n = 10) %>%
  ungroup()

# plot the marker genes in heatmaps
DoHeatmap(seumerged_joined, features = top10_regions$gene, group.by = "region", label = F, group.colors = col_regions)

# table with marker genes and regions

top10_regions_markers_only <- top10_regions %>% 
  select(cluster, gene) %>%
  rename(regions = cluster)

for (region in unique(top10_regions_markers_only$regions)){
  dir.create(file.path(fig_path,"04_marker_genes_regions",region), recursive = TRUE)
  markers <- top10_regions %>%
    filter(cluster == region) %>%
    pull(gene) #extract values of gene column
  for (g in markers){
    png(file.path(fig_path,"04_marker_genes_regions",region,paste0(g,"_vln_regions.png")))
    print(VlnPlot(seumerged_joined,features = g, group.by = 'region', pt.size = 0, cols = col_regions)) #needs to be explicitely printed
    dev.off()
  }
}
rm(region)

marker.features<- c('GRM2','NTNG2','KIT','WIF1','GRM4','GAD2','RELN', 'MAG','TGFBI','TP73')

#i nstead of defining before features can also be added as list in function
# and group by whatever metadata you want
p <- DotPlot(seumerged_joined, features = marker.features, group.by = 'region_md') +
  scale_color_gradientn(colors = c("cyan", "white", "red")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.position = "top",
        #legend.box = "vertical",
        legend.justification = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))

plotpath <- file.path(fig_path,"04_marker_genes_regions","dotplots_markergenes_regions_v8.png")
ggsave(filename = plotpath,
       plot = p,
       dpi = 300,
       width = 6,
       height = 6
       )
