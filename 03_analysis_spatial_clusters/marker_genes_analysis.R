## marker genes analysis for the different clusters/regions

## Loading packages

library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)

## setting path
data_path <- str_remove(here(),"/xenium_epilepsy_repo")
fig_path <- file.path(data_path, "figures_v6")

## Loading files
megadata <- readRDS(file.path(data_path,"MegaData_v6.rds"))
seumerged <- readRDS(file.path(data_path,"seumerged_log30_v6.rds"))

# defining gradient of colors for the regions in subregions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2-CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

# Joining layers for seurat functions
seumerged_joined <- JoinLayers(seumerged)

# removing void (and NAs)
seumerged_joined <- subset(x = seumerged_joined, subset = regions != "void")

# Find all markers for precast_clusters, regions and subregions
markers_regions <- FindAllMarkers(seumerged_joined, group.by = "regions", only.pos = TRUE)

# Keeping the top 10 for each, with log2 fold change >1
top10_regions <- markers_regions %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC>1) %>% 
  slice_head(n = 10) %>%
  ungroup()

# plot the marker genes in heatmaps
DoHeatmap(seumerged_joined, features = top10_regions$gene, group.by = "regions", label = F, group.colors = col_regions)

# table with marker genes and regions

top10_regions_markers_only <- top10_regions %>% 
  select(cluster, gene) %>%
  rename(regions = cluster)

for (region in unique(top10_regions_markers_only$regions)){
  dir.create(file.path(fig_path,"marker_genes_regions",region), recursive = TRUE)
  markers <- top10_regions %>%
    filter(cluster == region) %>%
    pull(gene) #extract values of gene column
  for (g in markers){
    png(file.path(fig_path,"marker_genes_regions",region,paste0(g,"_vln_regions.png")))
    print(VlnPlot(seumerged_joined,features = g, group.by = 'regions', pt.size = 0, cols = col_regions)) #needs to be explicitely printed
    dev.off()
  }
}
rm(region)

for (gene in marker.features){
  png(file.path(fig_path,"marker_genes_regions",,paste0(gene,"_vln_regions.png")))
  print(VlnPlot(seumerged_joined,features = gene, group.by = 'regions', pt.size = 0, cols = col_regions)) #needs to be explicitely printed
  dev.off()
}
rm(gene)

seumerged_joined$regions_md <- factor(seumerged_joined$regions_md, levels = c("GCL","ML","CA2_CA4","CA1","SUB","GABA","SL.SR.SLM","WM","vascular","pia"))
marker.features<- c('GRM2','SEMA5A','KIT','WIF1','GRM4','GAD2','RELN', 'MAG','TGFBI','TP73')

#instead of defining before features can also be added as list in function
# and group by whatever metadata you want
p <- DotPlot(seumerged_joined, features = marker.features, group.by = 'regions_md') +
  scale_color_gradientn(colors = c("cyan", "white", "red")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.position = "top",
        #legend.box = "vertical",
        legend.justification = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))

pdf_file <- file.path(fig_path,"figure2_panels","panelC",paste0("panelC_dotplots_markergenes.pdf"))
pdf(pdf_file, height = 4.5, width = 4)
print(p)
dev.off()