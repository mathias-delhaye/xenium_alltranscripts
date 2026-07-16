#' Highlight each cluster for every sample and define a common identity for these clusters throughout the samples
#' Check the UMAP after defining the regions to check for potential artifacts
#' Mathias Delhaye, July 2025

# Load libraries
library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)

# Set paths
path <- here()
data_path <- str_remove(path,"xenium_epilepsy_repo")
fig_path <- file.path(data_path,"figures_v5/")

# Load files
megadata <- readRDS(file.path(data_path,"MegaData_v5.rds"))
seulist_log30_v5 <- readRDS(file.path(data_path,"seulist_log30_v5.rds"))
seuint_log30_k23_v5 <- readRDS(file.path(data_path,"seuint_log30_k23_v5.rds"))

# Create a function to locate the position of specific clusters in each tissue section taking as input
# seuint = integrated seurat object
# orig = sample name (needs " ")
# cluster_nb = number of the cluster which will be highlighted
# k_cluster = number of clusters seuint has
cluster_location <- function(seuint, orig, cluster_nb, k_cluster){
  sub <- subset(seuint, subset = orig.ident == orig) # subset to keep only the spots of the sample wanted
  sub@meta.data <- sub@meta.data %>% mutate(batch = orig) # rename batch with the sample name for plotting wit spaplot
  cols_sub <- brewer.blues(n = k_cluster) # k_cluster number of blue (gradient)
  cols_sub[cluster_nb] = "#d62728FF" # replace the blue from the gradient of blues located at the number of the cluster which will be highlighter in red
  sub %>% SpaPlot(item = "cluster", batch = NULL, cols = cols_sub, point_size = 3.5, nrow.legend = 26, title_name = "") # function from PRECAST package to plot the cluster in a 2D plot (scatterplot)
}

# Run the cluster_location through all clusters for all samples
sample_names <- names(seulist_log30_v5)
for (smpl in sample_names){
  dir.create(file.path(fig_path,"highlight_clusters_position_v5/k_23",smpl), showWarnings = FALSE) #create a folder with subfolders for each sample
  path = file.path(fig_path,"highlight_clusters_position_v5/k_23",smpl)
  for (i in (1:23)){ #need to use the number of k_cluster contained by seuint
    png(file.path(path,paste0(smpl,"_cluster",as.character(i),"_k23_sans_ILAE3.png")))
    print(cluster_location(seuint_log30_k23_v5,smpl,i, 23)) #needs to be explicitely printed
    dev.off()
  } 
}

# The user looks at each image for each sample to give an identity to all the clusters. 
# Usually it requires some back and for between 03_PRECAST_to_seuInt.R and here to test multiple number of k_clusters and choose the one that spatially cluster the best the dataset.

# Here k = 23 was chosen

# Add the regions and subregions information to the metadata of seuInt
seuint_log30_k23_v5@meta.data <- seuint_log30_k23_v5@meta.data %>%
  mutate(regions = case_when(
           cluster == 20 | cluster == 15 ~ "DG",
           cluster == 2 | cluster == 7 | cluster == 17 | cluster == 22 ~ "CA3/4_PL",
           cluster == 5 | cluster == 10 | cluster == 18 ~ "CA1_PL",
           cluster == 8 | cluster == 9 | cluster == 12 | cluster == 13 | cluster == 21 ~ "CA1/3/4_neuropil",
           cluster == 1 ~ "subiculum",
           cluster == 3 | cluster == 4 | cluster == 6 ~ "alveus",
           cluster == 14 ~ "pia",
           cluster == 11 | cluster == 23 ~ "bloodvessels",
           cluster == 16 | cluster == 19 ~ "void"),
         subregions = case_when(
           cluster == 20 ~ "DG_ML",
           cluster == 15 ~ "DG_GC",
           cluster == 2 ~ "CA3/4_PL_1",
           cluster == 7 ~ "CA3/4_PL_2",
           cluster == 17 ~ "CA3/4_PL_3",
           cluster == 22 ~ "CA3/4_PL_4",
           cluster == 5 ~ "CA1_PL1",
           cluster == 10 ~ "CA1_PL2",
           cluster == 18 ~ "CA1_PL3",
           cluster == 8 ~ "CA1/3/4_neuropil1",
           cluster == 9 ~ "CA1/3/4_neuropil2",
           cluster == 12 ~ "CA1/3/4_neuropil3",
           cluster == 13 ~ "CA1_neuropil1",
           cluster == 21 ~ "CA1_neuropil2",
           cluster == 1 ~ "subiculum",
           cluster == 3 ~ "alveus1",
           cluster == 4 ~ "alveus2",
           cluster == 6 ~ "alveus3",
           cluster == 14 ~ "pia",
           cluster == 11 | cluster == 23 ~ "bloodvessels",
           cluster == 16 | cluster == 19 ~ "void"))

# Define colors for the regions and subregions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2-CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

# Plot UMAP of the integrated seurat object coloring by ILAE_score & samples to check for any obvious artifact

# Generate a table with ILAE type for each sample
ILAE_table <- megadata@meta.data %>% 
  group_by(orig.ident,ILAE_score) %>% 
  tally() %>% 
  select(orig.ident,ILAE_score) %>% 
  mutate(ILAE_score = as.factor(ILAE_score))

# Add UMAP to seuint object
seuint_umap <- AddUMAP(seuint)
seuint_umap@meta.data <- left_join(seuint_umap@meta.data, ILAE_table, by = "orig.ident")

# Save umap by coloring by orig.ident
png(file.path(fig_path,"umap/k23",paste0("umap_k23_v5_cases.png")))
dimPlot(seuint_umap, item = "orig.ident")
dev.off()

# Save umap by coloring by subregions
png(file.path(fig_path,"umap/k23",paste0("umap_k23_v5_regions.png")))
dimPlot(seuint_umap, item = "regions", cols = col_regions)
dev.off()

# Save umap by coloring by ILAE_score
png(file.path(fig_path,"umap/k23",paste0("umap_k23_v5_ILAE.png")))
dimPlot(seuint_umap, item = "ILAE_score")
dev.off()