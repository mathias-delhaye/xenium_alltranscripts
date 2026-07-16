#' Create a merged Seurat object using the integrated Seurat object and the filter and log-normalized Seurat list
#' Save scatter plots for each sample of all clusters obtained.
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

# Merge all the elements of the list containing the seurat objects for each case
seumerged_log30_v5 <- merge(seulist_log30_v5[[1]], seulist_log30_v5[-1])

# Add the ILAE score to the metadata of seumerged
metadata_seuint <- seuint_log30_k23_v5@meta.data %>% 
  select(cluster, regions, subregions) %>% 
  rownames_to_column(var = "id")
metadata_seumerged <- seumerged_log30_v5@meta.data %>% 
  rownames_to_column(var = "id") %>% 
  left_join(ILAE_table, by = "orig.ident") %>% 
  left_join(metadata_seuint, by = "id") %>% 
  rename(cluster_precast = cluster) %>% 
  column_to_rownames(var = "id") 
seumerged_log30_v5@meta.data <- metadata_seumerged
rm(metadata_seuint, metadata_seumerged)

# Save seumerged with the regions and subregions
saveRDS(seumerged_log30_v5,file.path(data_path,"seumerged_log30_v5.rds"))

# Plot for each sample the different regions and subregions in a scatterplot with the coordinates of each spot as X and Y to check the regions defined above.

# Save the plots

# Regions
for (smpl in sample_names){
  png(file.path(fig_path,"cluster_spa_location/k_23/regions",paste0(smpl,"_k23_v5_regions.png")))
  
  df <- seumerged_log30_v5@meta.data %>% 
    filter(orig.ident == smpl)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)
  
  p <-  df %>% 
    ggplot(aes(x = col, y = row, color = regions))+
    geom_point(size = 2.75)+
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
    ) +
    labs(title = smpl,
         color = "Region") + 
    annotation_custom(
      grob = rectGrob(
        gp = gpar(col = "black", fill = NA, lwd = 1)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  print(p) #needs to be explicitely printed
  dev.off()
}

#subregions
for (i in sample_names){
  png(file.path(fig_path,"cluster_spa_location/k_23/subregions",paste0(i,"_k23_v5_subregions.png")))
  
  df <- seumerged_log30_v5@meta.data %>% 
    filter(orig.ident == i)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)
  
  p <-  df %>% 
    ggplot(aes(x = col, y = row, colour = subregions))+
    geom_point(size = 2.75)+
    scale_color_manual(values = cols_subregions) +
    scale_y_reverse() +
    coord_fixed() +  # keep aspect ratio equal
    theme_minimal() +
    theme(
      axis.text = element_blank(),         # Remove axis numbers
      axis.ticks = element_blank(),        # Remove axis ticks
      axis.title = element_blank(),        # Remove axis titles
      panel.grid = element_blank(),        # Remove grid lines
      panel.border = element_blank(),  # turn off default border
    ) +
    labs(title = i,
         color = "Subregion") + 
    annotation_custom(
      grob = rectGrob(
        gp = gpar(col = "black", fill = NA, lwd = 1)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  print(p) #needs to be explicitely printed
  dev.off()
}