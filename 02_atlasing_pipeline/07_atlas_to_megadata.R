#' Add regions, subregions and regions_md to the metadata of megadata.
#' Mathias Delhaye, July 2025

# Load libraries
library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)

# Set paths
data_path <- str_remove(here(),"/xenium_epilepsy_repo")
fig_path <- file.path(data_path, "figures_v5")

# Load files
megadata <- read_rds(file.path(data_path,"MegaData_v5.rds"))
seumerged_log30_v5 <- readRDS(file.path(data_path,"seumerged_log30_v5.rds"))

# create a tibble to store the genes + their coordinates
cell_coord <- tibble(cell_name = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

# Get the coordinates of each cells (centroid coordinate) for each sample
sample_names <- unique(seumerged_log30_v5$orig.ident)
for (smpl in sample_names){
  coord <- megadata@images[[smpl]]@boundaries[["centroids"]]@coords
  cell_names <- megadata@images[[smpl]]@boundaries[["centroids"]]@cells
  temp_df <- tibble(cell_names = cell_names, x = coord[,1], y = coord[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
  cell_coord <- rbind(cell_coord,temp_df) # joining the temp df with the tibble containing all genes
  rm(temp_df)
}

# Add column for the corresponding col and row coordinates by dividing the xy coordinates by 150(px) chosen to define the bins initially
# ceiling is used here as the bins are 150px by 150px and we want the location of each cell in the right col x row coord
cell_coord <- cell_coord %>% 
  mutate(col = ceiling(x/150), 
         row = ceiling(y/150))

# Extract the sample name for each cell using the metadata of megadata
metadata_megadata <- megadata@meta.data
cell_to_orig <- metadata_megadata %>% 
  rownames_to_column(var = "cell_names") %>% 
  select(cell_names,orig.ident) #keep only the cell name and the sample name
cell_coord <- left_join(cell_coord, cell_to_orig, by = "cell_names") # add the sample name to the cell_coord table
rm(cell_to_orig)


# Join cell_coord and metadata_seumerged by orig.ident, col and row to get the regions, subregions and regions md for each cell
# The sample name (orig.ident), col and row are needed to give each cell to the right region
metadata_seumerged <- seumerged_log30_v5@meta.data
cell_coord <- cell_coord %>%
  left_join(
    metadata_seumerged %>%
      select(orig.ident, col, row, regions, subregions, regions_md),
    by = c("orig.ident", "col", "row") 
  )

# Join the metadata of megadata and seumerged together using the cell_names
metadata_megadata <- metadata_megadata %>% 
  rownames_to_column(var = "cell_names") %>% 
  left_join(cell_coord %>% 
              select(cell_names, regions, subregions, regions_md),
            by = "cell_names") %>% 
  column_to_rownames(var = "cell_names")

saveRDS(metadata_megadata, file.path(data_path,"metadata_md.rds"))