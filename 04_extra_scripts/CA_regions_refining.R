## Run PRECAST on the CA regions only to try to refine and get CA2

library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)
library(grid)

## setting path
data_path <- str_remove(here(),"/xenium_epilepsy_repo")

## Loading seulist filtered
seulist_log30_v5 <- read_rds(file.path(data_path,"seulist_log30_v5.rds"))

## Initializing the list
seulist_log30_CA <- seulist_log30_v5

## getting metadata from seuint to keep only regions column with id of spots
metadata_seuint <- seuint_log30_k23_v5@meta.data %>% 
  select(regions) %>% 
  rownames_to_column(var = "id")

## defining the regions to use. Initially only PL regions but not enough, so added neuropil of CA
CA_regions <- c("CA1_PL","CA3/4_PL","CA1/3/4_neuropil")

## subsetting the list to keep only spots in CA regions
for (im in names(seulist_log30_CA)){
  seulist_log30_CAPL[[im]]@meta.data <- seulist_log30_CA[[im]]@meta.data %>% 
    rownames_to_column(var = "id") %>% 
    left_join(metadata_seuint, by = "id") %>% 
    column_to_rownames(var = "id")
  seulist_log30_CA[[im]] <- subset(seulist_log30_CA[[im]], subset = regions %in% PL_regions)
}
rm(metadata_seuint)

## Create PRECAST object
preobj_CA <- CreatePRECASTObject(seuList = seulist_log30_CA)
saveRDS(preobj_CA,file.path(data_path,"CA_refinement/preobj_CA.rds"))

## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
PRECASTObj_CA <- AddAdjList(preobj_CA, platform = "Other_SRT")
saveRDS(PRECASTObj_CA,file.path(data_path,"CA_refinement/PRECASTObj_CA.rds"))

## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj_CA <- AddParSetting(PRECASTObj_CA, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

## Try K = 15 and 20, initially there was 12 different regions defined
PRECASTObj_CA_k15 <- PRECAST(PRECASTObj_CA, K = 15)
saveRDS(PRECASTObj_CA_k15,file.path(data_path,"CA_refinement/PRECASTObj_CA_k15.rds"))

PRECASTObj_CA_k20 <- PRECAST(PRECASTObj_CA, K = 20)
saveRDS(PRECASTObj_CA_k20,file.path(data_path,"CA_refinement/PRECASTObj_CA_k20.rds"))

## Selecting model
seuint_CA_k15 <- PRECASTObj_CA_k15 %>% 
  SelectModel() %>% 
  IntegrateSpaData(species = "unknown")
saveRDS(seuint_CA_k15,file.path(data_path,"CA_refinement/seuint_CA_k15.rds"))

seuint_CA_k20 <- PRECASTObj_CA_k20 %>% 
  SelectModel() %>% 
  IntegrateSpaData(species = "unknown")
saveRDS(seuint_CA_k20,file.path(data_path,"CA_refinement/seuint_CA_k20.rds"))

seuint_CA_k15 <- read_rds(file.path(data_path,"CA_refinement/seuint_CA_k15.rds"))
seuint_CA_k20 <- read_rds(file.path(data_path,"CA_refinement/seuint_CA_k20.rds"))

# locate the position of specific clusters in each tissue section
cluster_location <- function(seuint, orig, cluster_nb, k_cluster){
  sub <- subset(seuint, subset = orig.ident == orig)
  sub@meta.data <- sub@meta.data %>% mutate(batch = orig)
  cols_sub <- brewer.blues(n = k_cluster)
  cols_sub[cluster_nb] = "#d62728FF"
  sub %>% SpaPlot(item = "cluster", batch = NULL, cols = cols_sub, point_size = 3.5, nrow.legend = 26, title_name = "")
}

# running through all clusters for all samples
fig_path <- file.path(data_path,"figures_v5")

for (smpl in names(seulist_log30_CA)){
  dir.create(file.path(fig_path,"highlight_CA_clusters_position_v5/k_15",smpl), showWarnings = FALSE)
  path = file.path(fig_path,"highlight_CA_clusters_position_v5/k_15",smpl)
  for (i in (1:15)){
    png(file.path(path,paste0(smpl,"_cluster",as.character(i),"_k15.png")))
    print(cluster_location(seuint_CA_k15,smpl,i, 15)) #needs to be explicitely printed
    dev.off()
  }
  dir.create(file.path(fig_path,"highlight_CA_clusters_position_v5/k_20",smpl), showWarnings = FALSE)
  path = file.path(fig_path,"highlight_CA_clusters_position_v5/k_20",smpl)
  for (i in (1:20)){
    png(file.path(path,paste0(smpl,"_cluster",as.character(i),"_k20.png")))
    print(cluster_location(seuint_CA_k20,smpl,i, 20)) #needs to be explicitely printed
    dev.off()
  }
}