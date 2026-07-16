#' Integrate the PRECAST object.
#' Save the integrated Seurat object (seuint)
#' Mathias Delhaye, July 2025

# Load libraries
library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)

# Set paths
path <- here()
data_path <- str_remove(path,"xenium_epilepsy_repo")

## Load PRECASTObj
PRECASTObj_log30_k23_v5 <- read_rds(file.path(data_path,"PRECASTObj_log30_k23.rds"))

## Selecting model
seuint_log30_k23_v5 <- PRECASTObj_log30_k23_v5 %>% 
  SelectModel() %>% 
  IntegrateSpaData(species = "unknown")

# exporting seuint_log30_k17 into RDS object
saveRDS(seuint_log30_k23_v5,file.path(data_path,"seuint_log30_k23_v5.rds"))