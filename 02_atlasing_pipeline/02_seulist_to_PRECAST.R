#' Run PRECAST on Seurat list filtered and log-normalized.
#' Save the PRECAST objects with and without the adjacent matrix, and PRECAST object with the spatial clusters.
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

# Load Seulist
seulist_log30_v5 <- read_rds(file.path(data_path,"seulist_log30_v5.rds"))

# Create PRECAST object
preobj_v5 <- CreatePRECASTObject(seuList = seulist_log30_v5)

# Save preobj
saveRDS(preobj_v5,file.path(data_path,"preobj_v5.rds"))

# Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
PRECASTObj_log30_v5 <- AddAdjList(preobj_v5, platform = "Other_SRT")

# Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj_log30_v5 <- AddParSetting(PRECASTObj_log30_v5, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

# Implement chosen number of desired clusters. In our case 23 works the best
PRECASTObj_log30_k23_v5 <- PRECAST(PRECASTObj_log30_v5, K = 23)

# Saving into an RDS object the PRECASTobj and PRECASTobj k-cluster = 13
saveRDS(PRECASTObj_log30_v5,file.path(data_path,"PRECASTObj_log30_v5.rds"))
saveRDS(PRECASTObj_log30_k23_v5,file.path(data_path,"PRECASTObj_log30_k17_v5.rds"))