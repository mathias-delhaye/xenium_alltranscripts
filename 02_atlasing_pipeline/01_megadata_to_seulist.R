#' Build a Seurat list containing all samples binned with the all the transcripts detected by xenium contained in each bin, which is filtered and log-normalized.
#' Save the Seurat lists with the raw counts and the filtered and log-normalized counts.
#' Mathias Delhaye, July 2025

# Load libraries

library(tidyverse)
library(Seurat)
library(ggplot2)
library(here)
library(ggrepel)

# Set paths

path <- here()
data_path <- str_remove(path,"xenium_epilepsy_repo")

# Load megadata
megadata = read_rds(file.path(data_path,"MegaData_v5.rds"))
# name of all samples in megadata in order assigned in megadata
sample_names <- names(megadata@images)
# all the genes analyzed with xenium (aka 266 panel + 100 epil specific)
gene_list <- Features(megadata)

## Building seulist with raw transcript count

#initiating seurat list
seulist <- list()

# set bin dimension arbitrarily - 150*150px works fine
bin_width <- 150
bin_height <- 150

for (im in sample_names){
  # field of view where all the genes are stored individually, containing the coordinates of all the transcripts
  fov <- megadata@images[[im]]@molecules$molecules
  
  # create a tibble to store the genes + their coordinates
  transcript_coord <- tibble(gene = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)
  
  for (n in (1:length(gene_list))){
    gene <- gene_list[n] #going through the list of genes
    gene_matrix <- fov[[gene]]@coords #selecting in the fov the gene
    temp_df <- tibble(gene = gene, x = gene_matrix[,1], y = gene_matrix[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
    transcript_coord <- rbind(transcript_coord,temp_df) # joining the temp df with the tibble containing all genes
    rm(temp_df)
  }
  
  # Assign each transcript to a spot
  transcript_coord$bin_x <- ceiling(transcript_coord$x / bin_width) # x coordinate / width of bin
  transcript_coord$bin_y <- ceiling(transcript_coord$y / bin_height) # y coordinate / height of bin
  
  # Find max number of bin per column and row
  number_bin_x <- ceiling(max(transcript_coord$x)/bin_width)
  number_bin_y <- ceiling(max(transcript_coord$y)/bin_height)
  
  # Gives the order of bins
  transcript_coord$bin_number <- (transcript_coord$bin_y-1)*number_bin_x+transcript_coord$bin_x # (bin_y coordinate - 1 aka number of full lines) * number of bins in X  + coordinate bin_x for the number of bins to add to the full lines
  
  # counting the number of transcripts for each bin and each gene
  transcript_genexbin <- transcript_coord %>% 
    group_by(bin_number, gene) %>% 
    tally() %>% 
    ungroup()
  
  # make the tibble wider to have nb of row = nb of genes and nb of column = nb of bins
  transcript_genexbin <- transcript_genexbin %>% 
    pivot_wider(names_from = bin_number, values_from = n) %>% 
    arrange(gene) %>% 
    column_to_rownames(var = "gene") %>% 
    rename_with(~ paste0(im,"_", .)) #change the name of the bin to add the sample name _ and not just having numbers
  
  # find the bins that we kept through previous operation, which are bins with at least 1 transcript
  bin_w_transcript <- colnames(transcript_genexbin)
  # find the missing bins and add empty columns with these missing bins
  missing_bin = setdiff(paste0(im,"_",as.character(1:(number_bin_x*number_bin_y))),bin_w_transcript)
  transcript_genexbin[,missing_bin] <- NA
  
  #convert tibble into matrix
  transcript_genexbin_mat <- as.matrix(transcript_genexbin)
  
  # replace NA by 0
  transcript_genexbin_mat[is.na(transcript_genexbin_mat)] <- 0
  
  # re-order the column to have them in the order of bins
  transcript_genexbin_mat <- transcript_genexbin_mat[,paste0(im,"_",as.character(1:(number_bin_x*number_bin_y)))]
  
  # Get the coordinates for each bin
  bin_grid <- expand_grid(row = 1:number_bin_y, col = 1:number_bin_x) %>% #create all combination possible with the value given (here being the number of bins in X and Y), with the second argument being the one used first to increase value (1-1, 1-2, 1-3, ... 1-n, 2-1, 2-2, ... 2-n, etc...)
    select(col,row) %>% # re-order the columns to have x first then y
    as.data.frame() # convert into a df
  # name each row with the corresponding bin
  row.names(bin_grid) <- paste0(im,"_",as.character(1:(number_bin_x*number_bin_y)))
  # adding the position windows for each bin in X & Y so that we can place segmented area in each bin based on the centroids position
  bin_grid$x_window <- bin_grid$col*bin_width
  bin_grid$y_window <- bin_grid$row*bin_height
  
  seu <- CreateSeuratObject(
    counts = transcript_genexbin_mat,
    assay = "RNA",
    meta.data = bin_grid
  )
  
  # change orig.ident to set it to sample name - needs to use level as it is considered as a factor
  levels(seu@meta.data$orig.ident) <- im
  
  seulist[[im]] <- seu
}

#saving seurat list
saveRDS(seulist,file.path(data_path,"seulist_v5.rds"))

## Building seulist_log30 filter with >30 counts/spot + log normalization

# Initializing the list
seulist_log30 <- list()

# Going through the list
for (im in sample_names){
  seulist_log30[[im]] <- seulist[[im]] %>% 
    subset(subset = nCount_RNA>30) %>% #filter out bins with <30 transcripts
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData()
}

saveRDS(seulist_log30,file.path(data_path,"seulist_log30_v5.rds"))