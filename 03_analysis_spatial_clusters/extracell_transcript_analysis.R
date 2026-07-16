
seulist <- readRDS(file.path(data_path,"seulist_v6.rds"))

incell_transcript_coord <- tibble(gene = character(), cell = character(), x = numeric(), stringsAsFactors = FALSE)

for (n in 1:length(sample_names)){
  xen_count <- megadata@assays$Xenium[n]
  temp_df <- as.data.frame(summary(xen_count))
  cell_names <- xen_count@Dimnames[[2]]
  temp_df$cell <- cell_names[temp_df$j]
  temp_df$gene <- gene_list[temp_df$i]
  temp_df <- temp_df[, c("gene", "cell", "x")]
  incell_transcript_coord <- rbind(incell_transcript_coord, temp_df)
}

metadata_seumerged <- seumerged@meta.data

transcript_per_bin_sansna$x[is.na(transcript_per_bin_sansna$x)] <- 0

df_bin_number <- metadata_seumerged %>% 
  rownames_to_column(var = "bin_number") %>% 
  select(orig.ident, col, row, bin_number)

rm(metadata_seumerged)

# create a tibble to store the genes + their coordinates
cell_coord <- tibble(cell = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

for (smpl in sample_names){
  coord <- megadata@images[[smpl]]@boundaries[["centroids"]]@coords
  cell_names <- megadata@images[[smpl]]@boundaries[["centroids"]]@cells
  temp_df <- tibble(cell = cell_names, x = coord[,1], y = coord[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
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
  rownames_to_column(var = "cell") %>% 
  select(cell,orig.ident) #keep only the cell name and the sample name
cell_coord <- left_join(cell_coord, cell_to_orig, by = "cell") # add the sample name to the cell_coord table
rm(cell_to_orig, metadata_megadata)

# Join cell_coord and metadata_seumerged by orig.ident, col and row to get the regions, subregions and regions md for each cell
# The sample name (orig.ident), col and row are needed to give each cell to the right region
cell_coord <- cell_coord %>%
  left_join(
    df_bin_number,
    by = c("orig.ident", "col", "row") 
  )
incell_transcript_coord <- incell_transcript_coord %>% 
  left_join(cell_coord %>% 
              select(cell,bin_number),
            by = "cell")
incell_transcript_coord_bin <- incell_transcript_coord %>% 
  group_by(bin_number,gene) %>% 
  summarise(x = sum(x)) %>% 
  ungroup()

saveRDS(incell_transcript_coord, file.path(data_path,"inboundaries_transcripts_cell.rds"))
saveRDS(incell_transcript_coord_bin, file.path(data_path,"inboundaries_transcripts_bin.rds"))

#initiating seurat list
seulist_extracell <- list()

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
  transcript_genexbin$bin_number <- as.character(transcript_genexbin$bin_number)
  transcript_genexbin$bin_number <- paste(im, transcript_genexbin$bin_number, sep = "_")
  
  transcript_genexbin <-  transcript_genexbin %>% 
    left_join(incell_transcript_coord_bin,
              by = c("bin_number","gene"))
  transcript_genexbin$x[is.na(transcript_genexbin$x)] <- 0
  transcript_genexbin <- transcript_genexbin %>% 
    mutate(n = n-x) %>% 
    select(bin_number,gene,n)
  
  # make the tibble wider to have nb of row = nb of genes and nb of column = nb of bins
  transcript_genexbin <- transcript_genexbin %>% 
    pivot_wider(names_from = bin_number, values_from = n) %>% 
    arrange(gene) %>% 
    column_to_rownames(var = "gene")
  
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
  
  seu <- CreateSeuratObject(
    counts = transcript_genexbin_mat,
    assay = "RNA",
    meta.data = bin_grid
  )
  
 seulist_extracell[[im]] <- seu
}

rm(seu, transcript_coord, transcript_genexbin, transcript_genexbin_mat)

# Initializing the list
seulist_extracell_log5 <- list()

# Going through the list
for (im in sample_names){
  seulist_extracell_log5[[im]] <- seulist_extracell[[im]] %>% 
    subset(subset = nCount_RNA>5) %>% #filter out bins with <30 transcripts
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData()
}

seumerged_extracell <- merge(seulist_extracell_log5[[1]], seulist_extracell_log5[-1])

metadata_seumerged_extracell <- seumerged_extracell@meta.data

df_ILAE_reg <- metadata_seumerged %>% 
  select(smpl_ILAE, ILAE_score, cluster_precast, regions, regions_md, regions_md_hilus) %>% 
  rownames_to_column(var = "bin_name")

metadata_seumerged_extracell <- metadata_seumerged_extracell %>% 
  rownames_to_column(var="bin_name") %>% 
  left_join(df_ILAE_reg,
            by = "bin_name") %>% 
  column_to_rownames(var = "bin_name")
  
seumerged_extracell@meta.data <- metadata_seumerged_extracell


seumerged_extracell <- subset(seumerged_extracell, subset = !is.na(regions))

saveRDS(seumerged_extracell, file.path(data_path,"seumerged_extraboundaries.rds"))
saveRDS(seulist_extracell, file.path(data_path,"seulist_extraboundaries.rds"))
saveRDS(seulist_extracell_log5, file.path(data_path,"seulist_extraboundaries_log5.rds"))