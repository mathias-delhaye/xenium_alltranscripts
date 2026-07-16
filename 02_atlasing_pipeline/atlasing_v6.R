## run PRECAST on seulist_log30 for k = 23 clusters

library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)
library(grid)
library(shiny)
library(plotly)
library(DT)

## setting path
data_path <- str_remove(here(),"/xenium_epilepsy_repo")
fig_path <- file.path(data_path, "figures_v6")

## Loading files
megadata <- read_rds(file.path(data_path,"MegaData_v6.rds"))

# name of all samples in megadata in order assigned in megadata
sample_names <- names(megadata@images)
# all the genes analyzed with xenium (aka 266 panel + 100 epil specific)
gene_list <- Features(megadata)

## 01 - megadata to seulist
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
saveRDS(seulist,file.path(data_path,"seulist_v6.rds"))

## Building seulist_log30 filter with >30 counts/spot + log normalization

# Initializing the list
seulist_log30 <- list()

# Going through the list
for (im in sample_names){
  seulist_log30_v6[[im]] <- seulist[[im]] %>% 
    subset(subset = nCount_RNA>30) %>% #filter out bins with <30 transcripts
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData()
}

saveRDS(seulist_log30,file.path(data_path,"seulist_log30_v6.rds"))

## 02 - seulist to PRECAST
## Create PRECAST object
preobj <- CreatePRECASTObject(seuList = seulist_log30)

saveRDS(preobj, file.path(data_path,"preobj_v6.rds"))

## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
PRECASTObj_log30 <- AddAdjList(preobj, platform = "Other_SRT")

## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj_log30 <- AddParSetting(PRECASTObj_log30, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

saveRDS(PRECASTObj_log30, file.path(data_path,"PRECASTObj_log30.rds"))

# k = 23
PRECASTObj_log30_k23 <- PRECAST(PRECASTObj_log30, K = 23)

saveRDS(PRECASTObj_log30_k23, file.path(data_path,"PRECASTObj_log30_k23_v6.rds"))

## 03 - PRECAST to seuint
## creating integrated seurat object

# Selecting model
seuint <- PRECASTObj_log30_k23 %>% 
  SelectModel() %>% 
  IntegrateSpaData(species = "unknown")

# locate the position of specific clusters in each tissue section
cluster_location <- function(seuint, orig, cluster_nb, k_cluster){
  sub <- subset(seuint, subset = orig.ident == orig)
  sub@meta.data <- sub@meta.data %>% mutate(batch = orig)
  cols_sub <- brewer.blues(n = k_cluster)
  cols_sub[cluster_nb] = "#d62728FF"
  sub %>% SpaPlot(item = "cluster", batch = NULL, cols = cols_sub, point_size = 3.5, nrow.legend = 26, title_name = "")
}

# running through all clusters for all samples
for (smpl in sample_names){
  dir.create(file.path(fig_path,"highlight_clusters_position_v6/k23_spaplot",smpl), recursive = TRUE) # allow to create subfolders if needed
  path = file.path(fig_path,"highlight_clusters_position_v6/k23_spaplot",smpl)
  for (i in (1:23)){
    png(file.path(path,paste0(smpl,"_cluster",as.character(i),"_k23_sans_ILAE3.png")))
    print(cluster_location(seuint_log30_k23_v6,smpl,i, 23)) #needs to be explicitely printed
    dev.off()
  } 
}

# adding the regions and subregions information to the metadata of seuInt
seuint@meta.data <- seuint@meta.data %>% 
  mutate(batch = orig.ident, #so I can use batch for the naming of the figures
         regions = case_when(
           cluster == 7 | cluster == 11 ~ "GCL",
           cluster == 23 ~ "ML",
           cluster == 12 | cluster == 18 | cluster == 20 ~ "CA2_CA4",
           cluster == 1 | cluster == 2 ~ "CA1",
           cluster == 3 | cluster == 9 | cluster == 17 ~ "SL.SR.SLM",
           cluster == 8 ~ "GABA",
           cluster == 5 | cluster == 13 ~ "SUB",
           cluster == 4 | cluster == 15 | cluster == 16 ~ "WM",
           cluster == 21 ~ "pia",
           cluster == 10 | cluster == 14 | cluster == 22 ~ "vascular",
           cluster == 6 | cluster == 19 ~ "void"))

saveRDS(seuint, file.path(data_path,"seuint_log30_k23_v6.rds"))

## 04 - seuint + seulist to seumerged

# defining gradient of colors for the regions and subregions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))


# seulist to seumerged
#merging all the elements of the list containing the seurat objects for each case
seumerged <- merge(seulist_log30[[1]], seulist_log30[-1])

# generate table with ILAE score
ILAE_table <- megadata@meta.data %>% 
  group_by(orig.ident,ILAE_score) %>% 
  tally() %>% 
  select(orig.ident,ILAE_score) %>% 
  mutate(ILAE_score = factor(ILAE_score, levels = c("P","0","1","2")))

# add the ilae score to the metadata + adding back the rownames
metadata_seuint <- seuint@meta.data %>% 
  select(cluster, regions) %>% 
  rownames_to_column(var = "id")
metadata_seumerged <- seumerged@meta.data %>% 
  rownames_to_column(var = "id") %>% 
  left_join(ILAE_table, by = "orig.ident") %>% 
  left_join(metadata_seuint, by = "id") %>% 
  rename(cluster_precast = cluster) %>% 
  column_to_rownames(var = "id") 
seumerged@meta.data <- metadata_seumerged

## Analysis to interactivaly change identity of spots on graph

metadata_seumerged$spot_id <- rownames(metadata_seumerged)

ui <- fluidPage(
  titlePanel("Interactive Spot Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("case", "Select Case:", choices = unique(metadata_seumerged$orig.ident)),
      actionButton("change_id", "Change Identity of Selected Points"),
      textInput("new_id", "New Identity:", value = "")
    ),
    mainPanel(
      plotlyOutput("scatterPlot"),
      DTOutput("selectedPoints")
    )
  )
)

server <- function(input, output, session) {
  metadata_reactive <- reactiveVal(metadata_seumerged)
  filteredData <- reactive({
    metadata_reactive() %>% filter(orig.ident == input$case)
  })
  
  output$scatterPlot <- renderPlotly({
    plot_ly(
      data = filteredData(),
      x = ~col,
      y = ~row,
      type = 'scatter',
      mode = 'markers',
      color = ~regions,
      colors = col_regions,
      text = ~paste("regions:", regions, "<br>Spot:", spot_id),
      key = ~spot_id,
      marker = list(size = 6),
      source = "scatter"
    ) %>% layout(yaxis = list(autorange = "reversed"))
  })
  
  selected_spots <- reactive({
    event_data("plotly_selected", source = "scatter")
  })
  
  output$selectedPoints <- renderDT({
    if (is.null(selected_spots())) return(NULL)
    spots <- selected_spots()$key
    filteredData() %>% filter(spot_id %in% spots)
  })
  
  observeEvent(input$change_id, {
    sel <- selected_spots()
    if (!is.null(sel) && input$new_id != "") {
      spots <- sel$key
      updated <- metadata_reactive() %>%
        mutate(regions = ifelse(spot_id %in% spots, input$new_id, regions))
      metadata_reactive(updated)  # Update reactiveVal
    }
  })
  
  session$onSessionEnded(function() {
    # Save metadata to global env (or disk) when session ends
    final_metadata <- isolate(metadata_reactive())  # Use isolate() to safely access it
    metadata_seumerged <<- final_metadata  # Update global variable
  })
}

shinyApp(ui, server)

regions_md_hilus <- metadata_seumerged %>%
  rename(regions_md_hilus = regions) %>% 
  select(spot_id,regions_md_hilus)

seumerged@meta.data <- seumerged@meta.data %>% 
  rownames_to_column(var = "spot_id") %>% 
  left_join(regions_md_hilus, by = "spot_id") %>% 
  column_to_rownames(var = "spot_id")

rm(metadata_seumerged, regions_md, metadata_seuint)

seumerged$regions <- factor(seumerged$regions, levels = c("GCL", "ML", "CA2_CA4", "CA1", "SUB", "GABA", "SL.SR.SLM", "WM", "vascular", "pia", "void"))
seumerged$regions_md <- factor(seumerged$regions_md, levels = c("GCL", "ML", "CA2_CA4", "CA1", "SUB", "GABA", "SL.SR.SLM", "WM", "vascular", "pia", "void"))
seumerged$regions_md_hilus <- factor(seumerged$regions_md_hilus, levels = c("GCL", "ML", "hilus","CA2_CA4", "CA1", "SUB", "GABA", "SL.SR.SLM", "WM", "vascular", "pia", "void"))
seumerged@meta.data <- seumerged@meta.data %>% 
  unite(smpl_ILAE, c(orig.ident, ILAE_score), sep = "_", remove = FALSE)
seumerged$smpl_ILAE <- factor(seumerged$smpl_ILAE, levels = c("L10_P", "L5_P", "M70_P", "B1_0", "F1_0", "H1_0", "I1_0", "J1_0", "L20_0", "C1_1", "C2_1", "E008_1", "L14_1", "L15_1", "L18_1", "L1_2", "L21_2", "L3_2"))

saveRDS(seumerged, file.path(data_path,"seumerged_log30_v6.rds"))

# plotting the regions
dir.create(file.path(fig_path,"cluster_spa_location/k_23/regions"), recursive = TRUE)

# regions
for (smpl in sample_names){
  png(file.path(fig_path,"cluster_spa_location/k_23/regions",paste0(smpl,"_v6_k23_regions.png")))
  
  df <- seumerged@meta.data %>% 
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

# regions_md
for (smpl in sample_names){
  png(file.path(fig_path,"cluster_spa_location/k_23/regions_md",paste0(smpl,"_v6_k23_regions_md.png")))
  
  df <- seumerged@meta.data %>% 
    filter(orig.ident == smpl)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)
  
  p <-  df %>% 
    ggplot(aes(x = col, y = row, color = regions_md))+
    geom_point(size = 2.75)+
    scale_color_manual(name = "Regions_md", values = col_regions) +
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
         color = "Region_md") + 
    annotation_custom(
      grob = rectGrob(
        gp = gpar(col = "black", fill = NA, lwd = 1)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  print(p) #needs to be explicitely printed
  dev.off()
}

seumerged_sansNA <- subset(seumerged, subset = !is.na(regions_md))
# Highlight clusters location with seumerged
for (smpl in sample_names){
  dir.create(file.path(fig_path,"highlight_clusters_position_v6/k23_seumerged",smpl), recursive = TRUE)
  
  # keep metadata for smpl only
  df <- seumerged_sansNA@meta.data %>% 
    filter(orig.ident == smpl)
  # get the precast clusters contained by smpl
  df$cluster_precast <- droplevels(df$cluster_precast)
  # finding missing clusters
  missing_clusters <- setdiff(as.character(1:23), levels(df$cluster_precast))
  
  for (i in levels(df$cluster_precast)){
    png(file.path(fig_path,"highlight_clusters_position_v6/k23_seumerged",smpl,paste0(smpl,"cluster_",i,"_v6_k23.png")))
    col_clusters <- brewer.blues(n = 23)
    col_clusters[as.numeric(i)] = "#d62728FF"
    # remove missing clusters
    if (length(missing_clusters)!=0){
      col_clusters <- col_clusters[-as.numeric(missing_clusters)]
    }
    
    # Get the plot limits
    x_min <- min(df$col)
    x_max <- max(df$col)
    y_min <- min(df$row)
    y_max <- max(df$row)
    
    p <-  df %>% 
      ggplot(aes(x = col, y = row, color = cluster_precast))+
      geom_point(size = 2.75)+
      scale_color_manual(name = "cluster_precast", values = col_clusters) +
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
      labs(title = smpl) + 
      guides(color = guide_legend(title = "precast_clusters",ncol = 1))+
      annotation_custom(
        grob = rectGrob(
          gp = gpar(col = "black", fill = NA, lwd = 1)
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      )
    print(p) #needs to be explicitely printed
    dev.off()
  }
}

## 05 - Regions to megadata
## adding regions and subregions to metadata of megadata

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
  rownames_to_column(var = "cell_names") %>% 
  select(cell_names,orig.ident) #keep only the cell name and the sample name
cell_coord <- left_join(cell_coord, cell_to_orig, by = "cell_names") # add the sample name to the cell_coord table
rm(cell_to_orig)

# Join cell_coord and metadata_seumerged by orig.ident, col and row to get the regions, subregions and regions md for each cell
# The sample name (orig.ident), col and row are needed to give each cell to the right region
metadata_seumerged <- seumerged@meta.data
cell_coord <- cell_coord %>%
  left_join(
    metadata_seumerged %>%
      select(orig.ident, col, row, regions, regions_md, regions_md_hilus),
    by = c("orig.ident", "col", "row") 
  )

# Join the metadata of megadata and seumerged together using the cell_names
metadata_megadata <- metadata_megadata %>% 
  rownames_to_column(var = "cell_names") %>% 
  left_join(cell_coord %>% 
              select(cell_names, regions, regions_md, regions_md_hilus),
            by = "cell_names") %>% 
  column_to_rownames(var = "cell_names")

# add column to metadata of megadata merging orig.ident and ILAE
metadata_megadata <- metadata_megadata %>% 
  unite(smpl_ILAE, c(orig.ident, ILAE_score), sep = "_", remove = FALSE)

saveRDS(metadata_megadata, file.path(data_path,"metadata_md_v6.rds"))