##################################################
## Script purpose: generate regional domains based on all transcripts detected by Xenium. run PRECAST on seulist_log30 for k = 23 clusters
## Date: 2026-07-16
## Update date: 2026-07-28
## Author: Mathias Delhaye
##################################################

##################################################
## Section 1: Libraries
##################################################

library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)
library(grid)
library(plotly)
library(DT)
library(shiny)

##################################################
## Section 2: Environment
##################################################

# Set paths
choose.files()
data_path <- choose.dir()
rds_path <- file.path(data_path, "01_rds_files")
fig_path <- file.path(data_path, "02_plots")

# Load files
megadata <- read_rds(file.path(rds_path,"01_MegaData_v8.rds"))

# Names of every sample in megadata in order assigned in megadata
samples_names <- names(megadata@images)

# Genes analyzed with xenium (266 panel + 100 epil specific)
genes_list <- Features(megadata)

##################################################
## Section 3: megadata to seulist
##################################################

# Build seulist with raw transcript count

# Initiating seurat list
seulist <- list()

# Set bin dimension arbitrarily - 150*150px works fine
bin_width <- 150
bin_height <- 150

for (im in samples_names){
  # Field of view where all the genes are stored individually, containing the coordinates of all the transcripts
  fov <- megadata@images[[im]]@molecules$molecules
  
  # Create a tibble to store the genes + their coordinates
  transcript_coord <- tibble(gene = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)
  
  for (n in (1:length(genes_list))){
    gene <- genes_list[n] # go through the list of genes
    gene_matrix <- fov[[gene]]@coords # select in the fov the gene
    temp_df <- tibble(gene = gene, x = gene_matrix[,1], y = gene_matrix[,2]) # temporary df to add all transcripts with their coordinates (transcript name being just the gene name)
    transcript_coord <- rbind(transcript_coord,temp_df) # join the temp df with the tibble containing all genes
    rm(temp_df)
  }
  
  # Assign each transcript to a spot
  transcript_coord$bin_x <- ceiling(transcript_coord$x / bin_width) # x coordinate / width of bin
  transcript_coord$bin_y <- ceiling(transcript_coord$y / bin_height) # y coordinate / height of bin
  
  # Find max number of bin per column and row
  number_bin_x <- ceiling(max(transcript_coord$x)/bin_width)
  number_bin_y <- ceiling(max(transcript_coord$y)/bin_height)
  
  # Give the order of bins
  transcript_coord$bin_number <- (transcript_coord$bin_y-1)*number_bin_x+transcript_coord$bin_x # (bin_y coordinate - 1 aka number of full lines) * number of bins in X  + coordinate bin_x for the number of bins to add to the full lines
  
  # Count the number of transcripts for each bin and each gene
  transcript_genexbin <- transcript_coord %>% 
    group_by(bin_number, gene) %>% 
    tally() %>% 
    ungroup()
  
  # Make the tibble wider to have nb of row = nb of genes and nb of column = nb of bins
  transcript_genexbin <- transcript_genexbin %>% 
    pivot_wider(names_from = bin_number, values_from = n) %>% 
    arrange(gene) %>% 
    column_to_rownames(var = "gene") %>% 
    rename_with(~ paste0(im,"_", .)) #change the name of the bin to add the sample name _ and not just having numbers
  
  # Find the bins that we kept through previous operation, which are bins with at least 1 transcript
  bin_w_transcript <- colnames(transcript_genexbin)
  # Find the missing bins and add empty columns with these missing bins
  missing_bin = setdiff(paste0(im,"_",as.character(1:(number_bin_x*number_bin_y))),bin_w_transcript)
  transcript_genexbin[,missing_bin] <- NA
  
  # Convert tibble into matrix
  transcript_genexbin_mat <- as.matrix(transcript_genexbin)
  
  # Replace NA by 0
  transcript_genexbin_mat[is.na(transcript_genexbin_mat)] <- 0
  
  # Re-order the column to have them in the order of bins
  transcript_genexbin_mat <- transcript_genexbin_mat[,paste0(im,"_",as.character(1:(number_bin_x*number_bin_y)))]
  
  # Get the coordinates for each bin
  bin_grid <- expand_grid(row = 1:number_bin_y, col = 1:number_bin_x) %>% # create all combination possible with the value given (here being the number of bins in X and Y), with the second argument being the one used first to increase value (1-1, 1-2, 1-3, ... 1-n, 2-1, 2-2, ... 2-n, etc...)
    select(col,row) %>% # re-order the columns to have x first then y
    as.data.frame() # convert into a df
  # Name each row with the corresponding bin
  row.names(bin_grid) <- paste0(im,"_",as.character(1:(number_bin_x*number_bin_y)))
  # Add the position windows for each bin in X & Y so that we can place segmented area in each bin based on the centroids position
  bin_grid$x_window <- bin_grid$col*bin_width
  bin_grid$y_window <- bin_grid$row*bin_height
  
  seu <- CreateSeuratObject(
    counts = transcript_genexbin_mat,
    assay = "RNA",
    meta.data = bin_grid
  )
  
  # Change orig.ident to set it to sample name - needs to use level as it is considered as a factor
  levels(seu@meta.data$orig.ident) <- im
  
  seulist[[im]] <- seu
}

# Save seurat list
saveRDS(seulist,file.path(rds_path,"02_seulist_v8.rds"))

seulist <- read_rds(file.path(rds_path,"02_seulist_v8.rds"))

# Building seulist_log30 filter with >30 counts/spot + log normalization

seulist_log30 <- list()

# Going through the list
for (im in samples_names){
  seulist_log30[[im]] <- seulist[[im]] %>% 
    subset(subset = nCount_RNA>30) %>% #filter out bins with <30 transcripts
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData()
}

saveRDS(seulist_log30,file.path(rds_path,"03_seulist_log30_v8.rds"))

##################################################
## Section 4: seulist to PRECAST
##################################################

# Load seulist_log30 if needed
seulist_log30 <- readRDS(file.path(rds_path,"03_seulist_log30_v8.rds"))

# Create PRECAST object
preobj <- CreatePRECASTObject(seuList = seulist_log30)

saveRDS(preobj, file.path(rds_path,"04_preobj_v8.rds"))

# Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
PRECASTObj <- AddAdjList(preobj, platform = "Other_SRT")

# Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

saveRDS(PRECASTObj, file.path(rds_path,"05_PRECASTobj_v8.rds"))

# Find range of optimal number of cluster
precast_list_k1to30 = list()
for (k in 1:30){
  precast_list_k1to30[[paste0("k",as.character(k))]] <- PRECAST(PRECASTObj, K = k)
}

saveRDS(precast_list_k1to30, file.path(rds_path,"06_precast_list_k1to30.rds"))

aicList<-list()
for(k in 1:30){
  PRECASTObj <- SelectModel(precast_list_k1to30[[k]],criteria='AIC')
  aicList[k]<-PRECASTObj@resList$icMat[2]
  print(aicList[k])
  rm(PRECASTObj)
  gc()
}

aic<-as.numeric(unlist(aicList))

# Adjust the indices to start from 2
aic_df <- data.frame(Index = 1:30, AIC = aic)

# Plot using ggplot2 with y-axis labels in scientific notation
p <- ggplot(aic_df, aes(x = Index, y = AIC)) +
  geom_line() + # Add line
  geom_point()+
  scale_y_continuous(labels = scales::scientific)+
  geom_vline(xintercept = 8, linetype = "dashed", color = "red", linewidth = 1)+
  geom_vline(xintercept = 22, linetype = "dashed", color = "red", linewidth = 1)+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,
                               angle = 45)
  )
  
ggsave(filename = file.path(fig_path,"01_precast_AICvskcluster.png"),
      plot = p,
      dpi = 300,
      width = 8,
      height = 6)

# Optimal k between 21 and 23

##################################################
## Section 5: PRECAST to seuint
##################################################

# Load seulist_log30 if needed
seulist_log30 <- readRDS(file.path(rds_path,"03_seulist_log30_v8.rds"))

# Create seumerged for later use
seumerged <- merge(seulist_log30[[1]], seulist_log30[-1])

df_temp <- seumerged@meta.data %>% 
  rownames_to_column(var = "spot_id") %>% 
  select(all_of(c("spot_id","col","row")))

# Loop through desired k to attribute region to each cluster

k_opt = c(21, 22,23)

for (k in k_opt) {
  # Extract PRECASTObj for desired k
  PRECASTObj <- precast_list_k1to30[[paste0("k", as.character(k))]]
  
  # Select model
  seuint <- PRECASTObj %>% 
    SelectModel() %>% 
    IntegrateSpaData(species = "unknown")
  
  # Extract metadata and merge coordinates
  metadata_seuint <- seuint@meta.data %>% 
    mutate(orig.ident = case_when(
      orig.ident == "M" ~ "M_56",
      TRUE ~ orig.ident
    )) %>% 
    rownames_to_column(var = "id") %>% 
    left_join(df_temp, by = "id")
  
  # Define all possible cluster levels for this k
  all_clusters <- levels(factor(metadata_seuint$cluster)) 
  
  # Set base palette and assign explicit names
  col_clusters_blues <- brewer.blues(n = k)
  names(col_clusters_blues) <- all_clusters
  
  for (smpl in samples_names) {
    df <- metadata_seuint %>% 
      filter(orig.ident == smpl)
    
    # Strictly enforce factor levels
    df$cluster <- factor(df$cluster, levels = all_clusters)
    
    # --- FIX FOR LEGEND BOXES ON 0-COUNT LEVELS ---
    missing_levels <- setdiff(all_clusters, unique(df$cluster[!is.na(df$cluster)]))
    if (length(missing_levels) > 0) {
      dummy_rows <- df[rep(1, length(missing_levels)), ]
      dummy_rows$cluster <- factor(missing_levels, levels = all_clusters)
      dummy_rows$col <- NA  # NA ensures point is skipped on the plot
      dummy_rows$row <- NA
      df <- rbind(df, dummy_rows)
    }
    # -----------------------------------------------
    
    path <- file.path(fig_path, "02_clusters_location_v8", paste0("k", k), smpl)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    
    for (i in seq_along(all_clusters)) {
      current_cluster <- all_clusters[i]
      
      # Copy named palette and highlight target cluster
      col_clusters_highlighted <- col_clusters_blues
      col_clusters_highlighted[current_cluster] <- "#d62728FF"
      
      max_dim <- max(c(df$col, df$row), na.rm = TRUE)
      point_size_pt <- (6 * 72) / max_dim
      
      p <- df %>% 
        ggplot(aes(x = col, y = row, color = cluster, fill = cluster)) +
        geom_point(shape = 15, size = point_size_pt, na.rm = TRUE) +
        scale_color_manual(
          values = col_clusters_highlighted,
          drop = FALSE,
          guide = "none"
        ) +
        scale_fill_manual(
          name = "Cluster",
          values = col_clusters_highlighted,
          drop = FALSE,
          guide = guide_legend(
            ncol = 1,
            override.aes = list(
              shape = 22,      # Square shape
              color = "black", # Black border outline
              stroke = 0.5,    # Outline thickness
              size = 4
            )
          )
        ) +
        scale_y_reverse(limits = c(max_dim, 1), expand = c(0, 0)) +
        scale_x_continuous(limits = c(1, max_dim), expand = c(0, 0)) +
        coord_fixed(expand = FALSE) +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          plot.background = element_rect(fill = "white", color = NA),
          plot.title.position = "panel",
          plot.title = element_text(
            size = 20,
            hjust = 0.01,
            margin = margin(t = 0, b = -25),
            face = "bold"
          ),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10, face = "bold"),
          legend.box.margin = margin(l = -5),
          legend.margin = margin(l = 0),
          legend.key = element_blank(),
          aspect.ratio = 1
        ) +
        labs(
          title = smpl, 
          fill = "Cluster"
        ) + 
        annotation_custom(
          grob = rectGrob(
            gp = gpar(col = "black", fill = NA, lwd = 2.5)
          ),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        )
      
      ggsave(
        filename = file.path(path, paste0(smpl, "_cluster", current_cluster, "_k", k, ".png")),
        plot = p,
        dpi = 300,
        width = 6,
        height = 6
      )
    }
  }
}

# Now plot case colored by region for each case and choose best k
# Colors for the regions
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

# Extract PRECASTObj for desired k
PRECASTObj <- precast_list_k1to30$k23

# Select model
seuint <- PRECASTObj %>% 
  SelectModel() %>% 
  IntegrateSpaData(species = "unknown")

# Extract metadata from seuint
metadata_seuint <- seuint@meta.data %>% 
  select(orig.ident,cluster) %>% 
  rownames_to_column(var = "spot_id") %>% 
  mutate(orig.ident = case_when(
    orig.ident == "M" ~ "M_56",
    TRUE ~ orig.ident
  ))

# Cells that survived PRECAST integration (i.e., got a cluster assignment)
kept_cells <- Cells(seuint)
length(kept_cells)

# Filter each per-sample object to only these cells, BEFORE merging
seulist_log30_filtered <- lapply(seulist_log30, function(seu) {
  cells_to_keep <- intersect(Cells(seu), kept_cells)
  subset(seu, cells = cells_to_keep)
})

# Re-compute seumerged to omit the missing spots
seumerged <- merge(seulist_log30_filtered[[1]], seulist_log30_filtered[-1])

metadata_seumerged <- seumerged@meta.data %>% 
  rownames_to_column(var = "spot_id") %>% 
  left_join(metadata_seuint, by = c("spot_id","orig.ident")) %>%
  filter(!is.na(cluster)) %>% 
  mutate(region = case_when(
           cluster == 7 | cluster == 21 ~ "GCL",
           cluster == 1 ~ "ML",
           cluster == 6 | cluster == 11 | cluster == 16 ~ "CA2_CA4",
           cluster == 4 | cluster == 14 | cluster == 19 ~ "CA1",
           cluster == 3 | cluster == 5 | cluster == 18 ~ "SL.SR.SLM",
           cluster == 9 ~ "GABA",
           cluster == 12 ~ "SUB",
           cluster == 15 | cluster == 17 ~ "WM",
           cluster == 23 ~ "pia",
           cluster == 13 | cluster == 22 ~ "vascular",
           cluster == 2 | cluster == 8 | cluster == 10 | cluster == 20 ~ "void"))

path <- file.path(fig_path, "03_regions_location_v8", "k_23")
dir.create(path, recursive = TRUE, showWarnings = FALSE)

for (smpl in samples_names) {
  df <- metadata_seumerged %>% 
    filter(orig.ident == smpl)
  max_dim <- max(c(df$col, df$row), na.rm = TRUE)
  point_size_pt <- (6 * 72) / max_dim
  p <- df %>% 
    ggplot(aes(x = col, y = row, color = region, fill = region)) +
    geom_point(shape = 15, size = point_size_pt, na.rm = TRUE) +
    scale_color_manual(
      values = col_regions,
      drop = FALSE,
      guide = "none"
    ) +
    scale_fill_manual(
      name = "Region",
      values = col_regions,
      drop = FALSE,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(
          shape = 22,      # Square shape
          color = "black", # Black border outline
          stroke = 0.5,    # Outline thickness
          size = 4
        )
      )
    ) +
    scale_y_reverse(limits = c(max_dim, 1), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1, max_dim), expand = c(0, 0)) +
    coord_fixed(expand = FALSE) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.title.position = "panel",
      plot.title = element_text(
        size = 20,
        hjust = 0.01,
        margin = margin(t = 0, b = -25),
        face = "bold"
      ),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      legend.box.margin = margin(l = -5),
      legend.margin = margin(l = 0),
      legend.key = element_blank(),
      aspect.ratio = 1
    ) +
    labs(
      title = smpl, 
      fill = "Region"
    ) + 
    annotation_custom(
      grob = rectGrob(
        gp = gpar(col = "black", fill = NA, lwd = 2.5)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  
  ggsave(
    filename = file.path(path, paste0(smpl, "_k_21.png")),
    plot = p,
    dpi = 300,
    width = 6,
    height = 6
  )
}

# from the visualization of the domains, k = 23 fits the most

##################################################
## Section 6: seuint + seulist to seumerged
##################################################

# generate table with ILAE score
ILAE_table <- megadata@meta.data %>% 
  group_by(orig.ident,ILAE_score) %>% 
  tally() %>% 
  select(orig.ident,ILAE_score) %>% 
  mutate(ILAE_score = factor(ILAE_score, levels = c("P","0","1","2","3")))

# add the ilae score to the metadata + adding back the rownames
metadata_seumerged <- metadata_seumerged %>% 
  mutate(orig.ident = case_when(
    orig.ident == "M_56" ~ "M56",
    TRUE ~ orig.ident
  )) %>% 
  left_join(ILAE_table, by = "orig.ident") %>% 
  rename(cluster_precast = cluster) %>% 
  mutate(region_md = region)

# Analysis to interactivaly change identity of spots on graph

ui <- fluidPage(
  titlePanel("Interactive Spot Viewer"),
  tags$script(HTML("
    $(document).on('keydown', '#new_id', function(e) {
      if (e.key === 'Enter' || e.keyCode === 13) {
        e.preventDefault();
        $('#change_id').click();
      }
    });
  ")),
  sidebarLayout(
    sidebarPanel(
      selectInput("case", "Select Case:", choices = unique(metadata_seumerged$orig.ident)),
      actionButton("change_id", "Change Identity of Selected Points"),
      textInput("new_id", "New Identity:", value = "")
    ),
    mainPanel(
      plotlyOutput("scatterPlot", height = "85vh", width = "100%"),
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
    df <- filteredData()
    x_range <- range(df$col, na.rm = TRUE)
    
    plot_ly(
      data = df,
      x = ~col,
      y = ~row,
      type = 'scatter',
      mode = 'markers',
      color = ~region_md_hilus,
      colors = col_regions,
      text = ~paste("region:", region_md_hilus, "<br>Spot:", spot_id),
      key = ~spot_id,
      marker = list(size = 6, symbol = "square"),
      source = "scatter"
    ) %>%
      layout(
        autosize = TRUE,
        xaxis = list(
          scaleanchor = "y",
          scaleratio = 1,
          constrain = "domain",
          range = x_range
        ),
        yaxis = list(
          autorange = "reversed",
          constrain = "domain"
        ),
        dragmode = "pan"   # pan is now the default drag behavior
      ) %>%
      config(
        scrollZoom = TRUE,
        displaylogo = FALSE
      ) %>%
      onRender("
        function(el, x) {

          // ---- marker size scaling on zoom ----
          el.on('plotly_relayout', function(ev) {
            var xa = el.layout.xaxis;
            var rangeNow;
            if (ev['xaxis.range[0]'] !== undefined) {
              rangeNow = Math.abs(ev['xaxis.range[1]'] - ev['xaxis.range[0]']);
            } else if (ev['xaxis.autorange']) {
              rangeNow = Math.abs(xa.range[1] - xa.range[0]);
            } else {
              return;
            }
            el._baselineRange = el._baselineRange || rangeNow;
            var zoomFactor = el._baselineRange / rangeNow;
            var baseSize = 6;
            var newSize = Math.max(2, Math.min(40, baseSize * zoomFactor));
            Plotly.restyle(el, {'marker.size': newSize});
          });

          // ---- keyboard shortcut: space toggles lasso <-> pan ----
          var lassoActive = false;

          function isTyping(target) {
            var tag = target.tagName;
            return tag === 'INPUT' || tag === 'TEXTAREA' || target.isContentEditable;
          }

          document.addEventListener('keydown', function(e) {
            if (isTyping(e.target)) return;
            if (e.repeat) return; // ignore key-repeat

            if (e.code === 'Space') {
              e.preventDefault(); // stop page scroll
              lassoActive = !lassoActive;
              var newMode = lassoActive ? 'lasso' : 'pan';
              console.log('Space pressed -> switching to', newMode);
              Plotly.relayout(el, {dragmode: newMode});
            }
          });
        }
      ")
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
        mutate(region_md_hilus = ifelse(spot_id %in% spots, input$new_id, region_md_hilus))
      metadata_reactive(updated)
    }
  })
  
  session$onSessionEnded(function() {
    final_metadata <- isolate(metadata_reactive())
    metadata_seumerged <<- final_metadata
  })
}

shinyApp(ui, server)

# Add extra column for hilus - re-run the shiny app 
metadata_seumerged <- metadata_seumerged %>%
  mutate(region_md_hilus = region_md)

seumerged@meta.data <- metadata_seumerged %>% 
  column_to_rownames(var = "spot_id")

rm(metadata_seumerged, metadata_seuint)

seumerged$region <- factor(seumerged$region, levels = c("CA1", "SL.SR.SLM", "CA2_CA4", "SUB", "GCL", "ML", "GABA", "WM", "vascular", "pia", "void"))
seumerged$region_md <- factor(seumerged$region_md, levels = c("CA1", "SL.SR.SLM", "CA2_CA4", "SUB", "GCL", "ML", "GABA", "WM", "vascular", "pia", "void"))
seumerged$region_md_hilus <- factor(seumerged$region_md_hilus, levels = c("CA1", "SL.SR.SLM", "CA2_CA4", "SUB", "GCL", "ML", "GABA", "WM", "vascular", "pia", "void"))
seumerged@meta.data <- seumerged@meta.data %>% 
  unite(smpl_ILAE, c(orig.ident, ILAE_score), sep = "_", remove = FALSE)
seumerged$smpl_ILAE <- factor(seumerged$smpl_ILAE, levels = c("L5_P", "L10_P", "M56_P", "M70_P", "P28_P", "P51_P", "P60_P", "P71_P", "B1_0", "D8_0", "D11_0", "F1_0", "H1_0", "I1_0", "J1_0", "L20_0", "C1_1", "C2_1", "E008_1", "L14_1", "L15_1", "L18_1", "D12_2", "D13_2", "L1_2", "L3_2", "L21_2", "D9_3", "D14_3", "E015_3", "L2_3"))

saveRDS(seumerged, file.path(rds_path,"07_seumerged_v8.rds"))

# Plot region_md and region_md_hilus

path <- file.path(fig_path, "03_regions_location_v8", "k_23_region_md_hilus")
dir.create(path, recursive = TRUE, showWarnings = FALSE)

col_regions_hilus <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300","grey"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void", "hilus"))

for (smpl in samples_names) {
  df <- seumerged@meta.data %>% 
    mutate(orig.ident = case_when(
      orig.ident == "M56" ~ "M_56",
      TRUE ~ orig.ident
    )) %>% 
    filter(orig.ident == smpl)
  max_dim <- max(c(df$col, df$row), na.rm = TRUE)
  point_size_pt <- (6 * 72) / max_dim
  p <- df %>% 
    ggplot(aes(x = col, y = row, color = region_md_hilus, fill = region_md_hilus)) +
    geom_point(shape = 15, size = point_size_pt, na.rm = TRUE) +
    scale_color_manual(
      values = col_regions,
      drop = FALSE,
      guide = "none"
    ) +
    scale_fill_manual(
      name = "Region",
      values = col_regions_hilus,
      drop = FALSE,
      guide = guide_legend(
        ncol = 1,
        override.aes = list(
          shape = 22,      # Square shape
          color = "black", # Black border outline
          stroke = 0.5,    # Outline thickness
          size = 4
        )
      )
    ) +
    scale_y_reverse(limits = c(max_dim, 1), expand = c(0, 0)) +
    scale_x_continuous(limits = c(1, max_dim), expand = c(0, 0)) +
    coord_fixed(expand = FALSE) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.title.position = "panel",
      plot.title = element_text(
        size = 20,
        hjust = 0.01,
        margin = margin(t = 0, b = -25),
        face = "bold"
      ),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      legend.box.margin = margin(l = -5),
      legend.margin = margin(l = 0),
      legend.key = element_blank(),
      aspect.ratio = 1
    ) +
    labs(
      title = smpl, 
      fill = "Region"
    ) + 
    annotation_custom(
      grob = rectGrob(
        gp = gpar(col = "black", fill = NA, lwd = 2.5)
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  
  ggsave(
    filename = file.path(path, paste0(smpl, "_k_21.png")),
    plot = p,
    dpi = 300,
    width = 6,
    height = 6
  )
}

##################################################
## Section 7: Regions to megadata
##################################################

# create a tibble to store the genes + their coordinates
cell_coord <- tibble(cell = character(), x = numeric(), y = numeric(), stringsAsFactors = FALSE)

for (smpl in samples_names){
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
rm(cell_to_orig)

# Join cell_coord and metadata_seumerged by orig.ident, col and row to get the regions, subregions and regions md for each cell
# The sample name (orig.ident), col and row are needed to give each cell to the right region
metadata_seumerged <- seumerged@meta.data
cell_coord <- cell_coord %>%
  left_join(
    metadata_seumerged %>%
      select(orig.ident, col, row, region, region_md, region_md_hilus),
    by = c("orig.ident", "col", "row") 
  )

# Join the metadata of megadata and seumerged together using the cell_names
metadata_megadata <- metadata_megadata %>% 
  rownames_to_column(var = "cell") %>% 
  left_join(cell_coord %>% 
              select(cell, region, region_md, region_md_hilus),
            by = "cell") %>% 
  column_to_rownames(var = "cell")

saveRDS(metadata_megadata, file.path(rds_path,"08_megadata_metadata_md_v8.rds"))
