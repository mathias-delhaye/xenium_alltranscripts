#' ui for the user to manual curate the regions.
#' Mathias Delhaye, July 2025

# Load libraries
library(tidyverse)
library(Seurat)
library(here)
library(future)
library(scCustomize)
library(PRECAST)
library(pals)
library(shiny)
library(plotly)
library(DT)

# Set paths
path <- here()
data_path <- str_remove(path,"xenium_epilepsy_repo")

# Load files
seumerged_log30_v6 <- readRDS(file.path(data_path,"seumerged_log30_v6.rds"))

## Analysis to interactivaly change identity of spots on graph

metadata <- seumerged_log30_v6@meta.data
metadata$spot_id <- rownames(metadata)

ui <- fluidPage(
  titlePanel("Interactive Spot Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("case", "Select Case:", choices = unique(metadata$orig.ident)),
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
  metadata_reactive <- reactiveVal(metadata)
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
    metadata <<- final_metadata  # Update global variable
  })
}

shinyApp(ui, server)

regions_md <- metadata %>%
  rename(regions_md = regions) %>% 
  select(spot_id,regions_md)

metadata_seumerged <- seumerged_log30_v6@meta.data %>% 
  rownames_to_column(var = "spot_id") %>% 
  left_join(regions_md, by = "spot_id") %>% 
  column_to_rownames(var = "spot_id")

seumerged_log30_v6@meta.data <- metadata_seumerged
rm(metadata_seumerged, regions_md)
saveRDS(seumerged_log30_v6, file.path(data_path,"seumerged_log30_v6.rds"))

# Save the scatter plots of regions_md for the samples that were curated.
dir.create(file.path(fig_path,"cluster_spa_location/k_23/regions_md"), recursive = TRUE)

for (smpl in sample_names){
  png(file.path(fig_path,"cluster_spa_location/k_23/regions_md",paste0(smpl,"_v6_k23_regions_md.png")))
  
  df <- seumerged_log30_v6@meta.data %>% 
    filter(orig.ident == smpl)
  
  # Get the plot limits
  x_min <- min(df$col)
  x_max <- max(df$col)
  y_min <- min(df$row)
  y_max <- max(df$row)

  p <-  df %>% 
    ggplot(aes(x = col, y = row, color = regions_md))+
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