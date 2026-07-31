##################################################
## Script purpose: Run Differential Expression Analysis using Pseudobulked Data between ILAE scores
## Date: 2026-07-28
## Update date: 2026-07-28
## Author: Mathias Delhaye
###################################################

##################################################
## Section 1: Libraries
##################################################

library(EnhancedVolcano)
library(tidyverse)
library(DESeq2)
library(here)
library(FSA)
library(ggpubr)
library(grid)
library(colorspace)
library(Matrix)
library(future)
library(future.apply)
library(Seurat)
library(ggh4x)

##################################################
## Section 2: Environment
##################################################

# Set paths
choose.files()
data_path <- choose.dir()
rds_path <- file.path(data_path, "01_rds_files")
fig_path <- file.path(data_path, "02_plots")

# Load files
megadata <- read_rds(file.path(rds_path, "01_MegaData_v8.rds"))
seulist <- readRDS(file.path(rds_path, "02_seulist_v8.rds"))
samples_names <- names(seulist)

col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))
color_ILAE <- setNames(c("#C67BFF","#F8766C","#7BAD00","#00C0C5","#E69F00"),
                       c("P", "ILAE0","ILAE1","ILAE2","ILAE3"))

# Genes analyzed with xenium (266 panel + 100 epil specific)
genes_list <- Features(megadata)

# Merge all per-sample layers into a single counts layer
megadata[["Xenium"]] <- JoinLayers(megadata[["Xenium"]])

##################################################
## Section 3: Extract in-boundaries transcripts for each bin 
##################################################

# Pull the full counts matrix once
xen_count <- LayerData(megadata, assay = "Xenium", layer = "counts")

# Build the transcript coordinate table in one shot
temp_df <- as.data.frame(summary(xen_count))
cell_names <- colnames(xen_count)     # use colnames/rownames directly, not gene_list
gene_names <- rownames(xen_count)

incell_transcript_coord <- tibble(
  gene = gene_names[temp_df$i],
  cell = cell_names[temp_df$j],
  n_count    = temp_df$x
)

metadata_seumerged <- seumerged@meta.data

df_bin_number <- metadata_seumerged %>% 
  rownames_to_column(var = "bin_number") %>% 
  select(orig.ident, col, row, bin_number)

cell_coord <- map_dfr(samples_names, function(smpl) {
  coord <- megadata@images[[smpl]]@boundaries[["centroids"]]@coords
  cell_names <- megadata@images[[smpl]]@boundaries[["centroids"]]@cells
  tibble(cell = cell_names, x = coord[,1], y = coord[,2])
})

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
  summarise(n_count = sum(n_count)) %>% 
  ungroup()

saveRDS(incell_transcript_coord, file.path(rds_path,"09_inboundaries_transcripts_percell.rds"))
saveRDS(incell_transcript_coord_bin, file.path(rds_path,"10_inboundaries_transcripts_perbin.rds"))

##################################################
## Section 4: Isolate extra-boundaries only transcripts
##################################################

build_extracell_seurat <- function(im, gene_list, megadata,
                                   incell_bin_counts,
                                   bin_width = 150, bin_height = 150) {
  
  fov <- megadata@images[[im]]@molecules$molecules
  
  # Pull every gene's coords once, combine with a single rbind (not one per gene)
  coords_list <- lapply(gene_list, function(g) fov[[g]]@coords)
  n_per_gene  <- vapply(coords_list, nrow, integer(1))
  all_coords  <- do.call(rbind, coords_list)
  gene_vec    <- rep(gene_list, n_per_gene)
  
  bin_x <- ceiling(all_coords[, 1] / bin_width)
  bin_y <- ceiling(all_coords[, 2] / bin_height)
  
  number_bin_x <- ceiling(max(all_coords[, 1]) / bin_width)
  number_bin_y <- ceiling(max(all_coords[, 2]) / bin_height)
  n_bins <- number_bin_x * number_bin_y
  
  bin_number <- (bin_y - 1) * number_bin_x + bin_x
  gene_idx   <- match(gene_vec, gene_list)
  bin_names  <- paste0(im, "_", seq_len(n_bins))
  
  # Total transcripts per gene x bin (sparse; duplicate i,j pairs are summed automatically)
  mat_total <- sparseMatrix(
    i = gene_idx, j = bin_number, x = 1,
    dims = c(length(gene_list), n_bins),
    dimnames = list(gene_list, bin_names)
  )
  
  # Matching in-cell counts, aligned to the same gene/bin dimensions
  incell_sub <- incell_bin_counts[incell_bin_counts$bin_number %in% bin_names, ]
  mat_incell <- sparseMatrix(
    i = match(incell_sub$gene, gene_list),
    j = match(incell_sub$bin_number, bin_names),
    x = incell_sub$n_count,   # <-- confirm this is the right column name
    dims = c(length(gene_list), n_bins),
    dimnames = list(gene_list, bin_names)
  )
  
  mat_extracell <- mat_total - mat_incell
  mat_extracell[mat_extracell < 0] <- 0  # safety guard, see note below
  
  # col fastest, row slowest -- matches your original bin_number formula
  bin_grid <- expand.grid(col = 1:number_bin_x, row = 1:number_bin_y)
  rownames(bin_grid) <- bin_names
  
  CreateSeuratObject(counts = mat_extracell, assay = "RNA", meta.data = bin_grid)
}

seulist_extracell <- setNames(
  lapply(samples_names, build_extracell_seurat,
         gene_list = genes_list,
         megadata = megadata,
         incell_bin_counts = incell_transcript_coord_bin),
  samples_names
)

saveRDS(seulist_extracell, file.path(rds_path,"11_seulist_extraboundaries_v8.rds"))

plan(multisession, workers = parallel::detectCores() - 1)  # adjust worker count as needed

seulist_extracell_log5 <- future_lapply(seulist_extracell[samples_names], function(seu) {
  seu %>%
    subset(subset = nCount_RNA > 5) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE)
}, future.seed = TRUE)

names(seulist_extracell_log5) <- samples_names
plan(sequential)  # reset back to sequential when done
saveRDS(seulist_extracell_log5, file.path(rds_path,"12_seulist_extraboundaries_log5_v8.rds"))

seumerged_extracell <- merge(seulist_extracell_log5[[1]], seulist_extracell_log5[-1])

metadata_seumerged_extracell <- seumerged_extracell@meta.data

df_ILAE_reg <- metadata_seumerged %>% 
  select(smpl_ILAE, ILAE_score, cluster_precast, region, region_md, region_md_hilus) %>% 
  rownames_to_column(var = "bin_name")

metadata_seumerged_extracell <- metadata_seumerged_extracell %>% 
  rownames_to_column(var="bin_name") %>% 
  left_join(df_ILAE_reg,
            by = "bin_name") %>% 
  column_to_rownames(var = "bin_name")
  
seumerged_extracell@meta.data <- metadata_seumerged_extracell

saveRDS(seumerged_extracell, file.path(rds_path,"13_seumerged_extraboundaries_v8.rds"))

##################################################
## Section 5: Run pseudobulking 
##################################################

# Load seumerged_extrabound if needed
seumerged_extrabound <- readRDS(file.path(rds_path,"13_seumerged_extraboundaries_v8.rds"))

# Join layers for seurat functions
seumerged_extrabound_joined <- JoinLayers(seumerged_extrabound)

# removing void (and NAs)
seumerged_extrabound_joined <- subset(x = seumerged_extrabound_joined, subset = region != "void")

#you use ILAE_score or ILAE_LK which will have Mathias ILAE scores for FF samples
pseudo_seu_extrabound <- AggregateExpression(seumerged_extrabound_joined, assays = "RNA", return.seurat = T, 
                                             group.by = c("ILAE_score", "orig.ident","region_md")) 

pseudo_seu_extrabound@meta.data <- pseudo_seu_extrabound@meta.data %>% 
  mutate(region_md = case_when(
    region_md == "CA2-CA4" ~ "CA2_CA4",
    TRUE ~ region_md
  ),
  ILAE_score = case_when(
    ILAE_score == "g0" ~ "ILAE0",
    ILAE_score == "g1" ~ "ILAE1",
    ILAE_score == "g2" ~ "ILAE2",
    ILAE_score == "g3" ~ "ILAE3",
    TRUE ~ ILAE_score
  ),
  ILAE_score = factor(ILAE_score, c( "P", "ILAE0", "ILAE1", "ILAE2", "ILAE3")))

# Check data to see if grouping makes sense
tail(Cells(pseudo_seu_extrabound))

# Add metadata column which will be a combination of region and ILAE grade
pseudo_seu_extrabound$ILAE_reg <- paste(pseudo_seu_extrabound$ILAE_score, pseudo_seu_extrabound$region_md, sep = "_")

# Set Idents() to the new ILAE_clust combo to test DE genes within one cluster between ILAE grades
Idents(pseudo_seu_extrabound) <- "ILAE_reg"

##################################################
## Section 6: run DEG analysis on pseudobulk data
##################################################

# helper function to compute the DEG table once
compute_deg_table <- function(pseudo,
                              cat = "ILAE_score",
                              cluster = "region_md",
                              regions = NULL,
                              group_sep = "_",
                              fc_cutoff = 1,
                              padj_cutoff = 0.05,
                              verbose = TRUE) {
  
  ILAE_score  <- levels(pseudo[[cat]][, 1])
  if (is.null(regions)) regions <- unique(pseudo[[cluster]][, 1])
  comb_idents <- combn(ILAE_score, 2)
  
  results_list <- list()
  
  for (r in regions) {
    temp_mat <- apply(comb_idents, 2, function(col) paste(col, r, sep = group_sep))
    
    for (i in seq_len(ncol(comb_idents))) {
      ident1 <- temp_mat[1, i]
      ident2 <- temp_mat[2, i]
      comparison <- paste0(comb_idents[1, i], "_vs_", comb_idents[2, i])
      
      if (verbose) message("DE: ", ident1, " vs ", ident2)
      
      de <- tryCatch(
        FindMarkers(object = pseudo, ident.1 = ident1, ident.2 = ident2,
                    test.use = "DESeq2"),
        error = function(e) {
          warning("Skipped ", ident1, " vs ", ident2, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(de)) next
      
      de <- de %>%
        tibble::rownames_to_column(var = "gene") %>%
        mutate(
          region     = r,
          group1     = ident1,
          group2     = ident2,
          comparison = comparison,
          direction  = case_when(
            avg_log2FC >  fc_cutoff & p_val_adj < padj_cutoff ~ "up",
            avg_log2FC < -fc_cutoff & p_val_adj < padj_cutoff ~ "down",
            TRUE ~ "ns"
          ),
          significant = direction != "ns"
        )
      
      results_list[[paste(r, ident1, ident2, sep = "|")]] <- de
    }
  }
  
  bind_rows(results_list)
}

saveRDS(df_deg, file = file.path(rds_path,"14_degtable_v8.rds"))

# helper function to generate volcano plots and upload them to a specific directory
plot_deg_volcanoes <- function(deg_table, output.dir,
                               pCutoff = 0.05, labSize = 4, w = 6, h = 6) {
  
  for (r in unique(deg_table$region)) {
    output_save <- file.path(output.dir, r)
    dir.create(output_save, recursive = TRUE, showWarnings = FALSE)
    
    combos <- df_deg %>% filter(region == r) %>% distinct(group1, group2)
    
    for (j in seq_len(nrow(combos))) {
      g1 <- combos$group1[j]; g2 <- combos$group2[j]
      
      sub_de <- df_deg %>%
        filter(region == r, group1 == g1, group2 == g2) %>%
        column_to_rownames("gene")
      
      p <- EnhancedVolcano(sub_de,
                           lab = rownames(sub_de),
                           x = "avg_log2FC", y = "p_val_adj",
                           drawConnectors = TRUE, 
                           subtitle = NULL,
                           captionLabSize = 12,
                           legendLabSize = 10,
                           titleLabSize = 14,
                           pCutoff = 0.05) +
        ggtitle(paste(g1, "vs", g2, " - ", r))
      p
      filename <- file.path(output_save, paste0("Volcano_", g1, "_vs_", g2, "_region_", r,".png"))
      
      ggsave(
        filename = filename,
        plot = p,
        dpi = 300,
        width = 6,
        height = 6
      )
    }
  }
}

# helper function to export CSV of significant genes
export_deg_csv <- function(deg_table, output.dir) {
  for (r in unique(deg_table$region)) {
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    
    df <- deg_table %>%
      filter(region == r, significant) %>%
      select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, group1, group2)
    
    write.csv(df, file.path(output.dir, paste0("table_deg_", r, ".csv")), row.names = FALSE)
  }
}

# helper function to summarize deg table for stacked plot
summarize_deg_stackbox <- function(deg_table, regions = NULL) {
  if (is.null(regions)) regions <- unique(deg_table$region)
  comparisons <- unique(deg_table$comparison)
  
  full_grid <- expand.grid(region = regions, comparison = comparisons,
                           FC = c("down", "up"), stringsAsFactors = FALSE)
  
  counts <- deg_table %>%
    filter(significant) %>%
    mutate(FC = direction) %>%
    group_by(region, comparison, FC) %>%
    tally() %>%
    ungroup()
  
  full_grid %>%
    left_join(counts, by = c("region", "comparison", "FC")) %>%
    mutate(n = replace_na(n, 0),
                  region_FC = paste0(region, "_", FC))
}

# run functions on pseudo_seu_extrabound
df_deg <- compute_deg_table(
  pseudo = pseudo_seu_extrabound,
  cat = "ILAE_score",
  cluster = "region_md"
  )

csv_path <- file.path(data_path,"03_csv_deg")
export_deg_csv(df_deg,csv_path)

output_dir <- file.path(fig_path, "06_deg_pseudobulk_extracell_transcripts")
plot_deg_volcanoes(df_deg, output.dir = output_dir)

df_deg_neuropil_summarized <- summarize_deg_stackbox(df_deg, neuropil_region) %>%
  mutate(region = factor(region, levels = rev(neuropil_region)),
         comparison = factor(comparison, levels = c("P_vs_ILAE0", "P_vs_ILAE1", "P_vs_ILAE2", "P_vs_ILAE3",
                                                    "ILAE0_vs_ILAE1", "ILAE0_vs_ILAE2", "ILAE0_vs_ILAE3",
                                                    "ILAE1_vs_ILAE2", "ILAE1_vs_ILAE3", "ILAE2_vs_ILAE3")),
         n = case_when(FC == "down" ~ n * -1, TRUE ~ n)) %>%
  separate(comparison, into = c("group1", "group2"), sep = "_vs_", remove = FALSE) %>%
  mutate(group1 = factor(group1, levels = c("P", "ILAE0", "ILAE1", "ILAE2")),
         group2 = factor(group2, levels = c("ILAE0", "ILAE1", "ILAE2", "ILAE3")))


# lighten color palette of region
col_region_neuropil <- col_regions[neuropil_region]
col_region_np_light <- lighten(col_region_neuropil, 0.3)

# combine regular and lighten color palette for display down and up regulated gene counts
fill_palette <- c()
for (ct in neuropil_region) {
  fill_palette[paste0(ct, "_up")] <- col_region_np_light[ct]
  fill_palette[paste0(ct, "_down")] <- col_region_neuropil[ct]
}

df_deg_neuropil_summarized <- summarize_deg_stackbox(df_deg, neuropil_region) %>%
  mutate(region = factor(region, levels = rev(neuropil_region)),
         comparison = factor(comparison, levels = c("P_vs_ILAE0", "P_vs_ILAE1", "P_vs_ILAE2", "P_vs_ILAE3",
                                                    "ILAE0_vs_ILAE1", "ILAE0_vs_ILAE2", "ILAE0_vs_ILAE3",
                                                    "ILAE1_vs_ILAE2", "ILAE1_vs_ILAE3", 
                                                    "ILAE2_vs_ILAE3")),
         n = case_when(FC == "down" ~ n * -1, TRUE ~ n))

# Each letter below is one panel, assigned in the order of `comparison`'s factor
# levels. "#" = no panel drawn at all (blank space, not even a grey background).
#
#        ILAE0 ILAE1 ILAE2 ILAE3
#   P      A     B     C     D
# ILAE0    #     E     F     G
# ILAE1    #     #     H     I
# ILAE2    #     #     #     J
design <- "
ABCD
#EFG
##HI
###J
"

p <- ggplot(df_deg_neuropil_summarized, aes(x = n, y = region, fill = region_FC)) +
  geom_bar(stat = "identity") +
  facet_manual(vars(comparison), design = design) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = fill_palette, guide = "none") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = c("#cccccc3b"), colour = NA),
    strip.text = element_text(size = 8),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
  )

ggsave(
  filename = file.path(output_dir,"degcount_stack.png"),
  plot = p,
  width = 8,
  height = 6
)

p <- ggplot(df, aes(x = n, y = region, fill = region_FC)) +
  geom_bar(stat = "identity") +
  facet_wrap(~comparison, nrow = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = fill_palette) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = c("#cccccc3b"), colour = NA),
    strip.text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
  )

##################################################
## Section 7: Analysis of specific genes
##################################################

ILAE_comparison_pseudobulk <- function(pseudo, region_col = "region_md", loc, gene, color, l.width = 0.2){
  
  # get the data for the wanted gene only
  df <- FetchData(pseudo, vars = c(gene, region_col, "ILAE_score"), assay = Assays(pseudo))
  
  # isolate for the specific region
  df_sub <- filter(df, region_md == loc)
  names(df_sub)[1] <- "expression"
  
  kw <- kruskal.test(df_sub[,1], g = df_sub[,3])
  dunn <- dunnTest(expression ~ ILAE_score, data = df_sub, method = "bh")
  dunn <- dunn$res %>% 
    mutate(significance = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01 ~ "**",
      P.adj < 0.05 ~ "**",
      TRUE ~ "ns"
    )) %>% 
    separate(Comparison, into = c("group1", "group2"), sep = " - ")
  dunn_plot <- dunn %>% 
    filter(P.adj<0.05)
  
  max_df <-  max(df[,1])
  y_limit <- max_df+0.3*4 #should accommodate 4 significant comparisons
  kw_pval <- paste("p.val",as.character(round(kw$p.value, 4)))
  
  p <- ggplot(df_sub, aes(x = ILAE_score, y = df_sub[,1], fill = ILAE_score))+
    geom_boxplot(lwd = l.width, colour = "black",outlier.shape = NA)+ #can be adjusted
    theme_minimal()+
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black",
                                     size = 9),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_blank(),
          title = element_text(colour = "black",
                                     size = 9),
    )+
    scale_fill_manual(values = color)+
    geom_jitter(color="black", size=0, alpha=0.9, width = 0.15)+
    labs(title = paste0(loc, " - ",kw_pval)) +
    coord_cartesian(ylim = c(0,y_limit), clip = "off")
  
  if (nrow(dunn_plot)>0){
    position_comparison = seq(from = max_df + 0.3, to = max_df + nrow(dunn_plot)*0.3, by = 0.3)
    p <- p +
      stat_pvalue_manual(dunn_plot, label = "significance", tip.length = 0.01, y.position = position_comparison)
  }
  
  result_list <- list(
    p = p ,
    kw = kw,
    dunn = dunn
  )
  return(result_list)
}

test <- ILAE_comparison_pseudobulk(pseudo_seu_extrabound,loc = "CA2_CA4",gene = "GABBR1", color = color_ILAE)
test$p
ggsave(
  filename = file.path(fig_path,"test.png"),
  plot = test$p,
  width = 2,
  height = 4,
  dpi = 300
)

## run ILAE_comparison_pseudobulk for the genes defined in ML, SL.SR.SLM and CA2_CA4 
export_gene_comparison <- function(list_genes,w = 2, h = 4, output_dir){
  for (region in names(list_genes)){
    save_path <- file.path(output_dir,region)
    dir.create(path = save_path, recursive = T, showWarnings = F)
    genes <- list_genes[[region]]
    for (g in genes){
      gene_comp <- ILAE_comparison_pseudobulk(pseudo = pseudo_seu_extrabound, loc = region, gene = g, color = color_ILAE)
      filename = file.path(save_path, paste0("expression_", g, "_", region,".png"))
      ggsave(
        filename = filename,
        plot = gene_comp$p,
        width = w,
        height = h,
        dpi = 300
      )
    }
  }
}

## list of gene to compare between ILAE score from ML, SL.SR.SLM and CA2_CA4 based on volcano plot and table above
list_genes_to_compare <- list(
  ML = c("ANO3","GABBR1", "GABRG1", "GPNMB", "GRIN2C","NTS","TRIL","TSHZ2","CALB1","SLC17A6"),
  SL.SR.SLM = c("ANO3", "GABBR1", "GABRA2", "GABRG1", "GRIN2C","GRM5", "MCTP2", "MEPE", "MYO5B", "SLC17A6", "TOP2A", "TRHDE", "TRIL", "TRPC6", "TTYH1"),
  CA2_CA4 = c("ABCC9", "ANO3","CALCRL", "GABBR1", "GABRA2", "GABRG1","GRM5", "GRIN1","GRIN2B","GRIN2C", "MCTP2", "PVALB", "RIT2","SLC24A3", "TRHDE", "TRIL", "TRPC6", "TSHZ2", "TTYH1"),
  GABA = c("ANO3","GRIN2C","GABRA2", "GABRG1","MGST1","MCTP2","RIMS1", "RIMS2","TRIL","CALB2")
)
output_dir <- file.path(fig_path,"07_deg_pergene")
export_gene_comparison(list_genes = list_genes_to_compare, output_dir = output_dir)
