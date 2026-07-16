#' Perform Differential Expression Analysis using Pseudobulked Data between ILAE scores
#' Perform kruskal-wallis as a non-parametric statistical test followed by pairwise Dunn test for specific genes
#' 

# Load required packages
library(Seurat)
library(EnhancedVolcano)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(here)
library(FSA)
library(ggpubr)
library(grid)
library(colorspace)


## setting path
data_path <- str_remove(here(),"/xenium_epilepsy_repo")
fig_path <- file.path(data_path, "figures_v6")

## color palette
col_regions <- setNames(c("#3D0F99","#C7B8E6","#0F8299","#B8DEE6","#54990F","#CFE6B8","#333333","#990F26","#CCAA7A","#99600F","#FFD300"),
                        c("GCL","ML","CA2_CA4","CA1","SL.SR.SLM","GABA","SUB","WM","pia","vascular","void"))

## loading files
seumerged_extrabound <- read_rds(file.path(data_path,"seumerged_extraboundaries.rds"))
## Find all DEG between ILAE score

# joining layers for seurat functions
seumerged_extrabound_joined <- JoinLayers(seumerged_extrabound)

# removing void (and NAs)
seumerged_extrabound_joined <- subset(x = seumerged_extrabound_joined, subset = regions != "void")

#you use ILAE_score or ILAE_LK which will have Mathias ILAE scores for FF samples
pseudo_seu_extrabound <- AggregateExpression(seumerged_extrabound_joined, assays = "RNA", return.seurat = T, 
                                            group.by = c("ILAE_score", "orig.ident","regions_md")) 
pseudo_seu_extrabound@meta.data <- pseudo_seu_extrabound@meta.data %>% 
  mutate(regions_md = case_when(
    regions_md == "CA2-CA4" ~ "CA2_CA4",
    TRUE ~ regions_md
  ),
  ILAE_score = case_when(
    ILAE_score == "g0" ~ "ILAE 0",
    ILAE_score == "g1" ~ "ILAE 1",
    ILAE_score == "g2" ~ "ILAE 2",
    TRUE ~ ILAE_score
  ),
  ILAE_score = factor(ILAE_score, c( "P", "ILAE 0", "ILAE 1", "ILAE 2")))

#check data to see if grouping makes sense
tail(Cells(pseudo_seu_extrabound))

# Optional: look into how the cases clusters for each region using PCA only - expect neuropil/neurons/rest to cluster together
# pseudo_PCA <- pseudo_seu %>% 
#   FindVariableFeatures() %>% 
#   RunPCA(npcs = 50)

#add metadata column which will be a combination of region and ILAE grade
pseudo_seu_extrabound$ILAE_reg <- paste(pseudo_seu_extrabound$ILAE_score, pseudo_seu_extrabound$regions_md, sep = "_")

#Set Idents() to the new ILAE_clust combo to test DE genes within one cluster between ILAE grades
Idents(pseudo_seu_extrabound) <- "ILAE_reg"

# function to find deg test each ILAE_score against all other ILAE score and export a volcano plot as PNG and or PDF
run_FindMarkers_volcplot <- function(pseudo.data, cat, cluster, output.dir, PNG = F, PDF = T){
  library(EnhancedVolcano)
  ILAE_score <- unique(pseudo.data[[cat]][,1])
  regions <- unique(pseudo.data[[cluster]][,1])
  comb_idents <- combn(ILAE_score, 2)
  for (r in regions){
    output_save <- file.path(output.dir,r)
    dir.create(output_save, recursive = TRUE)
    temp_mat <- apply(comb_idents, 2, function(col) paste(col, r,sep = "_"))
    for (i in (1:ncol(comb_idents))){
      bulk.GC.de <- FindMarkers(object = pseudo.data,
                                ident.1 = temp_mat[1,i],
                                ident.2 = temp_mat[2,i],
                                test.use = "DESeq2")
      p <- EnhancedVolcano(bulk.GC.de,
                      lab = rownames(bulk.GC.de),  # Label genes
                      x = "avg_log2FC",  # X-axis: log2 fold change
                      y = "p_val_adj",  # Y-axis: adjusted p-value
                      drawConnectors = TRUE,  # Draw lines to labels
                      labSize = 4, pCutoff = 0.05) 
      p <- p + ggtitle(paste("DE between", temp_mat[1,i], "and", temp_mat[2,i], "in region", r, "extracell transcripts"))
      
      # Define the filename
      filename <- file.path(output_save, paste0("Volcano_", temp_mat[1,i], "_vs_", temp_mat[2,i], "_region_", r))
      
      # Save the volcano plot as a PDF
      if (PDF == T){
        pdf(file.path(filename,".pdf"), width = 10, height = 8)
        print(p)
        dev.off()
      }
      
      if (PNG == T){
        png(file.path(filename,".png"))
        print(p)
        dev.off()
      }
    }
  }
}

# function to find deg between each ILAE_score, and export CSV of genes with log2FC<-1 or log2FC>1 and padj<0.05
run_FindMarkers_regions <- function(pseudo, cat, cluster, r  = "ML", output.dir){
  output_save <- file.path(output.dir,r)
  ILAE_score <- levels(pseudo[[cat]][,1])
  regions <- unique(pseudo[[cluster]][,1])
  comb_idents <- combn(ILAE_score, 2)
  df <- tibble(gene = character(), p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), group1 = character(), group2 = character())
  temp_mat <- apply(comb_idents, 2, function(col) paste(col, r,sep = "_"))
  for (i in (1:ncol(comb_idents))){
    bulk.GC.de <- FindMarkers(object = pseudo,
                              ident.1 = temp_mat[1,i],
                              ident.2 = temp_mat[2,i],
                              test.use = "DESeq2")
    bulk.GC.de <- bulk.GC.de %>% 
      mutate(group1 = temp_mat[1,i],
             group2 = temp_mat[2,i]) %>% 
      filter((avg_log2FC < -1 | avg_log2FC > 1) & p_val_adj < 0.05) %>% 
      rownames_to_column(var = "gene")
    df <- rbind(df, bulk.GC.de)
  }
  write.csv(df, file.path(output_save, paste0("table_deg_", r, ".csv")))
}

# select the output direction
output_dir <- file.path(fig_path, "deg_pseudobulk","extracell_transcripts")

run_FindMarkers_volcplot(pseudo.data = pseudo_seu_extrabound, "ILAE_score","regions_md", output.dir = output_dir)

run_FindMarkers_regions(pseudo = pseudo_seu_extrabound, "ILAE_score", "regions_md", r = "GABA", output.dir = output_dir)

# function finding deg numbers of pseudo bulk data with the category needed to be compared specified, here ILAE_score, as well as the type of regions and which ones we want to compare for the selected category
run_FindMarkers_stackbox <- function(pseudo, cat = "ILAE_score", cluster = "regions_md", regions){
  ILAE_score <- levels(pseudo[[cat]][,1])
  comb_idents <- combn(ILAE_score, 2)
  # table storing the deg count up and down regulated
  df_1 <- tibble(region = character(), comparison = character(), FC = character(), n = integer())
  # table with all combination of regions x comparison x down/up regulation
  df_2 <- merge(data.frame(region = regions),data.frame(comparison = paste0(comb_idents[1,],"_vs_",comb_idents[2,])))
  df_2 <- merge(df_2,data.frame(FC = c("down","up")))
  for (r in regions){
    temp_mat <- apply(comb_idents, 2, function(col) paste(col, r,sep = "_"))
    for (i in (1:ncol(comb_idents))){
      bulk.GC.de <- FindMarkers(object = pseudo,
                                ident.1 = temp_mat[1,i],
                                ident.2 = temp_mat[2,i],
                                test.use = "DESeq2")
      bulk.GC.de <- bulk.GC.de %>% 
        mutate(comparison = paste0(comb_idents[1,i],"_vs_",comb_idents[2,i])) %>% 
        filter((avg_log2FC < -1 | avg_log2FC > 1) & p_val_adj < 0.05) %>% 
        mutate(FC = case_when(
          avg_log2FC < -1 ~ "down",
          TRUE ~ "up"
        )) %>% 
        group_by(comparison, FC) %>% 
        tally() %>% 
        ungroup() %>% 
        mutate(region = r) %>% 
        select(region,c(-region))
      df_1 <- rbind(df_1, bulk.GC.de)
    }
  }
  df_2 <- df_2 %>% 
    left_join(df_1, by = c("region","comparison","FC")) %>% 
    mutate(across(n,~replace_na(.,0)),
           region_FC = paste0(region, "_",FC))
  return(df_2)
}

# change order level ILAE score to match correct comparison 
# pseudo_seu_extrabound@meta.data <- pseudo_seu_extrabound@meta.data %>% 
  # mutate(ILAE_score = factor(ILAE_score, c( "ILAE 1", "ILAE 2", "ILAE 0","P")))
neuropil_regions <- c("ML", "CA2_CA4", "GABA", "SL.SR.SLM")
df <- run_FindMarkers_stackbox(pseudo = pseudo_seu_extrabound, "ILAE_score", "regions_md", regions = neuropil_regions)
df <- df %>% 
  mutate(region = factor(region, levels = rev(neuropil_regions)),
         comparison = factor(comparison, levels = c("ILAE 0_vs_P", "ILAE 1_vs_P", "ILAE 2_vs_P", "ILAE 1_vs_ILAE 0", "ILAE 2_vs_ILAE 0", "ILAE 1_vs_ILAE 2")),
         n = case_when(
           FC == "down" ~ n * -1,
           TRUE ~ n
         ))
saveRDS(df, file.path(data_path,"degcount_region_ILAEcomp.rds"))

# lighten color palette of regions
col_regions_neuropil <- col_regions[neuropil_regions]
col_regions_np_light <- lighten(col_regions_neuropil, 0.3)

# combine regular and lighten color palette for display down and up regulated gene counts
fill_palette <- c()
for (ct in neuropil_regions) {
  fill_palette[paste0(ct, "_up")] <- col_regions_np_light[ct]
  fill_palette[paste0(ct, "_down")] <- col_regions_neuropil[ct]
}

p <- ggplot(df, aes(x = n, y = region, fill = region_FC)) +
  geom_bar(stat = "identity") +
  facet_wrap(~comparison, nrow = 3) +
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

pdf(file = file.path(fig_path,"figure2_panels","deg_number_ILAEcomp.pdf"), width = 5, height = 4)
print(p)
dev.off()

## Look into specific genes to run statistical tests

color_ILAE <- setNames(c("#C67BFF","#F8766C","#7BAD00","#00C0C5"),
                       c("P", "ILAE 0","ILAE 1","ILAE 2"))

ILAE_comparison_pseudobulk <- function(pseudo, region_col = "regions_md", loc, gene, color, adj = 0.05, l.width = 0.2){
  
  # get the data for the wanted gene only
  df <- FetchData(pseudo, vars = c(gene, region_col, "ILAE_score"), assay = Assays(pseudo))
  
  # isolate for the specific region
  df_sub <- filter(df, regions_md == loc)
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
  g <- grobTree(textGrob(kw_pval, x = 1, y = 1 + adj, hjust = 1))
  
  p <- ggplot(df_sub, aes(x = ILAE_score, y = df_sub[,1], fill = ILAE_score))+
    geom_boxplot(lwd = l.width, colour = "black")+ #can be adjusted
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none"
          )+
    scale_fill_manual(values = color)+
    geom_jitter(color="black", size=0, alpha=0.9, width = 0.15)+
    ggtitle(loc) +
    annotation_custom(g)+
    theme(axis.title = element_blank())+
    # stat_compare_means(method = "kruskal.test",
    #                    label.y = y_limit + 0.45,                     # place above data
    #                    label.x = 3.6) +.
    coord_cartesian(ylim = c(0,4), clip = "off") # can be changed with the next line
    # coord_cartesian(ylim = c(0,y_limit), clip = "off")
  
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

df <- ILAE_comparison_pseudobulk(pseudo = pseudo_seu_extrabound, loc = "GABA", gene = "GRIN2C", color = color_ILAE, adj = 0.04, , l.width = 0.2)

## run ILAE_comparison_pseudobulk for the genes defined in ML, SL.SR.SLM and CA2_CA4 
export_gene_comparison <- function(list_genes, PDF = T, w = 4, h = 5, PNG = F, output_dir){
  for (region in names(list_genes)){
    save_path <- file.path(output_dir,region)
    dir.create(path = save_path, recursive = T)
    genes <- list_genes[[region]]
    for (g in genes){
      gene_comp <- ILAE_comparison_pseudobulk(pseudo = pseudo_seu_extrabound, loc = region, gene = g, color = color_ILAE, adj = 0.04)
      f.name = file.path(save_path, paste0("expression_", g, "_", region))
      if (PDF == T){
        pdf(file = paste0(f.name,".pdf"), width = w, height = h)
        print(gene_comp$p)
        dev.off()
      }
      if(PNG == T){
        png(filename = paste0(f.name,".png"), width = 480*w, height = 480*h)
        print(gene_comp$p)
        dev.off()  
      }
    }
  }
}

## list of gene to compare between ILAE score from ML, SL.SR.SLM and CA2_CA4 based on volcano plot and table above
list_genes_to_compare <- list(
  ML = c("GABBR1", "GABRG1", "GPNMB", "GRIN2C","NTS","TRIL","TSHZ2"),
  SL.SR.SLM = c("ANO3", "GABBR1", "GABRA2", "GABRG1", "GRIN2C", "MCTP2", "MEPE", "MYO5B", "SLC17A6", "TOP2A", "TRHDE", "TRIL", "TRPC6", "TTYH1"),
  CA2_CA4 = c("ABCC9", "ANO3", "GABBR1", "GABRA2", "GABRG1", "GRIN1", "GRIN2C", "MCTP2", "PVALB", "RIT2", "TRHDE", "TRIL", "TRPC6", "TSHZ2", "TTYH1"),
  GABA = c("GABRA2", "GABRG1","MGST1","MCTP2","RIMS1", "RIMS2","TRIL")
)

list_GABRG1 <- list(
  ML = c("GABRG1"),
  SL.SR.SLM = c("GABRG1"),
  CA2_CA4 = c("GABRG1"),
  GABA = c("GABRG1")
)

list_ANO3 <- list(
  SL.SR.SLM = c("ANO3"),
  CA2_CA4 = c("ANO3")
)

list_MCTP2 <- list(
  SL.SR.SLM = c("MCTP2"),
  GABA = c("MCTP2")
)

output_dir <- file.path(fig_path,"figure2","panels_mainfigure","panelG")
export_gene_comparison(list_genes = list_MCTP2, output_dir = output_dir, w = 1, h = 2.5)
