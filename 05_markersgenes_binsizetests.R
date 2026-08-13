# Extract PRECASTObj for desired k
PRECASTObj <- precast_list_k1to30$k26

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
  filter(!is.na(cluster))

# add the ilae score to the metadata + adding back the rownames
metadata_seumerged <- metadata_seumerged %>% 
  mutate(orig.ident = case_when(
    orig.ident == "M_56" ~ "M56",
    TRUE ~ orig.ident
  )) %>% 
  rename(cluster_precast = cluster)

metadata_seumerged <- metadata_seumerged %>% 
  column_to_rownames(var = "spot_id")
seumerged@meta.data <- metadata_seumerged

n_count_cluster <- metadata_seumerged %>% 
  group_by(orig.ident,cluster_precast) %>% 
  summarise(mean_n_count = mean(nCount_RNA))
p <- n_count_cluster %>% 
  ggplot(aes(x = cluster_precast, y = mean_n_count))+
  stat_summary(fun = mean, geom = "bar", fill = "lightgray", color = "black") + 
  # 2. Draw the individual points
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_minimal()
p


# Join layers for seurat functions
seumerged_joined <- JoinLayers(seumerged)

# Remove void (and NAs)
seumerged_joined <- subset(x = seumerged_joined, subset = !(cluster_precast %in% c(1, 15)))

marker.features<- c('GRM2','NTNG2','KIT','WIF1','GRM4','GAD2','RELN', 'MAG','TGFBI','TP73')

p <- DotPlot(seumerged_joined, features = marker.features, group.by = 'cluster_precast') +
  scale_color_gradientn(colors = c("cyan", "white", "red")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.position = "top",
        #legend.box = "vertical",
        legend.justification = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))
