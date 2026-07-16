n_spots_regions <- metadata_seumerged %>% 
  pivot_longer(cols = starts_with("regions"), names_to = "spa_cluster", values_to = "regions") %>% 
  group_by(orig.ident, spa_cluster, regions) %>% 
  tally(name = "n_spot") %>% 
  ungroup() %>% 
  filter(regions != "void") %>%
  mutate(area_um2 = n_spot * 150 * 150) %>% 
  group_by(orig.ident, spa_cluster) %>% 
  mutate(ratio_regions = 100*n_spot/sum(n_spot))

saveRDS(n_spots_regions,file.path(data_path,"n_spots_regions_v6.rds"))

n_cells_regions <- metadata_megadata %>% 
  group_by(orig.ident,ILAE_score,type, celltype, regions_md) %>% 
  tally() %>% 
  ungroup() %>% 
  filter(regions_md != "void") %>% 
  left_join(n_spots_regions,
            by = c("orig.ident", "regions_md")) %>% 
  mutate(cell_dens = n/n_spot)

excn_per_regions <- n_cells_regions %>% 
  filter(celltype == "Oligo") %>% 
  group_by(orig.ident, ILAE_score, regions_md) %>% 
  summarise(avg_dens = mean(cell_dens)) %>% 
  ungroup() %>% 
  ggplot(aes(x = orig.ident, y = avg_dens, fill = regions_md))+
  geom_bar(stat = "identity") +
  ylab("excn/spot") +
  scale_fill_manual(name = "Regions", values = col_regions) +
  theme_minimal(base_size = 10)