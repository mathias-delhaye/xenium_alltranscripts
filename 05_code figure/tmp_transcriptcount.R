a <- seumerged@meta.data %>% 
  group_by(cluster_precast) %>% 
  summarise(mean_count = mean(nCount_RNA))

p <- seumerged@meta.data %>% 
  ggplot(aes(x = cluster_precast, y = nCount_RNA)) +
  geom_boxplot()

p <- a %>% 
  ggplot(aes(x = cluster_precast, y = mean_count))+
  geom_col()
