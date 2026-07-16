# load seulist_6 object
seulist <- read_rds(file.path(data_path,"seulist_v6.rds"))
seumerged <- read_rds(file.path(data_path,"seumerged_log30_v6.rds"))


# initiate table where each row will be a case & location identified bin, for a gene and it's number of transcript in the identified bin
transcripts_per_bin <- tibble(bin_number = character(), gene = character(), x = numeric(), stringsAsFactors = FALSE) #it has to be x

for (smpl in sample_names){
  RNA_count <- seulist[[smpl]]@assays$RNA$counts
  temp_df <- as.data.frame(summary(RNA_count)) # the RNA count is under the name x
  bin_names <- RNA_count@Dimnames[[2]]
  temp_df$bin_number <- bin_names[temp_df$j]
  temp_df$gene <- gene_list[temp_df$i]
  temp_df <- temp_df[, c("bin_number", "gene", "x")]
  transcripts_per_bin <- rbind(transcripts_per_bin, temp_df)
}

rm(RNA_count, temp_df, bin_names, smpl)
# load inboundaries_transcript_bin
inboundaries_transcripts_bins <- read_rds(file.path(data_path,"inboundaries_transcripts_bin.rds"))

transcripts_per_bin <- transcripts_per_bin %>% 
  mutate(transcripts_tot = x) %>% #here change x by transcript tot which makes more sense
  select(bin_number, gene, transcripts_tot) %>%
  left_join(inboundaries_transcripts_bins, by = c("bin_number", "gene"))
transcripts_per_bin$x[is.na(transcripts_per_bin$x)] <- 0

# get info about case, regions and regions_md that will be added the transcripts_per_bin
metadata_seumerged <- seumerged@meta.data
df_bin_to_regions <- metadata_seumerged %>% 
  rownames_to_column(var = "bin_number") %>% 
  select(bin_number, orig.ident, regions, regions_md)

transcripts_per_bin <- left_join(transcripts_per_bin, df_bin_to_regions, by = "bin_number")
transcripts_per_bin <- transcripts_per_bin %>% 
  filter(!is.na(regions)) %>% #only keep bins that were identified as a region
  mutate(inbound_transcript = x) %>% #rename x variable
  select(!x)

rm(metadata_seumerged)

# get ratio of transcript inbound for each case

transcript_ratio_case <- transcripts_per_bin %>%
  filter(regions!= "void") %>% 
  group_by(orig.ident) %>% 
  summarise(ratio_in = mean(inbound_transcript/transcripts_tot)) %>% 
  ungroup() %>% 
  mutate(ratio_out = 1-ratio_in)

# generate table with ILAE score
ILAE_table <- megadata@meta.data %>% 
  group_by(orig.ident,ILAE_score) %>% 
  tally() %>% 
  select(orig.ident,ILAE_score) %>% 
  mutate(ILAE_score = factor(ILAE_score, levels = c("P","0","1","2")))

transcript_ratio_case <- transcript_ratio_case %>%
  left_join(ILAE_table, by = "orig.ident") %>%
  mutate(orig_ILAE = paste(orig.ident, ILAE_score, sep = "_"),
         orig_ILAE = factor(orig_ILAE, levels = c("L5_P","L10_P","M70_P","B1_0","F1_0","H1_0","I1_0","J1_0","L20_0","C1_1","C2_1","E008_1","L14_1","L15_1","L18_1","L1_2","L3_2","L21_2")))

# Suppose your data frame is called df
# First, pivot the two ratio columns to a long format
df_long <- transcript_ratio_case %>%
  pivot_longer(cols = c(ratio_out, ratio_in),
               names_to = "ratio_type",
               values_to = "ratio_value")

# Create the stacked bar plot
p <- ggplot(df_long, aes(x = orig_ILAE, y = ratio_value, fill = ratio_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ratio_in" = "gold", "ratio_out" = "#54990F")) + # optional colors
  labs(x = "Sample", y = "Ratio", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
png(filename = "C:/Users/Mango/OneDrive - UBC/Presentation/MC lab/lab meetings/Fig_labmeeting_20251027/ratio_transcripts_cases.png", height = 300, width = 600)
print(p)
dev.off()

# get ratio of transcript inbound for each case
transcript_ratio_region <- transcripts_per_bin %>%
  filter(regions!= "void") %>% 
  group_by( regions) %>% 
  summarise(ratio_in = mean(inbound_transcript/transcripts_tot)) %>% 
  ungroup() %>% 
  mutate(ratio_out = 1-ratio_in)%>% 
  pivot_longer(cols = c(ratio_out, ratio_in),
               names_to = "ratio_type",
               values_to = "ratio_value")

p <- ggplot(transcript_ratio_region, aes(x = regions, y = ratio_value, fill = ratio_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ratio_in" = "gold", "ratio_out" = "#54990F")) + # optional colors
  labs(x = "Sample", y = "Ratio", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
png(filename = "C:/Users/Mango/OneDrive - UBC/Presentation/MC lab/lab meetings/Fig_labmeeting_20251027/ratio_transcripts_regions.png", height = 300, width = 450)
print(p)
dev.off()
