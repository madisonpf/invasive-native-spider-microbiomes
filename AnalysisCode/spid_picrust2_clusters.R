####### HOUSEKEEPING ######
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dbplyr)
library(pals)

ec_data <- read_tsv("./SERDP/Invasive_Spider_Fecal_Run/Strat_Out/Functional_Pathway_Files/pred_metagenome_contrib.tsv.gz")
ec_descrip <- read_tsv("./SERDP/Invasive_Spider_Fecal_Run/Strat_Out/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
head(ec_data)
head(ec_descrip)
as.data.frame(ec_data)
as.data.frame(ec_descrip)

merged_contrib_function <- merge(ec_data, ec_descrip, by = "function", all = TRUE)

columns_to_keep <- c("function", "sample", "taxon", "taxon_abun", "taxon_rel_abun", "genome_function_count", "taxon_function_abun", "taxon_rel_function_abun", "norm_taxon_function_contrib", "description") # Replace with your column names
filtered_merge_ec <- merged_contrib_function %>%
  select(all_of(columns_to_keep))

tax_data <- read.csv("./SERDP/Invasive_Spider_Fecal_Run/data/merged_16S_tax_clusters.csv")
names(tax_data)[names(tax_data) == 'X'] <- 'taxon'

merged_picrust <- merge(filtered_merge_ec, tax_data, by = "taxon", all = TRUE)



########################## Filtering by most abundant genera ################################

genera_of_interest <- c(" g__Acinetobacter", " g__Blastomonas", " g__Chryseobacterium", " g__Corynebacterium", " g__Pseudomonas",
                        " g__Ralstonia", " g__Reyranella", " g__Sediminibacterium", " g__Stenotrophomonas", " g__Wolbachia", " g__Rickettsia", " g__Rickettsiella", " g__Candidatus Cardinium")
filtered_merged_pi <- merged_picrust %>%
  filter(Genus %in% genera_of_interest)

top_10_within_genus <- filtered_merged_pi %>%
  group_by(Genus) %>%                    # Group by Genus
  arrange(desc(taxon_rel_function_abun), .by_group = TRUE) %>%  # Sort by Abundance within each group
  slice_head(n = 10)                    # Select top 10 rows per Genus

# Ungroup the dataframe (optional)
top_10_within_genus <- top_10_within_genus %>%
  ungroup()

write.csv(top_10_within_genus, "./SERDP/Invasive_Spider_Fecal_Run/data/picrust_filtered_results.csv", row.names = TRUE)


# Calculate mean abundances for each Genus-Function pair
mean_abundance <- picrust_filtered_results %>%
  group_by(Genus, taxon_function_abun, description) %>% # Replace "Group" with your grouping variable
  summarise(MeanAbundance = mean(taxon_function_abun, na.rm = TRUE), .groups = "drop")

# Plot grouped bar plot
ggplot(mean_abundance, aes(x = reorder(Genus, -MeanAbundance), y = MeanAbundance, fill = description)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = pals::brewer.bupu(18)) + # Adjust colors to match your groups
  labs(
    title = "",
    x = "Genus",
    y = "Mean Abundance",
    fill = "Functional Description"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),  # Increase y-axis text size
        axis.title.x = element_text(size = 15, face = "bold"),  # Increase x-axis title size
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12),   # Increase legend text size
        legend.title = element_text(size = 15, face = "bold")
  )
# Calculate mean abundances for each Genus-Function pair
mean_abundance_grouped <- picrust_filtered_results %>%
  group_by(Genus, taxon_function_abun, functional.grouping) %>% # Replace "Group" with your grouping variable
  summarise(MeanAbundance = mean(taxon_function_abun, na.rm = TRUE), .groups = "drop")

# Plot grouped bar plot
ggplot(mean_abundance_grouped, aes(x = reorder(Genus, -MeanAbundance), y = MeanAbundance, fill = functional.grouping)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = pals::brewer.brbg(9)) + # Adjust colors to match your groups
  
  labs(
    title = "",
    x = "Genus",
    y = "Mean Abundance",
    fill = "Functional Groups"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),  # Increase y-axis text size
        axis.title.x = element_text(size = 15, face = "bold"),  # Increase x-axis title size
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12),   # Increase legend text size
        legend.title = element_text(size = 15, face = "bold")
  )
