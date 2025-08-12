########### PRE-PROCESSING AND ANALYSIS SCRIPT - CLUSTERED SEQUENCES #######################

######################### PRE-PROCESSING ################################################

### HOUSEKEEEPING ####
#Installing packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
install.packages("VennDiagram")

#Loading packages
library(BiocManager)
library(phyloseq)
library(microbiome)
library(ggplot2) 
library(plyr)
library(VennDiagram)
library(dplyr)
library(grid)
library(vegan)
library(reshape2)
library(scales)
library(ggpubr)
library(cowplot)
library(compositions)
library(readr)
library(tibble)
library(DESeq2)
library(pairwiseAdonis)
library(car)
library(dplyr)
library(purrr)
library(emmeans)

setwd("~/Documents")

#LOADING IN DATA FILES - Available on Dryad
final_clusters <- read.csv("./SERDP/Invasive_Spider_Fecal_Run/data/filtered_clustered_ASV.csv", row.names = 1) #ASV file (OTUs)
final_clusters <- as.matrix(final_clusters) #Convert to matrix 

merged_16s_tax_data= read.csv("./SERDP/Invasive_Spider_Fecal_Run/data/merged_16S_tax_clusters.csv", header = TRUE, row.names=1, sep = ",") #Taxonomic Assignment file 

threshold_percentage <- 0.90 #Set threshold for filtering taxa to 90%

filtered_taxa <- merged_16s_tax_data[merged_16s_tax_data$Confidence >= threshold_percentage, ] #Filter taxa to high than 90% confidence in assignment

tax_16s = as.matrix(filtered_taxa) #Convert to matrix 

##LOAD IN METADATA
metadata_16s = read.csv("./SERDP/Invasive_Spider_Fecal_Run/spider_metadata_16s_only.csv") 
rownames(metadata_16s) <- metadata_16s$Sample_ID #Change rownames to match taxonomy file

#LOAD IN PHY_TREE
#Each phy tree from separate A,B,C runs loaded separately 
phy_tree_A = read_tree("./SERDP/Invasive_Spider_Fecal_Run/data/tree-cluster-A.nwk") 
phy_tree_B = read_tree("./SERDP/Invasive_Spider_Fecal_Run/data/tree-cluster-B.nwk")
phy_tree_C = read_tree("./SERDP/Invasive_Spider_Fecal_Run/data/tree-cluster-C.nwk")

merged_tree_data <- rbind(phy_tree_A, phy_tree_B, phy_tree_C) #phy tree files merged to match OTU and tax files 


#Make the files phyloseq objects
FT_16S = otu_table(final_clusters, taxa_are_rows = TRUE) 
TAX_16S = tax_table(tax_16s)
META_16S = sample_data(metadata_16s)

#Check names to ensure correct naming schemes
taxa_names(FT_16S) 

taxa_names(TAX_16S)

## Make sure files have the same sample names
sample_names(FT_16S) 
sample_names(META_16S)

## Merge into one phyloseq object
phy_cluster = phyloseq(FT_16S, TAX_16S, META_16S, merged_tree_data) 
phy_cluster ## Final phyloseq file 

# =======================================
###### CLEAN DATA - PRE-PROCESSING ######
# =======================================

# 1) Prune out blanks/zero‐count taxa
phy_pruned <- prune_samples(sample_sums(phy_cluster) > 0, phy_cluster)
phy_pruned <- prune_taxa(taxa_sums(phy_pruned) > 0, phy_pruned)

# 2) Extract the OTU count matrix, making sure samples are rows
otu_mat <- as(otu_table(phy_pruned), "matrix")
if(taxa_are_rows(phy_pruned)) {
  otu_mat <- t(otu_mat)
}

par(mfrow = c(1, 1))

# 3) Plot rarefaction curves
rarecurve(otu_mat,
          step     = 100,
          label    = FALSE,
          xlab     = "Reads per sample",
          ylab     = "Observed ASV richness",
          col      = "steelblue",
          cex      = 0.6,
          xlim     = c(0, 2000))

# 1) Prune out blanks and zero‐count taxa
blanks  <- c("BLANK_16S_1", "BLANK_16S_2", "PCR_BLANK_16S")
phy_real <- phy_pruned %>%
  prune_samples(!(sample_names(.) %in% blanks), .) %>%
  prune_taxa(taxa_sums(.)   > 0, .)

# 2) Examine library sizes
libs <- sample_sums(phy_pruned)
summary(libs)
min_reads <- min(libs)
cat("Min reads among real samples:", min_reads, "\n")
# e.g. Min reads among real samples: 551

# Decide on depth (either min_reads to keep all, or e.g. 800 to drop only lowest)
depth <- 1000

# 3) Rarefy
set.seed(123)
phy_rarefied <- rarefy_even_depth(
  phy_real,
  sample.size = depth,
  rngseed     = TRUE,
  replace     = FALSE,
  trimOTUs    = TRUE
)

# 4) Quick diagnostics
cat("Samples before rarefaction:", nsamples(phy_pruned), "\n")
cat("Samples after rarefaction: ", nsamples(phy_rarefied), "\n")
cat("Total reads before:", sum(libs), "\n")
cat("Total reads after: ", sum(sample_sums(phy_rarefied)), "\n")

dropped <- setdiff(sample_names(phy_real), sample_names(phy_rarefied))
cat("Dropped samples (insufficient depth):\n")
print(data.frame(
  Sample = dropped,
  Reads  = libs[dropped]
))

phy_cluster <- phy_rarefied #Assign the phyloseq object to be cleaned object for downstream

# ======================================
### SPECIES CONSISTENCY - ASV COUNTS ###
# ======================================

meta <- as(sample_data(phy_rarefied), "data.frame") %>%
  dplyr::select(Species.Identification, Island)

# get counts per genus × island
counts <- meta %>% 
  group_by(Species.Identification, Island) %>% 
  tally(name = "n")

# reshape so it’s ready for a table
counts_wide <- counts %>% 
  tidyr::pivot_wider(names_from = Island, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(-Species.Identification))) 

# add a “Total” row
grand_totals <- counts_wide %>% 
  summarise(across(-Species.Identification, sum)) %>%
  mutate(Species.Identification = "Total") %>%
  select(Species.Identification, everything())

sample_counts <- bind_rows(counts_wide, grand_totals) # Merge to get final sample counts 

print(sample_counts)


# Define your four genera of interest
target_genera <- c("Rickettsiella", "Wolbachia", "Rickettsia", "Candidatus Cardinium", "Serratia")

# Initialize a list to store results
strain_consistency_list <- list()

# Loop over each genus
for (genus in target_genera) {
  
  # Find matching taxa
  genus_taxa <- taxa_names(phy_cluster)[grepl(genus, tax_table(phy_cluster)[, "Genus"])]
  
  # Subset phyloseq object
  phy_genus <- prune_taxa(genus_taxa, phy_cluster)
  phy_genus <- prune_samples(sample_sums(phy_genus) > 0, phy_genus)
  
  # Get metadata
  metadata <- as(sample_data(phy_genus), "data.frame")
  metadata$Sample_ID <- rownames(metadata)
  
  # Melt OTU table
  otu_df <- as.data.frame(otu_table(phy_genus))
  otu_df$ASV <- rownames(otu_df)
  long_df <- pivot_longer(otu_df, cols = -ASV, names_to = "Sample_ID", values_to = "Abundance")
  
  # Merge and filter
  long_df <- merge(long_df, metadata, by = "Sample_ID")
  long_df <- long_df[long_df$Abundance > 0, ]
  
  # Summarize
  strain_summary <- long_df %>%
    group_by(Island, Species.Identification) %>%
    summarise(
      Genus = genus,
      Unique_ASVs = n_distinct(ASV),
      Total_Individuals = n_distinct(Sample_ID),
      .groups = "drop"
    )
  
  # Add to list
  strain_consistency_list[[genus]] <- strain_summary
}

# Combine all results into a single data frame
strain_consistency_all <- bind_rows(strain_consistency_list) ## Final product: number of ASVs per host and site combination

# View result
print(strain_consistency_all)

write.csv(strain_consistency_all, "./SERDP/Invasive_Spider_Fecal_Run/data/strain_consistency_all.csv", row.names = TRUE)


species_to_remove <- c("Scytodes", "Thelacantha") #Remove these species due to low abundance across locations 

## Prep for plotting 
strain_summary <- strain_consistency_all %>%
  arrange(Genus, Island, Species.Identification) %>%
  filter(!(Species.Identification %in% species_to_remove)) 

ggplot(strain_summary, aes(x = Species.Identification, y = Unique_ASVs, fill = Genus)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Island, scales = "free_x") +
  labs(
    x = "Spider Species",
    y = "Number of Unique ASVs",
    fill = "Endosymbiont Genus",
    title = "Strain Diversity of Endosymbionts by Spider Species and Island"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############################################### GENERATE RELATIVE ABUNDANCE TABLE ###############################################################################################
# Remove samples with zero counts
ps <- prune_samples(sample_sums(phy_cluster) > 0, phy_cluster)

# Remove taxa with zero counts
ps <- prune_taxa(taxa_sums(phy_cluster) > 0, phy_cluster)

# Function to calculate relative abundance
transform_to_relative_abundance <- function(phy_cluster) {
  physeq_rel <- transform_sample_counts(phy_cluster, function(x) x / sum(x))
  return(physeq_rel)
}

# Apply relative abundance transformation
ps_rel <- transform_to_relative_abundance(phy_cluster)


#For each taxonomic assignment, congregate the taxonomic data
physeq_glom_p <- tax_glom(ps_rel, taxrank = "Phylum")
PhyOTU <- otu_table(physeq_glom_p) #generate OTU table
PhyTAX <- tax_table(physeq_glom_p)[,"Phylum"] #Generate tax table 
clrtrans_p <- microbiome::transform(PhyOTU, 'clr') #Convert to clr to better visualize the data 

#Repeat the above steps for all taxonomic levels 
physeq_glom_c <- tax_glom(ps_rel, taxrank = "Class")
PhyOTUc <- otu_table(physeq_glom_c)
PhyTAXc <- tax_table(physeq_glom_c)[,"Class"]
clrtrans_c <- microbiome::transform(PhyOTUc, 'clr')

physeq_glom_o <- tax_glom(ps_rel, taxrank = "Order")
PhyOTUo <- otu_table(physeq_glom_o)
PhyTAXo <- tax_table(physeq_glom_o)[,"Order"]
clrtrans_o <- microbiome::transform(PhyOTUo, 'clr')

physeq_glom_f <- tax_glom(ps_rel, taxrank = "Family")
PhyOTUf <- otu_table(physeq_glom_f)
PhyTAXf <- tax_table(physeq_glom_f)[,"Family"]
clrtrans_f <- microbiome::transform(PhyOTUf, 'clr')

physeq_glom_g <- tax_glom(ps_rel, taxrank = "Genus")
PhyOTUg <- otu_table(physeq_glom_g)
PhyTAXg <- tax_table(physeq_glom_g)[,"Genus"]
clrtrans_g <- microbiome::transform(PhyOTUg, 'clr')

physeq_glom_s <- tax_glom(ps_rel, taxrank = "Species")
PhyOTUs <- otu_table(physeq_glom_s)
PhyTAXs <- tax_table(physeq_glom_s)[,"Species"]
clrtrans_s <- microbiome::transform(PhyOTUg, 'clr')

#Using these files, create dataframes for each taxonomic level 

# =====================
#PHYLUM
# =====================

#Get the tax tables ready
copy_PhyTAX <- as.data.frame(PhyTAX) #make it a data frame so it's easy to use
copy_PhyTAX <- tibble::rownames_to_column(copy_PhyTAX, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_clr <- as.data.frame(PhyOTU) #make it a data frame so it's easy to use
copy_clr <- tibble::rownames_to_column(copy_clr, "rn")

#Merging the CLR and TAX data
Phy_merge <- merge(copy_clr,copy_PhyTAX, by=c("rn")) #merging data 
Phy_merge2 <- Phy_merge %>%
  dplyr::select(Phylum, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Phy_merge_final <- as.data.frame(t(Phy_merge2)) #make rows columns and columns rows
Phy_merge_final <- tibble::rownames_to_column(Phy_merge_final, "testing") #fixing some names

colnames(Phy_merge_final) <- Phy_merge_final[1, ]; Phy_merge_final <- Phy_merge_final[- 1, ] #make first row the column names
names(Phy_merge_final)[names(Phy_merge_final) == 'Phylum'] <- 'sampleid' #rename a column
View(Phy_merge_final)

# =====================
#CLASS
# =====================

#Get the tax tables ready
copy_PhyTAXc <- as.data.frame(PhyTAXc) #make it a data frame so it's easy to use
copy_PhyTAXc <- tibble::rownames_to_column(copy_PhyTAXc, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_clrc <- as.data.frame(PhyOTUc) #make it a data frame so it's easy to use
copy_clrc <- tibble::rownames_to_column(copy_clrc, "rn")

#Merging the CLR and TAX data
Phy_merge_c <- merge(copy_clrc,copy_PhyTAXc, by=c("rn")) #merging data 
Phy_merge2_c <- Phy_merge_c %>%
  dplyr::select(Class, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Cla_merge_final <- as.data.frame(t(Phy_merge2_c)) #make rows columns and columns rows
Cla_merge_final <- tibble::rownames_to_column(Cla_merge_final, "testing") #fixing some names

colnames(Cla_merge_final) <- Cla_merge_final[1, ]; Cla_merge_final <- Cla_merge_final[- 1, ] #make first row the column names
names(Cla_merge_final)[names(Cla_merge_final) == 'Class'] <- 'sampleid' #rename a column
View(Cla_merge_final)

# =====================
#ORDER
# =====================

#Get the tax tables ready
copy_PhyTAXo <- as.data.frame(PhyTAXo) #make it a data frame so it's easy to use
copy_PhyTAXo <- tibble::rownames_to_column(copy_PhyTAXo, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_otuo <- as.data.frame(PhyOTUo) #make it a data frame so it's easy to use
copy_otuo <- tibble::rownames_to_column(copy_otuo, "rn")

#Merging the CLR and TAX data
Phy_merge_o <- merge(copy_otuo,copy_PhyTAXo, by=c("rn")) #merging data 
Phy_merge2_o <- Phy_merge_o %>%
  dplyr::select(Order, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Or_merge_final <- as.data.frame(t(Phy_merge2_o)) #make rows columns and columns rows
Or_merge_final <- tibble::rownames_to_column(Or_merge_final, "testing") #fixing some names

colnames(Or_merge_final) <- Or_merge_final[1, ]; Or_merge_final <- Or_merge_final[- 1, ] #make first row the column names
names(Or_merge_final)[names(Or_merge_final) == 'Order'] <- 'sampleid' #rename a column
View(Or_merge_final)

# =====================
#FAMILY
# =====================

#Get the tax tables ready
copy_PhyTAXf <- as.data.frame(PhyTAXf) #make it a data frame so it's easy to use
copy_PhyTAXf <- tibble::rownames_to_column(copy_PhyTAXf, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_otuf <- as.data.frame(PhyOTUf) #make it a data frame so it's easy to use
copy_otuf <- tibble::rownames_to_column(copy_otuf, "rn")

#Merging the CLR and TAX data
Phy_merge_f <- merge(copy_otuf,copy_PhyTAXf, by=c("rn")) #merging data 
Phy_merge2_f <- Phy_merge_f %>%
  dplyr::select(Family, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Fam_merge_final <- as.data.frame(t(Phy_merge2_f)) #make rows columns and columns rows
Fam_merge_final <- tibble::rownames_to_column(Fam_merge_final, "testing") #fixing some names

colnames(Fam_merge_final) <- Fam_merge_final[1, ]; Fam_merge_final <- Fam_merge_final[- 1, ] #make first row the column names
names(Fam_merge_final)[names(Fam_merge_final) == 'Family'] <- 'sampleid' #rename a column
View(Fam_merge_final)

# =====================
#GENUS
# =====================

#Get the tax tables ready
copy_PhyTAXg <- as.data.frame(PhyTAXg) #make it a data frame so it's easy to use
copy_PhyTAXg <- tibble::rownames_to_column(copy_PhyTAXg, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_otug <- as.data.frame(PhyOTUg) #make it a data frame so it's easy to use
copy_otug <- tibble::rownames_to_column(copy_otug, "rn")

#Merging the CLR and TAX data
Phy_merge_g <- merge(copy_otug,copy_PhyTAXg, by=c("rn")) #merging data 
Phy_merge2_g <- Phy_merge_g %>%
  dplyr::select(Genus, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Gen_merge_final <- as.data.frame(t(Phy_merge2_g)) #make rows columns and columns rows
Gen_merge_final <- tibble::rownames_to_column(Gen_merge_final, "testing") #fixing some names

colnames(Gen_merge_final) <- Gen_merge_final[1, ]; Gen_merge_final <- Gen_merge_final[- 1, ] #make first row the column names
names(Gen_merge_final)[names(Gen_merge_final) == 'Genus'] <- 'sampleid' #rename a column
View(Gen_merge_final)

# =====================
#SPECIES
# =====================

#Get the tax tables ready
copy_PhyTAXs <- as.data.frame(PhyTAXs) #make it a data frame so it's easy to use
copy_PhyTAXs <- tibble::rownames_to_column(copy_PhyTAXs, "rn") #make row names a column for ease of use

#FOR PERCENTAGES: get OTU table ready 
copy_otus <- as.data.frame(PhyOTUs) #make it a data frame so it's easy to use
copy_otus <- tibble::rownames_to_column(copy_otus, "rn")

#Merging the CLR and TAX data
Phy_merge_s <- merge(copy_otus,copy_PhyTAXs, by=c("rn")) #merging data 
Phy_merge2_s <- Phy_merge_s %>%
  dplyr::select(Species, everything()) %>%  # bring Phylum to front
  dplyr::select(-rn)                      # drop the rn column

Sps_merge_final <- as.data.frame(t(Phy_merge2_s)) #make rows columns and columns rows
Sps_merge_final <- tibble::rownames_to_column(Sps_merge_final, "testing") #fixing some names

colnames(Sps_merge_final) <- Sps_merge_final[1, ]; Sps_merge_final <- Sps_merge_final[- 1, ] #make first row the column names
names(Sps_merge_final)[names(Sps_merge_final) == 'Species'] <- 'sampleid' #rename a column
View(Sps_merge_final)

# =====================
###MERGING 
# =====================

final_final <- merge(Phy_merge_final, Cla_merge_final, by=c("sampleid")) #do this the first time, then:
final_final <- merge(final_final, Fam_merge_final, by=c("sampleid")) #do this the rest of the times (copy this row as many times as you need)
final_final <- merge(final_final, Or_merge_final, by=c("sampleid")) #do this the rest of the times (copy this row as many times as you need)
final_final <- merge(final_final, Gen_merge_final, by=c("sampleid")) #do this the rest of the times (copy this row as many times as you need)
final_final <- merge(final_final, Sps_merge_final, by=c("sampleid")) #do this the rest of the times (copy this row as many times as you need)

View(final_final) #This is the final object

final_final <- as.data.frame(final_final) #Make dataframe

write.csv(final_final, "./SERDP/Invasive_Spider_Fecal_Run/data/final_final_18Jul25.csv", row.names = TRUE) ## Save as file



################################### ANALYSIS ###########################################

##### SUBSETTING FOR ENDOSYMBIONTS FOR SIGNIFICANCE TESTING ######################################

#Remove columns of known endosymbionts based on literature and assume gut microbiota for the rest 

final_final <- read.csv("./SERDP/Invasive_Spider_Fecal_Run/data/final_final_18Jul25.csv")

subset_endo <- final_final %>%
  dplyr::select("sampleid", "g__Wolbachia","g__Candidatus.Cardinium","g__Rickettsia","g__Rickettsiella") #Subset based on endosymbionts

names(subset_endo)[names(subset_endo) == 'sampleid'] <- 'Sample_ID'

metadata_merge <- merge(subset_endo, metadata_16s, by="Sample_ID")

no_blank_merged <- metadata_merge[-c(1,2,3), ] #Remove blanks from dataframe

subset_microbiota <- final_final %>%
  dplyr::select (-"g__Wolbachia",-"g__Candidatus.Cardinium",-"g__Rickettsia",-"g__Rickettsiella") #Remove known endosymbionts to obtain file with remaining gut microbiota

names(subset_microbiota)[names(subset_microbiota) == 'sampleid'] <- 'Sample_ID'

gut_microbe_merge <- merge(subset_microbiota, metadata_16s, by="Sample_ID")

no_blanks <- gut_microbe_merge[-c(1,2,3), ] #Remove blanks form dataframe

###### STATISTICAL TESTING FOR KNOWN ENDOSYMBIONTS #######################################################################################

ggplot(no_blank_merged, aes(x=g__Wolbachia)) + 
  geom_histogram() ##Visualize Wolbachia abundance

# Notes on Wolbachia presence:
# Wolbachia is only found in the genera Tetragnatha, Scytodes, and Steatoda,
# and only in the Berkeley and Oahu islands.

# Remove species with low Wolbachia abundance across locations to avoid bias in models
species_to_remove <- c("Scytodes", "Thelacantha") #Remove these species due to low abundance across locations 
no_blank_merged <- no_blank_merged %>%
  filter(!(Species.Identification %in% species_to_remove)) 


###Species-specific ###
# create presence/absence columns
no_blank_merged <- no_blank_merged %>%
  mutate(
    Wolbachia_present    = ifelse(g__Wolbachia >  0, 1, 0),
    Cardinium_present   = ifelse(g__Candidatus.Cardinium > 0, 1, 0),
    Rickettsia_present  = ifelse(g__Rickettsia >   0, 1, 0),
    Rickettsiella_present= ifelse(g__Rickettsiella >0, 1, 0)
  )

# 1. Wolbachia ~ Island + Specimen Type
m_wol_simple <- glmmTMB(
  Wolbachia_present ~ Island + Type.of.Specimen,
  data   = no_blank_merged,
  family = binomial
)
summary(m_wol_simple)
diagnostics_wol <- simulateResiduals(m_wol_simple); plot(diagnostics_wol)

# Wolbachia ~ Species
m_wol_sp <- glmmTMB(
  Wolbachia_present ~ Species.Identification,
  data   = no_blank_merged,
  family = binomial
)
summary(m_wol_sp)
diagnostics_wol_sp <- simulateResiduals(m_wol_sp); plot(diagnostics_wol_sp)

# 2. Cardinium ~ Island + Specimen Type
m_cand_simple <- glmmTMB(
  Cardinium_present ~ Island + Type.of.Specimen,
  data   = no_blank_merged,
  family = binomial
)
summary(m_cand_simple)
diagnostics_cand <- simulateResiduals(m_cand_simple); plot(diagnostics_cand)

# Cardinium ~ Species
m_cand_sp <- glmmTMB(
  Cardinium_present ~ Species.Identification,
  data   = no_blank_merged,
  family = binomial
)
summary(m_cand_sp)
diagnostics_cand_sp <- simulateResiduals(m_cand_sp); plot(diagnostics_cand_sp)

# 3. Rickettsia ~ Island + Specimen Type
m_rik_simple <- glmmTMB(
  Rickettsia_present ~ Island + Type.of.Specimen,
  data   = no_blank_merged,
  family = binomial
)
summary(m_rik_simple)
diagnostics_rik <- simulateResiduals(m_rik_simple); plot(diagnostics_rik)

# Rickettsia ~ Species
m_rik_sp <- glmmTMB(
  Rickettsia_present ~ Species.Identification,
  data   = no_blank_merged,
  family = binomial
)
summary(m_rik_sp)
diagnostics_rik_sp <- simulateResiduals(m_rik_sp); plot(diagnostics_rik_sp)

# 4. Rickettsiella ~ Island + Specimen Type
m_rikella_simple <- glmmTMB(
  Rickettsiella_present ~ Island + Type.of.Specimen,
  data   = no_blank_merged,
  family = binomial
)
summary(m_rikella_simple)
diagnostics_rikella <- simulateResiduals(m_rikella_simple); plot(diagnostics_rikella)

# Rickettsiella ~ Species
tab <- table(no_blank_merged$Species.Identification,
             no_blank_merged$Rickettsiella_present)

# 2. Identify which species have both 0s and 1s
keep_species <- rownames(tab)[rowSums(tab > 0) == 2]

# 3. Subset
df_keep <- subset(no_blank_merged,
                  Species.Identification %in% keep_species)

# 4. Re‐fit
m_rikella_sp2 <- glmmTMB(
  Rickettsiella_present ~ Species.Identification,
  data   = df_keep,
  family = binomial
)
summary(m_rikella_sp2)
diagnostics_rikella_sp <- simulateResiduals(m_rikella_sp2); plot(diagnostics_rikella_sp)

# ============================================
##### INFECTION RATES FOR ENDOSYMBIONTS ###### 
# ============================================

# Create binary presence/absence again just in case
no_blank_merged$Wolbachia_present <- ifelse(no_blank_merged$g__Wolbachia > 0, 1, 0)
no_blank_merged$Cardinium_present <- ifelse(no_blank_merged$g__Candidatus.Cardinium > 0, 1, 0)
no_blank_merged$Rickettsia_present <- ifelse(no_blank_merged$g__Rickettsia > 0, 1, 0)
no_blank_merged$Rickettsiella_present <- ifelse(no_blank_merged$g__Rickettsiella > 0, 1, 0)

# Function to calculate infection rate per Island + Species combo
calc_infection_rate <- function(data, var){
  data %>%
    group_by(Island, Species.Identification) %>%
    summarise(InfectionRate = mean(.data[[var]]), 
              N = n()) %>%
    mutate(Symbiont = var)
}

# Apply for each endosymbiont
wol_rate <- calc_infection_rate(no_blank_merged, "Wolbachia_present")
card_rate <- calc_infection_rate(no_blank_merged, "Cardinium_present")
rick_rate <- calc_infection_rate(no_blank_merged, "Rickettsia_present")
rikella_rate <- calc_infection_rate(no_blank_merged, "Rickettsiella_present")

# Combine into one dataframe
infection_rates <- bind_rows(wol_rate, card_rate, rick_rate, rikella_rate)

# Plot infection rate
ggplot(infection_rates, aes(x = Island, y = InfectionRate, fill = Symbiont)) +
  geom_col(position = "dodge") +
  facet_wrap(~Species.Identification, scales = "free_y") +
  labs(title = "Infection Rates per Island and Species",
       y = "Proportion Infected", x = "Island") +
  theme_minimal()

# Round infection rates and reshape for table formatting
infection_table <- infection_rates %>%
  mutate(InfectionRate = round(InfectionRate, 2)) %>%
  tidyr::pivot_wider(names_from = Symbiont, values_from = InfectionRate) %>%
  arrange(Species.Identification, Island)

# Optional: Replace NA with 0s or dashes if you prefer
infection_table[is.na(infection_table)] <- 0

# View table
print(infection_table)

write.csv(infection_table, "./SERDP/Invasive_Spider_Fecal_Run/data/infection_table.csv", row.names = TRUE)


###### GUT MICROBE VS. ENDOSYMBIONTS COMPOSITION PLOTS ####################################

# =====================
### COMPOSITION PLOTS
# =====================

# Prune samples and taxa with zero counts
ps <- prune_samples(sample_sums(phy_cluster) > 0, phy_cluster)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Remove blank samples
samples_to_remove <- c("BLANK_16S_1", "BLANK_16S_2", "PCR_BLANK_16S")
ps <- prune_samples(!(sample_names(ps) %in% samples_to_remove), ps)

# Transform counts to relative abundance
transform_to_relative_abundance <- function(physeq_obj) {
  transform_sample_counts(physeq_obj, function(x) x / sum(x))
}
ps_rel <- transform_to_relative_abundance(ps)

# =====================
### ENDOSYMBIONT SUBSETTING ###
# =====================

tax_df <- as.data.frame(tax_table(ps_rel))
genera_endosymbionts <- c("g__Wolbachia", "g__Candidatus Cardinium", "g__Rickettsia", "g__Rickettsiella", "g__Serratia")

# Note: Remove leading spaces if any from genera names
tax_df$Genus <- trimws(tax_df$Genus)

ps_endosymbionts <- prune_taxa(tax_df$Genus %in% genera_endosymbionts, ps_rel)

# =====================
### GUT MICROBIOTA SUBSETTING ###
# =====================

# Define gut microbes as everything NOT in endosymbionts plus some additional genera
genera_additional <- c("g__Brachybacterium", "g__Methylobacterium")
genera_remove <- c(genera_endosymbionts, genera_additional)

ps_gut <- prune_taxa(!(tax_df$Genus %in% genera_remove), ps_rel)

# Flag unidentified genera in gut microbiota subset
tax_gut_df <- as.data.frame(tax_table(ps_gut))
tax_gut_df$Genus <- trimws(tax_gut_df$Genus)
unidentified_genera <- grepl("^g__$", tax_gut_df$Genus)

# =====================
### ALPHA DIVERSITY ANALYSIS ON GUT MICROBIOTA ###
# =====================

# Prune taxa from original dataset (not relative abundance) for gut microbes
og_gut_only <- prune_taxa(!(tax_df$Genus %in% genera_remove), phy_cluster)

# Remove unwanted species from the sample data
species_to_remove <- c("Scytodes", "Thelacantha")
phy_cluster_filtered <- subset_samples(phy_cluster, !(Species.Identification %in% species_to_remove))

# Calculate Shannon diversity
alpha_div <- estimate_richness(phy_cluster_filtered, measures = "Shannon")

# Merge with sample metadata
sample_metadata <- data.frame(sample_data(phy_cluster_filtered))
alpha_div_merged <- cbind(sample_metadata, alpha_div)

# Remove blank or control samples by index if known
# (Assuming rows 1,2,3 are blanks here — adjust if needed)
alpha_div_clean <- alpha_div_merged[-c(1, 2, 3), ]

# Create interaction factor for plotting and analysis
alpha_div_clean <- alpha_div_clean %>%
  mutate(Site_Species = interaction(Island, Species.Identification, sep = "."))

# Define color mapping for Site x Species combinations
color_mapping <- c(
  "Berkeley.Badumna longinqua"       = "#9b4fff",
  "Berkeley.Tetragnatha versicolor"  = "#9e9ac8",
  "Hawaii (Big Island).Badumna longinqua"         = "#1b7837",
  "Hawaii (Big Island).Steatoda grossa"           = "green",
  "Hawaii (Big Island).Tetragnatha acuta"         = "#a6dba0",
  "Maui.Steatoda grossa"             = "aquamarine",
  "Maui.Tetragnatha eurychasma"      = "cyan2",
  "Oahu.Steatoda grossa"             = "#a83b51",
  "Oahu.Tetragnatha spp."            = "#c28d98"
  
)

# Verify all Site_Species have assigned colors
missing_colors <- setdiff(unique(alpha_div_clean$Site_Species), names(color_mapping))
if (length(missing_colors) > 0) {
  stop(paste("Missing colors for:", paste(missing_colors, collapse = ", ")))
}

# =====================
### ANOVA TESTS ###
# =====================

# Use Type III sum of squares for unbalanced data
alpha_div_clean$Island <- factor(alpha_div_clean$Island)
alpha_div_clean$Species.Identification <- factor(alpha_div_clean$Species.Identification)
alpha_div_clean$Type.of.Specimen <- factor(alpha_div_clean$Type.of.Specimen)

## 1) Main additive model
alpha_model_main <- lm(Shannon ~ Island + Species.Identification + Type.of.Specimen,
                       data = alpha_div_clean)

# Diagnostics
par(mfrow=c(2,2))
plot(alpha_model_main)
par(mfrow=c(1,1))

# Type III ANOVA
Anova(alpha_model_main, type="III")

# DHARMa residual check
res_main <- simulateResiduals(fittedModel = alpha_model_main, plot = TRUE)


## 3) Interaction: Island × Type.of.Specimen
alpha_model_it <- lm(Shannon ~ Island * Type.of.Specimen,
                     data = alpha_div_clean)

# Type III ANOVA
Anova(alpha_model_it, type="III")

# Residual check
res_it <- simulateResiduals(fittedModel = alpha_model_it, plot = TRUE)

# 3. For Species‐within‐Island differences, run separate one‐way ANOVAs on each island
for(isl in unique(alpha_div_clean$Island)){
  df_i <- subset(alpha_div_clean, Island==isl)
  if(length(unique(df_i$Species.Identification)) > 1){
    cat("\nIsland:", isl, "\n")
    print(Anova(aov(Shannon ~ Species.Identification, data=df_i), type="II"))
  }
}

emm <- emmeans(alpha_model_main, ~ Island * Species.Identification)

# 4) get pairwise Tukey comparisons
pw <- contrast(emm, method = "pairwise", adjust = "tukey")

# 5) pull out just the significant ones
pw_df <- as.data.frame(summary(pw))
sig <- pw_df[pw_df$p.value < 0.05, ]

# 6) build the list of comparisons for geom_signif
#    note that emmeans prints contrasts as "A B - C D"
comparisons <- strsplit(sig$contrast, " - ")


alpha_div_clean <- alpha_div_clean %>%
  dplyr::mutate(
    Site_Species = interaction(Island, Species.Identification, sep = "."),
    Site_Species = factor(Site_Species, levels = names(color_mapping))
  )

# =====================
### PLOT ###
# =====================

ggplot(alpha_div_clean, aes(x = Site_Species, y = Shannon, fill = Site_Species)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  coord_flip() +
  scale_fill_manual(values = color_mapping, guide = "none") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = "Site × Species",
    y = "Shannon–Wiener Diversity Index"
  )



# Clean genus names by trimming leading/trailing spaces
tax_table_df <- as.data.frame(tax_table(ps_rel))
tax_table_df$Genus <- trimws(tax_table_df$Genus)

# Define endosymbiont genera and additional genera for exclusion
endosymbiont_genera <- c("g__Wolbachia", "g__Candidatus Cardinium", "g__Rickettsia", "g__Rickettsiella", "g__Serratia")
additional_genera <- c("g__Brachybacterium", "g__Methylobacterium")

# =====================
# GUT MICROBIOTA SUBSET
# =====================

# Filter taxa NOT in endosymbionts or additional genera to get gut microbes
genera_to_exclude <- c(endosymbiont_genera, additional_genera)
gut_taxa <- tax_table_df$Genus %in% genera_to_exclude
ps_gut <- prune_taxa(!gut_taxa, ps_rel)

# Remove unidentified genera labeled as "g__"
tax_table_gut <- as.data.frame(tax_table(ps_gut))
tax_table_gut$Genus <- trimws(tax_table_gut$Genus)
unidentified_genera <- tax_table_gut$Genus == "g__"
ps_gut_clean <- prune_taxa(!unidentified_genera, ps_gut)

# =====================
# ENDOSYMBIONT SUBSET
# =====================

ps_endosymbiont <- prune_taxa(tax_table_df$Genus %in% endosymbiont_genera, ps_rel)

# =====================
# REMOVE UNWANTED SPECIES
# =====================

species_to_remove <- c("Scytodes", "Thelacantha")

ps_gut_filtered <- subset_samples(ps_gut_clean, !Species.Identification %in% species_to_remove)
ps_endosymbiont_filtered <- subset_samples(ps_endosymbiont, !Species.Identification %in% species_to_remove)

# =====================
# AGGREGATE TO GENUS LEVEL
# =====================

ps_gut_genus <- tax_glom(ps_gut_filtered, taxrank = "Genus")
ps_endosymbiont_genus <- tax_glom(ps_endosymbiont_filtered, taxrank = "Genus")

# =====================
# CONVERT TO DATA FRAMES FOR PLOTTING
# =====================

gut_df <- psmelt(ps_gut_genus)
endosymbiont_df <- psmelt(ps_endosymbiont_genus)

# Replace "g__" with "Not Assigned" for clarity
gut_df$Genus <- ifelse(gut_df$Genus == "g__", "Not Assigned", gut_df$Genus)

# =====================
# GUT MICROBIOTA: TOP 10 GENERA PLOT
# =====================

# Identify top 10 genera by total abundance
top_genera <- gut_df %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(dplyr::desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Filter gut data to top 10 genera
gut_df_top <- gut_df %>%
  filter(Genus %in% top_genera)

# Calculate relative abundance within each Island and Species
gut_df_top <- gut_df_top %>%
  group_by(Island, Species.Identification) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# Create color palette including white for "Not Assigned"
unique_genera <- unique(gut_df_top$Genus)
palette_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(length(unique_genera))
names(palette_colors) <- unique_genera
if ("Not Assigned" %in% unique_genera) {
  palette_colors["Not Assigned"] <- "white"
}

clean_labels <- gsub("^g__", "", names(palette_colors_other))
clean_labels[clean_labels==""] <- "Not Assigned"

# Plot gut microbiota top 10 genera
ggplot(gut_df_top, aes(x = as.factor(Species.Identification), y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_colors, labels = clean_labels) +
  theme_classic() +
  labs(x = "Spider Genus", y = "Relative Abundance (%)") +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +
  theme(
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15),
    strip.background = element_rect(color = "white", fill = "white"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text = element_text(size = 13)
  )

# =====================
# GUT MICROBIOTA: TOP 10 + OTHER CATEGORY PLOT
# =====================
top_genera15 <- gut_df %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# 2) Lump everything else as “Other” and compute % abundance
gut_df_extended <- gut_df %>%
  mutate(
    Genus = ifelse(Genus %in% top_genera15, Genus, "Other")
  ) %>%
  group_by(Island, Species.Identification, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Island, Species.Identification) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# 2b) Make sure "Other" is the last factor level
gut_df_extended$Genus <- factor(
  gut_df_extended$Genus,
  levels = c(top_genera15, "Other")
)

palette_extended <- c(
  setNames(
    colorRampPalette(brewer.pal(9, "YlGnBu"))(length(top_genera15)),
    top_genera15
  ),
  Other = "grey80"
)

# 4) Build a named vector of cleaned labels
labels_map <- setNames(
  gsub("^g__", "", names(palette_extended)),  # strip the "g__"
  names(palette_extended)                     # keep the same names
)
labels_map[labels_map == ""] <- "Not Assigned"  # fix any blanks

# 5) Plot, using named labels_map
ggplot(gut_df_extended, aes(Species.Identification, RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = palette_extended,
    breaks = names(palette_extended),   # ensure all levels show
    labels = labels_map                 # map each original name → cleaned label
  ) +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(x = "Spider Genus", y = "Relative Abundance (%)") +
  theme(
    axis.title       = element_text(size = 20),
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 13),
    strip.text       = element_text(size = 15),
    strip.background = element_rect(color = "white", fill = "white"),
    axis.text.x      = element_text(angle = 90, hjust = 1),
    axis.text        = element_text(size = 13)
  )

# =====================
# ENDOSYMBIONT: RELATIVE ABUNDANCE PLOT
# =====================

# Clean up endosymbiont Genus column
endosymbiont_df$Genus <- trimws(endosymbiont_df$Genus)
endosymbiont_df$Genus <- sub("^g__", "", endosymbiont_df$Genus)

# Calculate relative abundance within each Island × Species group
endo_df_rel <- endosymbiont_df %>%
  group_by(Island, Species.Identification) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# Define your exact hex palette
symbiont_colors <- c(
  "Wolbachia"              = "#eeba0b",    # yellow
  "Rickettsiella"          = "#99d6ea",    # steelblue3 hex
  "Rickettsia"             = "#8B0000",    # darkred hex
  "Candidatus Cardinium"   = "#ff8fa3",    # hotpink hex
  "Serratia"               = "#800080"     # purple hex
)

#Plot endosymbiont abundance plot 
ggplot(endo_df_rel, 
       aes(x = as.factor(Species.Identification), 
           y = RelativeAbundance, 
           fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = symbiont_colors) +
  theme_classic() +
  labs(
    x = "Spider Species",
    y = "Endosymbiont Relative Abundance (%)",
    fill = "Endosymbiont Genus"
  ) +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +
  theme(
    axis.title       = element_text(size = 20),
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 13),
    strip.text       = element_text(size = 15),
    strip.background = element_rect(color = "white", fill = "white"),
    axis.text.x      = element_text(angle = 90, hjust = 1),
    axis.text        = element_text(size = 13)
  )

###### ANALYSIS FOR GUT MICROBIOTA ###################################################
# Remove blank and PCR control samples from the combined phyloseq object
samples_to_remove <- c("BLANK_16S_1", "BLANK_16S_2", "PCR_BLANK_16S")
phy_cluster_pruned <- prune_samples(!(sample_names(phy_cluster) %in% samples_to_remove), phy_cluster)

# Calculate Bray-Curtis distance matrix for all gut samples
combined_distance_total <- phyloseq::distance(phy_cluster_pruned, method = "bray")

# Extract and format metadata
perm_16s_metadata <- as.data.frame(as.matrix(sample_data(phy_cluster_pruned)))
perm_16s_metadata$Island <- as.factor(perm_16s_metadata$Island)
perm_16s_metadata$Species.Identification <- as.factor(perm_16s_metadata$Species.Identification)
perm_16s_metadata$Type.of.Specimen <- as.factor(perm_16s_metadata$Type.of.Specimen)

# Check input classes before PERMANOVA
print(class(combined_distance_total))
print(class(perm_16s_metadata))

# Run PERMANOVA for full dataset
permanova_all <- vegan::adonis2(
  combined_distance_total ~ Island + Species.Identification + Type.of.Specimen,
  data = perm_16s_metadata,
  by = "terms"
)
print(permanova_all)
disp <- betadisper(permanova_all, perm_16s_metadata$Species.Identification)

### PERMANOVA for gut microbiota only
beta_dist_gut <- phyloseq::distance(ps_gut_filtered, method = "bray")
metadata_gut <- as.data.frame(as.matrix(sample_data(ps_gut_filtered)))
disp <- betadisper(beta_dist_gut, metadata_gut)
anova(disp)
plot(disp)

# 1) ANOSIM for each main factor
anosim_island  <- anosim(beta_dist_gut, metadata_gut$Island, permutations = 999)
anosim_species <- anosim(beta_dist_gut, metadata_gut$Species.Identification, permutations = 999)
anosim_type    <- anosim(beta_dist_gut, metadata_gut$Type.of.Specimen, permutations = 999)

# Print summaries
cat("\n--- ANOSIM: Island ---\n");  print(anosim_island)
cat("\n--- ANOSIM: Species ---\n"); print(anosim_species)
cat("\n--- ANOSIM: Type.of.Specimen ---\n"); print(anosim_type)


# 2) Re‑run pairwise PERMANOVA for species and for islands
pairwise_species <- pairwise.adonis2(
  beta_dist_gut ~ Species.Identification,
  data = metadata_gut,
  permutations = 999,
  by = "terms"
)

pairwise_island <- pairwise.adonis2(
  beta_dist_gut ~ Island,
  data = metadata_gut,
  permutations = 999,
  by = "terms"
)

cat("\n--- Pairwise PERMANOVA: Species ---\n");  print(pairwise_species)
cat("\n--- Pairwise PERMANOVA: Island ---\n");   print(pairwise_island)

# function to run pairwise dispersion checks
check_pairwise_dispersion <- function(dist_mat, grouping, factor_name) {
  levs <- unique(grouping)
  for(i in seq_along(levs)[-length(levs)]) {
    for(j in (i+1):length(levs)) {
      pair <- levs[c(i,j)]
      sel   <- grouping %in% pair
      # subset the distance matrix
      sub_mat   <- as.dist(as.matrix(dist_mat)[sel, sel])
      sub_group <- grouping[sel]
      disp      <- betadisper(sub_mat, sub_group)
      cat("\nDispersion ANOVA for", factor_name, ":", pair[1], "vs", pair[2], "\n")
      print(anova(disp))
    }
  }
}

# 1) Check dispersion for Species.Identification
check_pairwise_dispersion(
  dist_mat = beta_dist_gut,
  grouping = metadata_gut$Species.Identification,
  factor_name = "Species"
)

# 2) Check dispersion for Island
check_pairwise_dispersion(
  dist_mat = beta_dist_gut,
  grouping = metadata_gut$Island,
  factor_name = "Island"
)


# 3) ANOSIM pairwise for a sanity check
# You can also do simple pairwise ANOSIM via vegan::anosim in a loop:
for(fac in c("Island","Species.Identification")){
  cat("\n--- Pairwise ANOSIM for", fac, "---\n")
  lvls <- unique(metadata_gut[[fac]])
  combs <- combn(lvls, 2, simplify = FALSE)
  for(pair in combs){
    sel <- metadata_gut[[fac]] %in% pair
    dm  <- as.dist(as.matrix(beta_dist_gut)[sel, sel])
    md  <- factor(metadata_gut[sel, fac], levels = pair)
    res <- anosim(dm, md, permutations = 999)
    cat(pair[1], "vs", pair[2],
        ": R =", round(res$statistic, 3),
        ", p =", res$signif, "\n")
  }
}

##############################################
# SIMPER Analysis for Gut Microbiota
##############################################

# Select top 100 most abundant ASVs
top_taxa <- names(sort(taxa_sums(ps_gut_filtered), decreasing = TRUE)[1:100])
phy_gut_subset <- prune_taxa(top_taxa, ps_gut_filtered)

# Subset metadata to match samples in the pruned phyloseq object
metadata_gut <- metadata_gut[sample_names(phy_gut_subset), ]
metadata_gut$Island <- as.factor(metadata_gut$Island)
metadata_gut$Species.Identification <- as.factor(metadata_gut$Species.Identification)

# Transpose OTU table for SIMPER analysis
otu_table_matrix <- as.matrix(otu_table(phy_gut_subset))
otu_table_transposed <- t(otu_table_matrix)

# Ensure metadata and OTU matrix rownames match
metadata_gut <- metadata_gut[rownames(metadata_gut) %in% rownames(otu_table_transposed), ]
stopifnot(identical(rownames(otu_table_transposed), rownames(metadata_gut)))  # Should return TRUE

# Run SIMPER for Island and Species
simper_result <- vegan::simper(otu_table_transposed, metadata_gut$Island, permutations = 999)
simper_result_species <- vegan::simper(otu_table_transposed, metadata_gut$Species.Identification, permutations = 999)

# Summarize SIMPER results
simper_summary <- summary(simper_result)
simper_species_sum <- summary(simper_result_species)

# Convert SIMPER results to data frames
simper_df_list <- lapply(simper_summary, as.data.frame)
simper_sps_list <- lapply(simper_species_sum, as.data.frame)

# Combine all group comparisons into single dataframes
combined_simper_df <- do.call(rbind, simper_df_list)
combined_simper_sps <- do.call(rbind, simper_sps_list)

# Add ASV names as a column
combined_simper_df$ASV <- rownames(combined_simper_df)
combined_simper_sps$ASV <- rownames(combined_simper_sps)

##############################################
# Merge SIMPER results with Taxonomy
##############################################

# Extract and convert taxonomy table
taxonomy_table <- tax_table(phy_gut_subset)
taxonomy_df <- as.data.frame(taxonomy_table)
taxonomy_df$ASV <- rownames(taxonomy_df)

# Merge taxonomy info with SIMPER results
merged_results <- merge(combined_simper_df, taxonomy_df, by = "ASV", all.x = TRUE)
merged_sps_simper <- merge(combined_simper_sps, taxonomy_df, by = "ASV", all.x = TRUE)

# View merged outputs
head(merged_results)
head(merged_sps_simper)

# Save merged SIMPER + taxonomy tables
write.csv(merged_results, "./SERDP/Invasive_Spider_Fecal_Run/data/merged_results_simper.csv", row.names = FALSE)
write.csv(merged_sps_simper, "./SERDP/Invasive_Spider_Fecal_Run/data/merged_sps_simper.csv", row.names = FALSE)



###### VENN DIAGRAM FOR OVERLAPPING ASVS ######################################################################

# Check the variable used to group samples by species
sample_data(phy_cluster)$Species.Identification 

# Remove samples from species not of interest (Scytodes, Thelacantha, NA entries)
species_to_remove <- c("Scytodes", "Thelacantha", "NA")
venn_filtered <- subset_samples(phy_cluster, !Species.Identification %in% species_to_remove)

# Remove negative control and blank samples
samples_to_remove <- c("BLANK_16S_1", "BLANK_16S_2", "PCR_BLANK_16S")
venn_filtered <- subset_samples(phy_cluster, !Sample_ID %in% samples_to_remove)

# Clean and standardize the species identification variable
sample_data(venn_filtered)$Species.Identification <- as.character(sample_data(venn_filtered)$Species.Identification)
sample_data(venn_filtered)$Species.Identification <- trimws(sample_data(venn_filtered)$Species.Identification)

# Extract unique species levels to guide Venn diagram construction
species_levels <- unique(sample_data(venn_filtered)$Species.Identification)

# Initialize a list to store ASVs for each species
asv_list <- list()

# Subset samples for each species and extract non-zero ASVs
# 1. Badumna
ps_species_Badumna <- subset_samples(venn_filtered, Species.Identification == "Badumna")
if (length(sample_names(ps_species_Badumna)) > 0) {
  asv_Badumna <- taxa_names(prune_taxa(taxa_sums(ps_species_Badumna) > 0, ps_species_Badumna))
  asv_list[["Badumna"]] <- asv_Badumna
} else {
  print("No samples found for Badumna")
}

# 2. Steatoda
ps_species_SpiderB <- subset_samples(venn_filtered, Species.Identification == "Steatoda")
if (length(sample_names(ps_species_SpiderB)) > 0) {
  asv_SpiderB <- taxa_names(prune_taxa(taxa_sums(ps_species_SpiderB) > 0, ps_species_SpiderB))
  asv_list[["Steatoda"]] <- asv_SpiderB
} else {
  print("No samples found for Steatoda")
}

# 3. Tetragnatha
ps_species_SpiderC <- subset_samples(venn_filtered, Species.Identification == "Tetragnatha")
if (length(sample_names(ps_species_SpiderC)) > 0) {
  asv_SpiderC <- taxa_names(prune_taxa(taxa_sums(ps_species_SpiderC) > 0, ps_species_SpiderC))
  asv_list[["Tetragnatha"]] <- asv_SpiderC
} else {
  print("No samples found for Tetragnatha")
}

# Print the list of ASVs for each species
print(asv_list)

# Create and customize a Venn diagram showing ASV overlap between species
venn.plot <- venn.diagram(
  x = asv_list,
  category.names = names(asv_list),
  filename = NULL,
  output = TRUE,
  fill = c("#40916c", "#90e0ef", "#023e8a"),  # Fill colors for each circle
  alpha = 0.5,                                # Circle transparency
  cex = 2,                                    # Font size for numbers
  cat.cex = 2,                                # Font size for category labels
  fontfamily = "sans",                        # General font
  cat.fontface = "italic",                    # Italicize category labels
  cat.fontfamily = "sans",
  lwd = 2                                     # Line width for circle borders
)

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)

# -----------------------------------------------------
# Repeat the same steps above, this time for **Island**
# -----------------------------------------------------

# Clean and standardize the Island variable
sample_data(venn_filtered)$Island <- as.character(sample_data(venn_filtered)$Island)
sample_data(venn_filtered)$Island <- trimws(sample_data(venn_filtered)$Island)

# Extract island names to guide Venn creation
island_levels <- unique(sample_data(venn_filtered)$Island)

# Initialize a list to store ASVs for each island
asv_list_isl <- list()

# 1. Oahu
ps_is_oahu <- subset_samples(venn_filtered, Island == "Oahu")
if (length(sample_names(ps_is_oahu)) > 0) {
  asv_oahu <- taxa_names(prune_taxa(taxa_sums(ps_is_oahu) > 0, ps_is_oahu))
  asv_list_isl[["Oahu"]] <- asv_oahu
} else {
  print("No samples found for Oahu")
}

# 2. Maui
ps_is_maui <- subset_samples(venn_filtered, Island == "Maui")
if (length(sample_names(ps_is_maui)) > 0) {
  asv_maui <- taxa_names(prune_taxa(taxa_sums(ps_is_maui) > 0, ps_is_maui))
  asv_list_isl[["Maui"]] <- asv_maui
} else {
  print("No samples found for Maui")
}

# 3. Hawaii (Big Island)
ps_is_hi <- subset_samples(venn_filtered, Island == "Hawaii (Big Island)")
if (length(sample_names(ps_is_hi)) > 0) {
  asv_hi <- taxa_names(prune_taxa(taxa_sums(ps_is_hi) > 0, ps_is_hi))
  asv_list_isl[["Hawaii (Big Island)"]] <- asv_hi
} else {
  print("No samples found for Hawaii (Big Island)")
}

# 4. Berkeley
ps_is_berk <- subset_samples(venn_filtered, Island == "Berkeley")
if (length(sample_names(ps_is_berk)) > 0) {
  asv_berk <- taxa_names(prune_taxa(taxa_sums(ps_is_berk) > 0, ps_is_berk))
  asv_list_isl[["Berkeley"]] <- asv_berk
} else {
  print("No samples found for Berkeley")
}

# Print the ASV list for each island
print(asv_list_isl)

# Create and customize a Venn diagram showing ASV overlap between islands
venn.plot <- venn.diagram(
  x = asv_list_isl,
  category.names = names(asv_list_isl),
  filename = NULL,
  output = TRUE,
  fill = c("#f7b538", "#780116", "#db7c26", "#f25c54"),  # Fill colors
  alpha = 0.5,                                # Transparency
  cat.cex = 2,                                # Label font size
  cex = 2,                                    # Number font size
  fontfamily = "sans",
  cat.fontface = "italic",                    # Italicized labels
  cat.fontfamily = "sans",
  lwd = 2                                     # Line width
)

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)