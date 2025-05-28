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
Phy_merge2 <- Phy_merge %>% select(Phylum, everything()) #move last column to make it first column for ease of use
Phy_merge2 <- subset(Phy_merge2, select = -c(rn)) #removing column we don't need

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
Phy_merge2_c <- Phy_merge_c %>% select(Class, everything()) #move last column to make it first column for ease of use
Phy_merge2_c <- subset(Phy_merge2_c, select = -c(rn)) #removing column we don't need

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
Phy_merge2_o <- Phy_merge_o %>% select(Order, everything()) #move last column to make it first column for ease of use
Phy_merge2_o <- subset(Phy_merge2_o, select = -c(rn)) #removing column we don't need

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
Phy_merge2_f <- Phy_merge_f %>% select(Family, everything()) #move last column to make it first column for ease of use
Phy_merge2_f <- subset(Phy_merge2_f, select = -c(rn)) #removing column we don't need

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
Phy_merge2_g <- Phy_merge_g %>% select(Genus, everything()) #move last column to make it first column for ease of use
Phy_merge2_g <- subset(Phy_merge2_g, select = -c(rn)) #removing column we don't need

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
Phy_merge2_s <- Phy_merge_s %>% select(Species, everything()) #move last column to make it first column for ease of use
Phy_merge2_s <- subset(Phy_merge2_s, select = -c(rn)) #removing column we don't need

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

#write.csv(final_final, "./SERDP/Invasive_Spider_Fecal_Run/data/final_final_13Dec24.csv", row.names = TRUE) ## Save as file



################################### ANALYSIS ###########################################


##### SUBSETTING FOR ENDOSYMBIONTS FOR SIGNIFICANCE TESTING ######################################

#Remove columns of known endosymbionts based on literature and assume gut microbiota for the rest 

final_final <- read.csv("./SERDP/Invasive_Spider_Fecal_Run/data/final_final_13Dec24.csv")

subset_endo <- final_final %>%
  select("sampleid", "g__Wolbachia","g__Candidatus.Cardinium","g__Rickettsia","g__Rickettsiella") #Subset based on endosymbionts

names(subset_endo)[names(subset_endo) == 'sampleid'] <- 'Sample_ID'

metadata_merge <- merge(subset_endo, metadata_16s, by="Sample_ID")

no_blank_merged <- metadata_merge[-c(1,2,3), ] #Remove blanks from dataframe

subset_microbiota <- final_final %>%
  select (-"g__Wolbachia",-"g__Candidatus.Cardinium",-"g__Rickettsia",-"g__Rickettsiella") #Remove known endosymbionts to obtain file with remaining gut microbiota

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

# Add a very small pseudo-count to zero values to avoid issues with log or arcsin sqrt transformations
pseudo_count <- 1e-6
no_blank_merged$g__Wolbachia_no_zero <- no_blank_merged$g__Wolbachia + pseudo_count
no_blank_merged$g__Candidatus.Cardinium_no_zero <- no_blank_merged$g__Candidatus.Cardinium + pseudo_count
no_blank_merged$g__Rickettsia_no_zero <- no_blank_merged$g__Rickettsia + pseudo_count
no_blank_merged$g__Rickettsiella_no_zero <- no_blank_merged$g__Rickettsiella + pseudo_count

# Transform the abundance data to meet assumptions of normality for linear modeling:
# - Wolbachia abundance is transformed using arcsine square root (common for proportional data)
# - Cardinium, Rickettsia, and Rickettsiella abundances are log10 transformed
no_blank_merged$g__Wolbachia_arcsin <- asin(sqrt(no_blank_merged$g__Wolbachia_no_zero))
no_blank_merged$g__Candidatus.Cardinium_log10 <- log10(no_blank_merged$g__Candidatus.Cardinium_no_zero)
no_blank_merged$g__Rickettsia_log10 <- log10(no_blank_merged$g__Rickettsia_no_zero)
no_blank_merged$g__Rickettsiella_log10 <- log10(no_blank_merged$g__Rickettsiella_no_zero)

# Scale Cardinium and Rickettsia log10-transformed values (mean=0, sd=1) for modeling
no_blank_merged$g__Candidatus.Cardinium_scaled <- scale(no_blank_merged$g__Candidatus.Cardinium_log10)
no_blank_merged$g__Rickettsia_scaled <- scale(no_blank_merged$g__Rickettsia_log10)

# Fit a generalized linear model (GLM) for Wolbachia abundance as a function of Island, Species, and Specimen type
glm_wol <- glm(g__Wolbachia_arcsin ~ Island + Species.Identification + Type.of.Specimen, 
               data = no_blank_merged, family = gaussian)
summary(glm_wol)

# Check model residuals with Q-Q plot to assess normality assumption
qqnorm(residuals(glm_wol))
qqline(residuals(glm_wol))


# Cardinium is exclusively found in Tetragnatha on Oahu
# Convert Cardinium abundance column to numeric (sometimes it can be factor/character)
no_blank_merged$'g__Candidatus.Cardinium' <- as.numeric(no_blank_merged$'g__Candidatus.Cardinium')

# GLM for scaled Cardinium abundance with Island, Specimen type, and Species as predictors
glm_cand <- glm(g__Candidatus.Cardinium_scaled ~ Island + Type.of.Specimen + Species.Identification, 
                data = no_blank_merged, family = gaussian)
summary(glm_cand)

# Q-Q plot for Cardinium model residuals
qqnorm(residuals(glm_cand))
qqline(residuals(glm_cand))


# Rickettsia is exclusively found in Tetragnatha on the Big Island (Hawaii)
# Visualize distribution of Rickettsia (log10 transformed)
ggplot(no_blank_merged, aes(x = log10(g__Rickettsia_no_zero))) +
  geom_histogram()

# GLM for scaled Rickettsia abundance with Species, Specimen type, and Island as predictors
glm_rik_type <- glm(g__Rickettsia_scaled ~ Species.Identification + Type.of.Specimen + Island, 
                    data = no_blank_merged, family = gaussian)
summary(glm_rik_type)

# Q-Q plot for Rickettsia model residuals
qqnorm(residuals(glm_rik_type))
qqline(residuals(glm_rik_type))


# Rickettsiella is present in Steatoda and Scytodes on Maui and Oahu
# Note: Make sure arcsin sqrt transformed column for Rickettsiella exists (create if needed)
no_blank_merged$g__Rickettsiella_arcsin <- asin(sqrt(no_blank_merged$g__Rickettsiella_no_zero))

# GLM for Rickettsiella abundance using Species and Specimen type
glm_rikella <- glm(g__Rickettsiella_arcsin ~ Species.Identification + Type.of.Specimen, 
                   data = no_blank_merged, family = gaussian)
summary(glm_rikella)

# Q-Q plot for residuals of Rickettsiella model
qqnorm(residuals(glm_rikella))
qqline(residuals(glm_rikella))

# GLM for Rickettsiella abundance using Island and Specimen type
glm_rikella_site <- glm(g__Rickettsiella_arcsin ~ Island + Type.of.Specimen, 
                        data = no_blank_merged, family = gaussian)
summary(glm_rikella_site)

# Q-Q plot for residuals of site-based Rickettsiella model
qqnorm(residuals(glm_rikella_site))
qqline(residuals(glm_rikella_site))

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
  "Berkeley.Badumna" = "#005f73",
  "Hawaii (Big Island).Badumna" = "#94d2bd",
  "Hawaii (Big Island).Steatoda" = "#e9d8a6",
  "Maui.Steatoda" = "#ffba08",
  "Oahu.Steatoda" = "#f48c06",
  "Berkeley.Tetragnatha" = "#e01e37",
  "Hawaii (Big Island).Tetragnatha" = "#c71f37",
  "Maui.Tetragnatha" = "#a71e34",
  "Oahu.Tetragnatha" = "#85182a"
)

# Verify all Site_Species have assigned colors
missing_colors <- setdiff(unique(alpha_div_clean$Site_Species), names(color_mapping))
if (length(missing_colors) > 0) {
  stop(paste("Missing colors for:", paste(missing_colors, collapse = ", ")))
}

# =====================
### ANOVA TESTS ###
# =====================

anova_island <- aov(Shannon ~ Island, data = alpha_div_clean)
anova_type <- aov(Shannon ~ Type.of.Specimen, data = alpha_div_clean)
anova_species <- aov(Shannon ~ Species.Identification, data = alpha_div_clean)

# Check residuals visually for normality (repeat for each if desired)
qqnorm(residuals(anova_island)); qqline(residuals(anova_island), col = "red")

# ANOVAs within each island (species effects)
anova_by_island <- lapply(split(alpha_div_clean, alpha_div_clean$Island), function(df) {
  aov(Shannon ~ Species.Identification, data = df)
})

anova_summaries <- lapply(anova_by_island, summary)

print(anova_summaries$Berkeley)
print(anova_summaries$`Hawaii (Big Island)`)
print(anova_summaries$Maui)
print(anova_summaries$Oahu)

# Two-way ANOVA for Island x Species interaction
anova_interaction <- aov(Shannon ~ interaction(Island, Species.Identification), data = alpha_div_clean)
qqnorm(residuals(anova_interaction)); qqline(residuals(anova_interaction), col = "red")

# Tukey HSD posthoc test
tukey_results <- TukeyHSD(anova_interaction)
significant_pairs <- data.frame(tukey_results$`interaction(Island, Species.Identification)`) %>%
  dplyr::filter(p.adj < 0.05)

# Prepare comparisons for plotting
comparisons <- strsplit(rownames(significant_pairs), "-")

# =====================
### PLOT ###
# =====================

ggplot(alpha_div_clean, aes(x = interaction(Island, Species.Identification), y = Shannon, fill = Site_Species)) +
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
    x = "Site x Species Interaction",
    y = "Shannon-Weiner Diversity Index"
  ) +
  geom_signif(
    comparisons = comparisons,
    map_signif_level = TRUE,
    textsize = 3.5,
    tip_length = 0.02
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

# Plot gut microbiota top 10 genera
ggplot(gut_df_top, aes(x = as.factor(Species.Identification), y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_colors) +
  theme_classic() +
  labs(x = "Spider Genus", y = "Percentage of Abundance") +
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

# Categorize genera outside top 10 as "Other"
gut_df_other <- gut_df %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other"))

# Calculate relative abundance including "Other"
gut_df_other <- gut_df_other %>%
  group_by(Island, Species.Identification) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# Define palette with "Other" as white
palette_colors_other <- c(
  setNames(
    colorRampPalette(brewer.pal(9, "YlGnBu"))(length(top_genera)),
    top_genera
  ),
  Other = "white"
)

# Plot gut microbiota with "Other"
ggplot(gut_df_other, aes(x = as.factor(Species.Identification), y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_colors_other) +
  theme_classic() +
  labs(x = "Spider Genus", y = "Percentage of Abundance") +
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
# ENDOSYMBIONT ABUNDANCE PLOT
# =====================

ggplot(endosymbiont_df, aes(x = Species.Identification, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "OrRd") +
  theme_classic() +
  labs(x = "Species within Sites", y = "Abundance") +
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

### PERMANOVA for gut microbiota only
beta_dist_gut <- phyloseq::distance(ps_gut_filtered, method = "bray")
metadata_gut <- as.data.frame(as.matrix(sample_data(ps_gut_filtered)))

# Run PERMANOVA on gut microbiota
permanova_gut <- vegan::adonis2(
  beta_dist_gut ~ Island + Species.Identification + Type.of.Specimen,
  data = metadata_gut,
  by = "terms"
)

# Pairwise PERMANOVA for Species
pairwise_results <- pairwiseAdonis::pairwise.adonis2(
  beta_dist_gut ~ Species.Identification,
  data = metadata_gut,
  permutations = 999,
  by = "terms"
)

# Output results
print(permanova_gut)
print(pairwise_results)

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