## Raw Data for PhyloSeq Inputs

Here is the raw data file outputs from QIIME2 and inputs to the phyloseq analysis code. 

*filtered_clustered_ASV.csv* -  Amplicon Sequence Variant (ASV) count table  
- **Rows:** ASV IDs (Feature IDs)  
- **Columns:** Sample IDs  
- **Values:** Abundance of each ASV per sample

*merged_16S_tax_clusters.csv* -  Taxonomic annotations corresponding to clustered ASVs  
- **Rows:** ASV IDs (Feature IDs)  
- **Columns:** Taxonomic assignments from Kingdom to Species  
- Includes a “Confidence” column indicating match quality from SILVA

*spider_metadata_16s_only.csv* -  Sample metadata file  
- **Columns include:**  
  - `SampleID`  
  - `Species`  
  - `Site`  
  - `Region` (island/mainland)  
  - `Invasion_Status` (native/invasive)  
  - Additional sequencing prep info

*tree-cluster-A.nwk*, *tree-cluster-B.nwk*, *tree-cluster-C.nwk*
Phylogenetic tree files (Newick format) for three different clustering runs (A, B, and C), used for diversity and evolutionary analyses.
