# README

This repository was made to store data and code for the publication entitled "Do invasive spiders escape their native microbiomes? Patterns of microbial variation in native and invasive spiders in Hawai‘i" 

## Usage
This repository contains the following key files and folders:
- `filtered_clustered_ASV.csv`  
  Amplicon Sequence Variant (ASV) count table with abundance values for each sample.
- `merged_16S_tax_clusters.csv`  
  Taxonomic annotations for each ASV, with SILVA v138 confidence scores.
- `spider_metadata_16S_only.csv`  
  Metadata including sample IDs, spider species, site, region, invasion status, and sequencing info.
- `tree-cluster-A.nwk`, `tree-cluster-B.nwk`, `tree-cluster-C.nwk`  
  Phylogenetic tree files corresponding to different clustering runs.
- `cluster_analysis_figures.R`  
  R script for preprocessing microbial data, building phyloseq objects, statistical analyses, and figure creation.
- `spid_picrust2_clusters.R`  
  R script to process PICRUSt2 functional predictions and merge data for downstream analysis.

---

## File Organization

Data is organized in the following fashion: 
- **RawData/** — Data inputs (.csv) to phyloseq
- **SharedData/** — links to Dryad DOI (contains raw data (.FASTQ) and metadata files) 
- **AnalysisCode/** — R scripts for phyloseq analysis and PICRUSt2 processing (public)

---

## Contributing

Pull requests and suggestions are welcome! For major changes or questions, please open an issue or contact Madison Pfau at [mapfau@ucsc.edu](mailto:mapfau@ucsc.edu).

---

## Contact

If you have questions about the data or analysis scripts, please reach out to:  
Madison J. Pfau — [mapfau@ucsc.edu](mailto:mapfau@ucsc.edu)

---

*This README was generated to promote transparency and reproducibility in microbial ecology research.*

