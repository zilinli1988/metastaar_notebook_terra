# metastaar_notebook_terra
Notebook Script for MetaSTAAR in Terra.
## Description
These are the notebook scripts for rare variant meta-analysis of WGS data based on the Meta-analysis of 
variant-Set Test for Association using Annotation infoRmation (MetaSTAAR) method, 
including (1) generating summary statistics of rare variants (RVs) and (2) performing rare variant meta-analysis. 
## Usage
Container image: zilinli/metastaar_dockerimage_terra_notebook.

MetaSTAARWorker_Cov_Generation.R and MetaSTAAR_Wroker_Score_Generation.R are the notebook scripts for generating summary statistics of rare variants based on MetaSTAARWorker.  
Specifically, MetaSTAARWorker_Cov_Generation.R is used to generate the sparse weighted linkage disequilibrium (LD) matrices. 
MetaSTAARWorker_Score_Generation.R is used to generate the individual variant summary statistics of all variants across the genome.

MetaSTAAR_RVAT_Sliding_Window.R and MetaSTAAR_Individual_Analysis.R are the notebook scripts for performing rare variant meta-analysis using the summary statistics generated by MetaSTAARWorker.
Specifically, MetaSTAAR_RVAT_Sliding_Window.R is used to perform sliding window analysis across the genome. 
MetaSTAAR_Individual_Analysis.R is used to perform individual analysis of common variants across the genome. 

For the access of TOPMed aGDS files in Terra, please send an email to li@hsph.harvard.edu for access.
TOPMed Freeze 5 aGDS files are in the workspace "MetaSTAAR_Lipids_F8".
TOPMed Freeze 8 aGDS files are in the workspace "topmed_agds". 
## Computation Cost
MetaSTAARWorker_Estimated_Cost.pdf is the benchmarked cost of MetaSTAARWorker in Terra.
