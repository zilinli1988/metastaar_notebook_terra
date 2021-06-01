# metastaarworker_notebook_terra
Notebook scripts for generating summary-level data of WGS for TOPMed data in Terra.
## Description
These are the R scripts for generating summary-level data of WGS based on the Meta-analysis of variant-Set Test for Association using Annotation infoRmation (MetaSTAAR) method, including the LD matrix of rare variants (RVs) and the summary statistics. 
## Usage
Container image: zilinli/metastaar_dockerimage_terra_notebook.

aGDS file: send an email to li@hsph.harvard.edu for access.

MetaSTAAR_Worker_Cov_Generation.R is used to generate the variance-covariance matrices to represent the linkage disequilibrium (LD) structure among RVs. 

MetaSTAAR_Wroker_Score_Generation.R is used to generate the summary statistics of all variants across the genome.

