# 2023-Golconda_A19_paper
Computational analyses published in Golconda and Oberstein, 2023 (A19 manuscript).

To use these scripts, make sure to unzip all gzipped files after downloading the git repo.

These are required for the analyses in folders 3, 4, and 5 (clustering, plotting, and MLSeq):

* ./tx2gene_vector_TB40-nonCanonical_TB40_gencode34.transcripts.gz

* ./ttg.tab.gz

* ./2-limmaGene_GeTMM/GeTMM_allData_TMM_expValues_cpm.tab.gz

These are required for normalization of the gene expression matrix in folder "2-limmaGene_GeTMM":

* ./GSE*/SRR*/abundance.tsv.gz

## Organization:
* Processing scripts are located in numbered folders, e.g. "2-limmaGene_GeTMM"
  
* GSEXXXXXXXX prefixed folders contain alignment logs and kallisto output files for each study. These are only required to generate the normalized gene expression matrix using the script in "2-limmaGene_GeTMM."
  
* All subsequent analyses use the normalized gene expression matrix found in folder "2-limmaGene_GeTMM" named "GeTMM_allData_TMM_expValues_cpm.tab."
