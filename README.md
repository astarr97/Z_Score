# Z_Score
Code related to the Z Score project and for the paper "Accounting for cis-regulatory constraint prioritizes genes likely to affect species-specific traits".  

Prepare_For_MWU_Avg_Reads.ipynb contains the code used to preprocess the data and compare the human ASE distribution to the human/chimpanzee ASE distribution.
Pre_Post_Enrichment.ipynb contains code to process the resulting p-values, analyze enrichments, and make plots.
Used_For_Enrichment_Cur.ipynb contains the code for enrichment testing with GSEAPY
Single_Cell.ipynb contains code for processing and plotting the Kanton et al. single cell dataset.
Plot_Rachel.ipynb contains code to plot data from Agoglia et al. 2021.
Plotting_Expression.ipynb contains code to plot data from a variety of publications.
DEA_No_Function.R was used for testing for differential expression with DESeq2 and to compute TPM.
Compute_TPM.R also computes TPM.
Finally, the various Snakefiles were used to map data and count reads.  No dedup doesn't deduplicate reads, TrimGalore uses trimgalore to trim, and Hybrid is for hybrid data.  

There is no restriction on the use of this code whatsoever.
