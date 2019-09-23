# Networks
Networks for RNASeq data

The Rscript "coexpression_Spearman.R" can be used to create a coexpression network based on gene expression. 

Methodology: 

1. Calculates Spearman correlation for each pair of genes
2. Save the results in a large matrix 
3. You need to choose a threshold in order to get the most co-expressed genes (i.e: > 80 of Spearman)
4. Save the results in a Rdata file, then you don't need to run everything again 
5. Generate heatmaps
6. Create the final file with the gene source, gene target and weights. 

Now, you can open it on Cytoscape and have fun! 



