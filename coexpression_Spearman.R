### Katia Lopes
### Ricardo Vialle 
### July 19, 2019
### Create a co-expression network from RNASeq data 

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

date()

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if(!require("limma")) BiocManager::install("limma"); library("limma")
if(!require("edgeR")) BiocManager::install("edgeR"); library("edgeR")
if(!require("amap")) BiocManager::install("amap"); library("amap")
if(!require("flashClust")) BiocManager::install("flashClust"); library("flashClust")
if(!require("plyr")) BiocManager::install("plyr"); library("plyr")
if(!require("sva")) BiocManager::install("sva"); library("sva")
if(!require("pamr")) BiocManager::install("pamr"); library("pamr")
if(!require("DESeq2")) BiocManager::install("DESeq2"); library("DESeq2")

work_plots = "/sc/orga/projects/pd-omics/katia/chicken/"

############################# This part you just need to run once! 
############################# INPUT expression table starts here. 

gene_tpm = read.table(paste0(work_plots, "E-MTAB-2797_chicken.tsv"), header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
gene_tpm = gene_tpm[,-1] # Exclude the 1st column of gene_ids because we'll use the ensembl ids 
dim(gene_tpm)
gene_tpm[is.na(gene_tpm)] <- 0 # change NAs for 0. Check if you can do it with your dataset! 
# rownames(gene_tpm) <- gene_tpm[,1] # To set a different rowname 

#In case you need to remove genes not expressed. All table was filtered before. 
#gene_tpm_nozero = gene_tpm[rowSums(gene_tpm)!=0,]

head(gene_tpm)
datExpr = t(gene_tpm) #datExpr is the new table of expression 

############################# INPUT expression table is done.
# function to crossvalidate correlation among all genes and calculate mean.
# Iterative mean to not overload RAM.
cv.coexp_mean <- function(datExpr,n.sampling,iter,cor.thresh)
{
  # data = rows are samples, columns are genes
  p <- matrix(0, nrow = dim(datExpr)[2], ncol = dim(datExpr)[2])
  rownames(p) <- colnames(datExpr)
  colnames(p) <- colnames(datExpr)
  mean_mx <- p
  thresh_mx <- p
  if (n.sampling > dim(datExpr)[1]){
    n.sampling = round(dim(datExpr)[1]*.9)
    print("n.sampling is higher than the number of samples. Using 90% instead.")
  }
  
  for(j in 1:iter)
  {
    random_samples <- sample(seq(1:dim(datExpr)[1]),size = n.sampling)
    gene_gene_cor = cor(datExpr[random_samples,],method = "spearman")
    mean_mx <- (mean_mx*(j-1) + gene_gene_cor)/j
    thresh_mx <- thresh_mx + (abs(gene_gene_cor) >= cor.thresh)*1
    cat("\r",j,"iterations.")
  }
  
  high_cor = as.data.frame(which(thresh_mx==iter, arr.ind=TRUE)) # High threshold in all interations 
  rownames(high_cor) = c()
  high_cor = unique(high_cor[!high_cor$row==high_cor$col,]) # Remove diagonal
  high_cor_mean_mx = mean_mx[unique(high_cor$row),unique(high_cor$row)]
  
  ret_obj = list()
  ret_obj$mean_mx = mean_mx #mean of the correlation after iterations 
  ret_obj$thresh_mx = thresh_mx #number of times higher than the threshold
  ret_obj$high_cor_mean_mx = high_cor_mean_mx #correlation table only with pairs higher than the threshold in all iterations 
  
  return(ret_obj)
}

cor_mean = cv.coexp_mean(datExpr, 5, 10, 0.85) #number of samples to keep, iterations, threshold. USE at least 100 iterations! 
dim(cor_mean$mean_mx)

save(cor_mean,  file = paste0(work_plots, "cor_mean.Rdata")) 

############################ After save the data you just need to load for the 
# downstream analysis! 

work_plots = "/Users/katia/OneDrive/Documentos/Katia_github/Networks/"

load(paste0(work_plots, "cor_mean.Rdata"))

library(pheatmap)

#Plots 
png(paste0(work_plots, "HM_mean_mx.png"), width = 30, height = 30, res = 300, units = "in")
pheatmap(cor_mean$mean_mx, legend = F, labels_col = F, labels_row = F, show_rownames = F, show_colnames = F)
dev.off()

dim(cor_mean$high_cor_mean_mx)

png(paste0(work_plots, "HM_high_mean_mx.png"), width = 10, height = 10, res = 300, units = "in")
pheatmap(cor_mean$high_cor_mean_mx, legend = F, labels_col = F, labels_row = F, show_rownames = F, show_colnames = F)
dev.off()

high_cor = as.data.frame(which(abs(cor_mean$mean_mx)>0.85, arr.ind=TRUE)) # High threshold in all interations 
high_cor = high_cor[!high_cor$row==high_cor$col,] # Remove diagonal
genes4net = cor_mean$mean_mx[unique(high_cor$row),unique(high_cor$row)] #Select only genes with high correlation to show in the network 
dim(genes4net)

png(paste0(work_plots, "HM_genes4net.png"), width = 30, height = 30, res = 300, units = "in")
pheatmap(genes4net, legend = F, labels_col = F, labels_row = F, show_rownames = F, show_colnames = F)
dev.off()

library(reshape2) 
genes4net[upper.tri(genes4net)] <- NA
file4cytoscape = setNames(melt(genes4net, na.rm = TRUE), c('source', 'target', 'Spearman_corr'))
file4cytoscape = file4cytoscape[which( !file4cytoscape$source==file4cytoscape$target ),]
file4cytoscape_high = file4cytoscape[which(abs(file4cytoscape$Spearman_corr)>0.85),]

write.table(file4cytoscape_high, file = paste0(work_plots, "coexpression_net4cytos.txt"), sep = "\t", quote = F, row.names = F)



