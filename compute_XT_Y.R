# compute Beta using matrix multiplication
# code for input and shared functions is in input.R

library(reshape2)

snps.txt.file<-"GTEx_data/Lung10k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 

X<-importX(snps.txt.file)
Y<-importY(expr.txt.file)

n<-nrow(X); m<-ncol(X); k<-ncol(Y)
# X^T = m x n, Y = n x k

Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
if (FALSE){write.table(Beta, file = "Beta.Rmatrix.tsv", sep="\t")}

# find sigma, pvals, thresholds
SigmaHat<-findSigmaHats(X,Y)
Pvals<-findPvals(Beta, SigmaHat, n)
Thresh_Beta<-findThresholds(alpha=10^-8, n=nrow(X), SigmaHat=SigmaHat)
SigPairs <- findSigPairs(Beta, Thresh_Beta) # list of sig snp-pheno pairs

# make long dataframe of snp, pheno, beta, pval
Beta_long<-melt(Beta, id.vars=1)
Pvals_long<-melt(Pvals, id.vars=1)
Results_xty <- cbind(Beta_long, Pvals_long[,3])
colnames(Results_xty) <- c("snpID", "gene", "beta_true", "pval_true")
rm("Beta_long, Pvals_long") # rm once they are in Results_xty

# sort snp-pheno pairs by pval
# commented out bc sort the combined results instead
if (FALSE){
Results_xty_sorted<-Results_xty[order(Results_xty$pval),]
rownames(Results_xty_sorted)<-NULL
}

## examine first gene ##
if (FALSE) {
  Beta_123 <- as.matrix(Beta[,1])
  rownames(Beta_123)<-rownames(Beta); colnames(Beta_123)<-colnames(Beta)[1] 
  
  S_123<-Beta_123[,1] * sqrt(nrow(X)) / SigmaHat[,1]
  pvals_123 <- findPvals(Beta_123[,1], SigmaHat[,1], nrow(X))
  Beta_123 <- cbind(Beta_123, pvals_123)
  
  snpNo <- c(1:nrow(Beta_123));
  Beta_123 <- cbind(snpNo, Beta_123) # add row index pre-sorting
  #Beta_123 <- Beta_123[order(Beta_123[,2]),] # sort by pval
  
  if (FALSE){write.table(Beta_123, file = "Beta_first_gene.Rmatrix.tsv", sep="\t")}
}