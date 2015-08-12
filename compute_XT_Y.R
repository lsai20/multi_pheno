# compute Beta using matrix multiplication
# code for input and shared functions is in input.R
### SCRIPT ###
snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung30.expr.txt" 

X<-importX(snps.txt.file)
Y<-importY(expr.txt.file)

n<-nrow(X); m<-ncol(X); k<-ncol(Y)
# X^T = m x n, Y = n x k

Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
if (FALSE){write.table(Beta, file = "Beta.Rmatrix.tsv", sep="\t")}

SigmaHat<-findSigmaHats(X,Y)
Pvals<-findPvals(Beta, SigmaHat, n)
Thresh_Beta<-findThresholds(alpha=10^-8, n=nrow(X), SigmaHat=SigmaHat)
SigPairs <- findSigPairs(Beta, Thresh_Beta) # list of sig snp-pheno pairs


## examine first gene ##
Beta_123 <- as.matrix(Beta[,1])
rownames(Beta_123)<-rownames(Beta); colnames(Beta_123)<-colnames(Beta)[1] 

S_123<-Beta_123[,1] * sqrt(nrow(X)) / SigmaHat[,1]
pvals_123 <- findPvals(Beta_123[,1], SigmaHat[,1], nrow(X))
Beta_123 <- cbind(Beta_123, pvals_123)

snpNo <- c(1:nrow(Beta_123));
Beta_123 <- cbind(snpNo, Beta_123) # add row index pre-sorting
#Beta_123 <- Beta_123[order(Beta_123[,2]),] # sort by pval

if (FALSE){write.table(Beta_123, file = "Beta_first_gene.Rmatrix.tsv", sep="\t")}

