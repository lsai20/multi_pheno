# compute Beta using matrix multiplication
# code for input and shared functions is in input.R
### SCRIPT ###

X,Y<-importData()

# X^T = m x n, Y = n x k
Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
if (FALSE){write.table(Beta, file = "Beta.Rmatrix.tsv", sep="\t")}

SigmaHat<-findSigmaHats(X,Y)
Thresh_Beta<-findThresholds(alpha=10^-8, n=nrow(X), SigmaHat=SigmaHat)
sigPairsResults <- findSigPairs(Beta, Thresh_Beta) # list of sig snp=pheno pairs

Beta_123 <- as.matrix(Beta[,1])
rownames(Beta_123)<-rownames(Beta)
colnames(Beta_123)<-colnames(Beta)[1] 
S<-Beta_123[,1] * sqrt(nrow(X)) / SigmaHat[,1]
#phi_inverse <- pnorm(Beta_123[,1] / sqrt(nrow(X)) * SigmaHat[,1]) 
# pval !=  phi^-1(     Beta_ih * sqrt(n) / Sigma_ih) 
# pval = 2*phi^-1(abs(Beta_ih) * sqrt(n) / Sigma_ih)
pvals <- 2*pnorm(abs(Beta_123[,1]) * sqrt(nrow(X)) / SigmaHat[,1], lower.tail=FALSE)
Beta_123 <- cbind(Beta_123, pvals)
snpNo <- c(1:nrow(Beta_123))
Beta_123 <- cbind(snpNo, Beta_123) # add row index pre-sorting
#Beta_123 <- Beta_123[order(Beta_123[,2]),] # sort by pval

if (FALSE){
  write.table(Beta_123, file = "Beta_first_gene.Rmatrix.tsv", sep="\t")
}



if (FALSE){
  #DONE matrixify/normalize G
  #DONE convert Y.df to matrix, keep probe labels
  #DONE find probes and snps where X^T*Y > thresh
  #check whether probes are in correct region
}

if (FALSE) {
  # A = U S V^T
  # U, columns are eigenvectors of A*A^T
  # S, Sigma, diagonal matrix of eigenvalues
  # V, columns are eigenvectors of A^T*A
}