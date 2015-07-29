# say there are n indiv, m snp, k probes. for real data, m >> k

# G, data frame of 0/1/2 genotypes, n x m
# X, matrix of normalized genotypes, n x m
# Y, matrix of normalized phenotypes, n x k

setwd('~/Github/multi_pheno')

### FUNCTIONS ###
# find estimated std error (sigmaHat) for each snp, pheno pair
# return m x k matrix of sigmaHat's
findSigmaHats <- function(X,Y){
  n<-nrow(X)
  m<-ncol(X)
  k<-ncol(Y)
  Beta <- crossprod(X,Y)/n  # m x k matrix of betas
  SigmaHat <- matrix(data=NA, nrow=m, ncol=k)
  # nested for loop probably bad R
  # for one h,i pair
  for (i in 1:m){
    for (h in 1:k){
      e_hi <- Y[,h] - Beta[i,h]*X[,i] # 1 x n vector of residuals 
        # ^also note mean of Y_h is 0
      SigmaHat[i,h] <- sqrt( crossprod(e_hi) / (n-2) )
    }
  }
  rownames(SigmaHat)<-colnames(X)
  colnames(SigmaHat)<-colnames(Y)
  return (SigmaHat)
}


# find thresholds for Beta_ih (as opposed to S_ih)
# find t s.t. |S| > phi-1(alpha/2) iff | beta | = | X^T*Y[snp,pheno]/n_ih | >  t
# return m x k matrix of t's
# note: input SigmaHat is an m x k matrix of estimated std err
findThresholds <- function(alpha, n, SigmaHat){
  # different sigmaHat for each snp-phenotype pair
  Sthresh <- qnorm(alpha/2, lower.tail=FALSE) # the threshold if using standardized assoc stat S
  Thresh_Beta <- Sthresh * SigmaHat / sqrt(n)
  rownames(Thresh_Beta)<-rownames(SigmaHat)
  colnames(Thresh_Beta)<-colnames(SigmaHat)
  return (Thresh_Beta)
}


### SCRIPT ###
# first 1k snps are all from chr 1
# chr 1 has 2k genes, chr 2 has 3k so hopefully result within first 5k 

snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung30.expr.txt" 

Gt<-read.table(snps.txt.file, header=TRUE, sep="", row.names=1)
Yt.df<-read.table(expr.txt.file, header=TRUE, sep="", row.names=1)

G<-t(Gt)
X<-scale(G)  # standardize columns (genotypes) to mean 0, var 1
Y<-t(data.matrix(Yt.df))
# X^T = m x n, Y = n x k
Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
if (FALSE){
  write.table(Beta, file = "Beta.Rmatrix.tsv", sep="\t")
}


# TODO no pairs cross threshold for alpha=10^-8 in first 30 pheno
#     and out of 5k phenos, only 50 cross threshold for alpha=10^-8.
#     does this match expected?

SigmaHat<-findSigmaHats(X,Y)
Thresh_Beta<-findThresholds(alpha=10^-8, n=nrow(X), SigmaHat=SigmaHat)
sigPairs <- which(abs(Beta) > Thresh_Beta, arr.ind=TRUE) # indices of sig snp-pheno pairs

# add column names, beta threshold, beta value 
sigPairs_results <- cbind(sigPairs, 
                          colnames(Beta)[sigPairs[,2]],
                          Thresh_Beta[sigPairs],
                          Beta[sigPairs]
                          )
colnames(sigPairs_results)<-c("row", "col", "pheno", "threshold_ih","beta_ih")


Beta_123 <- as.matrix(Beta[,1])
rownames(Beta_123)<-rownames(Beta)
colnames(Beta_123)<-colnames(Beta)[1]
# TODO is this right? phi^-1(Beta_ih * sqrt(n) / Sigma_ih)
pvals <- pnorm(Beta_123[,1] * sqrt(nrow(X)) / SigmaHat[,1])
Beta_123 <- cbind(Beta_123, pvals)
Beta_123 <- cbind(Beta_123, c(1:nrow(Beta_123))) # add row index pre-sorting

#TODO don't think i'm getting orig snp no - has beta already been sorted somewhere?
origSnpNo <- order(Beta[,2]) # orig snp no
Beta_123 <- Beta_123[order(Beta_123[,2]),] # sort by pval
Beta_123 <- cbind(Beta_123, origSnpNo)

colnames(Beta_123)

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