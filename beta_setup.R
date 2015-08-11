# this file has input
# and other functions that are shared between methods

# say there are n indiv, m snp, k probes. for real data, m >> k

# G, data frame of 0/1/2 genotypes, n x m
# X, matrix of normalized genotypes, n x m
# Y, matrix of normalized phenotypes, n x k

setwd('~/Github/multi_pheno')

### FUNCTIONS ###

importData<-function(){
  # first 1k snps are all from chr 1
  # chr 1 has 2k genes, chr 2 has 3k so hopefully result within first 5k 
  snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
  expr.txt.file<-"GTEx_data/Lung30.expr.txt" 
  
  Gt<-read.table(snps.txt.file, header=TRUE, sep="", row.names=1)
  Yt.df<-read.table(expr.txt.file, header=TRUE, sep="", row.names=1)
  
  G<-t(Gt)
  X<-scale(G)  # standardize columns (genotypes) to mean 0, var 1
  Y<-t(data.matrix(Yt.df))
  return (X,Y)
}

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
      e_hi <- Y[,h] - Beta[i,h] * X[,i] # 1 x n vector of residuals
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



# get list of snp-phenotype which cross threshold
findSigPairs <- function(Beta, Thresh_Beta){
  sigPairs <- which(abs(Beta) > Thresh_Beta, arr.ind=TRUE) # indices of sig snp-pheno pairs
  
  # add column names, beta threshold, beta value 
  sigPairs <- cbind(sigPairs,
              colnames(Beta)[sigPairs[,2]],
              Thresh_Beta[sigPairs],
              Beta[sigPairs]
              )
  colnames(sigPairs_results)<-c("row", "col", "pheno", "threshold_ih","beta_ih")
  return sigPairs_results
}


# find pvals (using Beta and SigmaHat to find ncp)
findPvals <- function(Beta, SigmaHat, n){
  # Note assoc stat S =  Beta * sqrt(n)/SigmaHat
  # Note pval = 2 * phi(|Beta|*sqrt(n)/SigmaHat)
  Pvals <- 2*pnorm(abs(Beta) * sqrt(n) / SigmaHat, lower.tail=FALSE)
  return (Pvals) # return matrix of pvals
}

