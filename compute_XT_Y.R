# say there are n indiv, m snp, k probes. for real data, m >> k

# G, data frame of 0/1/2 genotypes, n x m
# X, matrix of normalized genotypes, n x m
# Y, matrix of normalized phenotypes, n x k

# first 1k snps are all from chr 1
# chr 1 has 2k genes, chr 2 has 3k so hopefully result within first 5k 


# TODO why not use stdev(vals)/(N) for std err?

# find estimated std error (sigmaHat) for each snp, pheno pair
# return m x k matrix of sigmaHat's

findSigmaHats <- function(X,Y, N){
  if (false){ # for one h,i pair
  beta_hi = X_i^T*Y_h/N
  e_hi = Y_h - beta_hi*X_i # mean of Y_h is 0
  sigmaHat_hi = sqrt(t(e_hi) %*% e_hi)*sqrt(1/(N-2))
  }
  
  N = 119
  betas = (t(X) %*% Y)/N  # m x k
  
  # TODO what dim should es have
  es = Y - betas %*% X
  sigmaHats = sqrt()
}


# find t s.t. |S| > phi-1(alpha/2) iff | X^T*Y[snp,pheno] | >  t
# return m x k matrix of t's
findThesholds <- function(alpha, N, sigmaHats){
  # different sigmaHat for each snp-phenotype pair
  Sthresh <- qnorm(alpha/2)
  thresh <- sqrt(N) * Sthresh * sigmaHat
  return thresh
}

qnorm(0.025, lower.tail=FALSE)

snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 

Gt<-read.table(snps.txt.file, header=TRUE, sep="", row.names=1)
Yt.df<-read.table(expr.txt.file, header=TRUE, sep="", row.names=1)

G<-t(Gt)
X<-scale(G)  # standardize genotypes to mean 0, var 1
Y<-t(data.matrix(Yt.df))

# X^T = m x n, Y = n x k

Xt_Y <- t(X) %*% Y  # same as crossprod(X,Y)
temp <- which(Xt_Y > 40, arr.ind=TRUE)
temp_colnames <- colnames(Xt_Y)[temp[,2]]
temp_full <- cbind(temp, temp_colnames)


if (false){
TODO matrixify/normalize G
TODO convert Y.df to matrix, keep probe labels
TODO find probes and snps where X^T*Y > thresh
check whether probes are in correct region
}

if (false) {
# A = U S V^T
# U, columns are eigenvectors of A*A^T
# S, Sigma, diagonal matrix of eigenvalues
# V, columns are eigenvectors of A^T*A
}