# say there are n indiv, m snp, k probes. for real data, m >> k

# THIS IS A COPY OF COMPUTE_XT_Y WITH A FEW MODS FOR SVD
# merge the two files or separate out the shared parts later

# G, data frame of 0/1/2 genotypes, n x m
# X, matrix of normalized genotypes, n x m
# Y, matrix of normalized phenotypes, n x k

setwd('~/Github/multi_pheno')

### FUNCTIONS ###
# find estimated std error (sigmaHat) for each snp, pheno pair
# return m x k matrix of sigmaHat's
# TODO Beta now passed in, not computed from X,Y in function body
findSigmaHats <- function(X,Y, Beta){
  n<-nrow(X)
  m<-ncol(X)
  k<-ncol(Y)
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
### instead of cross prod, to get Beta, use SVD ###
#Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
svdX <- svd(X, nu=110, nv=70)
#dim(svdX$u); dim(diag(svdX$d)); dim( t(svdX$v))
#Dx <- diag(svdX$d)
Dx <- t(svdX$u) %*% X %*% svdX$v # D = U' X V
X2 <- svdX$u %*% Dx %*% t(svdX$v) #  X = U D V'
rownames(X2)<-rownames(X2); colnames(X2)<-colnames(X)
Beta <- crossprod(X2, Y)
Beta_true <- crossprod(X,Y)


SigmaHat<-findSigmaHats(X,Y,Beta)
pvals <- 2*pnorm(abs(Beta) * sqrt(nrow(X)) / SigmaHat[,1], lower.tail=FALSE)
pvals_true <- 2*pnorm(abs(Beta_true) * sqrt(nrow(X)) / SigmaHat[,1], lower.tail=FALSE)

svdY <- svd(Y, nu=7,nv=9)
dim(svdY$u); dim(diag(svdY$d)); dim(svdY$v)
Dy <- diag(svdY$d)
Y2 <- svdY$u %*% Dy %*% t(svdY$v)

u1 <- svdX$u
v1 <- Dx %*% t(svdX$v)
u2 <- svdY$u
v2 <- Dy %*% t(svdY$v)

dim(u1); dim(v1); dim(u2); dim(v2)

Beta_true = crossprod(X,Y)
Beta_1 = t(v1) %*% t(u1) %*% u2 %*% v2
Beta_2 = crossprod(X2, Y2)
rownames(Beta_2) <- rownames(Beta_true); colnames(Beta_2) <- colnames(Beta_true)

if (FALSE){
  write.table(Beta, file = "Beta.Rmatrix.tsv", sep="\t")
}


# TODO no pairs cross threshold for alpha=10^-8 in first 30 pheno
#     and out of 5k phenos, only 50 cross threshold for alpha=10^-8.
#     does this match expected?

SigmaHat<-findSigmaHats(X,Y,Beta)
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