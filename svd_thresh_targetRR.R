# find 99% thresh 
# i.e. thresh for pval_svd st. pval_svd < thresh contains any 99% of true
# (so if 100 true vals, want lowest thresh that gets any 99, not top 99)
library(reshape2)

RR_target = 0.99    # target recall rate
thresh_true = 10^-6 # desired threshold

# finds minimum pval_svd needed to achieve recall rate, for given rank
findThreshSvd_oneRank <- function(topTrue, RR_target, thresh_true){
  numTrue <- nrow(topTrue)
  numRecalled_target <- numTrue * RR_target # min number to recall
  
  thresh_svd <- topTrue$pval_svd[1] # worst case, take max pval_svd in top values
  for (i in 1:numTrue){
    # reduce thresh_svd to this lower thresh - do enough values pass new thresh_svd?
    new_thresh_svd <- topTrue$pval_svd[i]
    newRecalled <- topTrue[which(topTrue$pval_svd <= new_thresh_svd),]
    newNumRecalled <- length(which(newRecalled$pval_true <= thresh_true))
    #print(c(i, new_thresh_svd, newNumRecalled))
    if (newNumRecalled < numRecalled_target){ # if new svd thresh is too strict
      break # break without changing thresh_svd from last iteration
    }
    
    thresh_svd <- new_thresh_svd
  }
  
  return (thresh_svd)
}


# Find threshold and correlation for multiple ranks
# (Note: if only examining top true pvals, discard other values early to save time)
# (cbind of all snp-pheno take much longer than cbind for only a few snp-pheno)
findThreshSvd_multiRank <- function(RR_target, thresh_true){
  # find XT_Y values 
  snps.txt.file<-"GTEx_data/Lung10k.snps.txt" 
  expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 
  X<-importX(snps.txt.file)
  Y<-importY(expr.txt.file)
  n<-nrow(X); m<-ncol(X); k<-ncol(Y) # X^T = m x n, Y = n x k
  
  Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
  
  # find sigma, pvals, thresholds
  SigmaHat<-findSigmaHats(X,Y)
  Pvals<-findPvals(Beta, SigmaHat, n)
  Thresh_Beta<-findThresholds(alpha=10^-8, n=nrow(X), SigmaHat=SigmaHat)
  
  # make long dataframe of snp, pheno, beta, pval
  Beta_long<-melt(Beta, id.vars=1)
  Pvals_long<-melt(Pvals, id.vars=1)
  Results_xty <- cbind(Beta_long, Pvals_long[,3])
  colnames(Results_xty) <- c("snpID", "gene", "beta_true", "pval_true")
 
  topTrueInds <- which(Results_xty$pval_true <= thresh_true) 
  #Results_xty <- Results_xty[topTrueInds,] # use if only want top snp-pheno
  
  # finds SVD values for each rank and then finds thresh_svd
  doStuffRank <- function(rank){
    # run svd with given rank, find betas and pvals
    svdX <- svd(X, nu=rank, nv=rank) 
    
    if (rank > 1) {
      Dx <- diag(svdX$d[1:rank])
    } else {
      Dx <- svdX$d[1] # if rank 1, diagonal is just one entry/scalar
    }
    
    X2 <- svdX$u %*% Dx %*% t(svdX$v)  # X = U D V'
    rownames(X2) <- rownames(X); colnames(X2) <- colnames(X)
    Beta_svd <- crossprod(X2, Y)/n
    SigmaHat_svd<-findSigmaHats(X2,Y)
    Pvals_svd<-findPvals(Beta_svd, SigmaHat_svd, n)
    
    # make long dataframe of snp, pheno, beta, pval
    Beta_long<-melt(Beta_svd, id.vars=1)
    Pvals_long<-melt(Pvals_svd, id.vars=1)
    Results_svd <- cbind(Beta_long, Pvals_long[,3])
    colnames(Results_svd) <- c("snpID", "gene", "beta_svd", "pval_svd")
    # use if don't need snp-pheno that don't pass thresh
    #Results_svd <- Results_svd[topTrueInds,] 
    
    Results_all<-cbind(Results_svd, 
                       Results_xty[,3:4],
                       rank(Results_svd$pval_svd)/(m*k),
                       rank(Results_xty$pval_true)/(m*k)
    )
    colnames(Results_all)[7] <- "percentile_pval_svd"
    colnames(Results_all)[8] <- "percentile_pval_true"
    
    topTrue <- Results_all[topTrueInds,]
    
    threshSvd <- findThreshSvd_oneRank(topTrue, RR_target, thresh_true)
    corr_top <- cor(topTrue$pval_svd, topTrue$pval_true)
    corr_all <- cor(Results_all$pval_svd, Results_all$pval_true)
    
    return (rbind(threshSvd, corr_top, corr_all))
  }

  ThreshSvd_Corr <- vapply(1:119, doStuffRank, FUN.VALUE=double(3))
  ThreshSvd_Corr <- t(ThreshSvd_Corr)
  colnames(ThreshSvd_Corr) <- c("threshSvd", "corr_top", "corr_all")
  return (ThreshSvd_Corr)
}

ThreshSvd_Corr <- findThreshSvd_multiRank(RR_target = RR_target,
                                     thresh_true = thresh_true)

ThreshSvd_Corr <- cbind(1:119, ThreshSvd_Corr)
colnames(ThreshSvd_Corr)[1] <- "svdX_rank"

rownames(ThreshSvd_Corr)<-NULL
write.table(ThreshSvd_Corr, "thresh_corr_by_svdXrank_10ksnp_5kpheno.Rtable.txt")

ThreshSvd_Corr<-as.data.frame(ThreshSvd_Corr)
# TODO plot Results_1k_30
# TODO plot Results_10k_5k
# can also plot correlation for each rank, for all points

dev.off()

par(mfrow=c(1,3))
attach(ThreshSvd_Corr)
plot(svdX_rank, threshSvd, ylim=c(0,1), cex.lab=1.5, main="svd threshold for RR = 99%, alpha = 10^-6")
plot(svdX_rank, corr_top, cex.lab=1.5, main="correlation of pval_svd, pval_true, for snp-pheno's with pval_true >= 10^-6")
plot(svdX_rank, corr_all, cex.lab=1.5, main="correlation of pval_svd, pval_true, for all snp-pheno's")


dev.off()

par(mfrow=c(1,3))
attach(ThreshSvd_Corr)
plot(svdX_rank, threshSvd, cex.lab=1.5, main="svd threshold for RR = 99%, alpha = 10^-6", log="y")
plot(svdX_rank, corr_top, cex.lab=1.5, main="correlation of pval_svd, pval_true, for snp-pheno's with pval_true >= 10^-6", log="y")
plot(svdX_rank, corr_all, cex.lab=1.5, main="correlation of pval_svd, pval_true, for all snp-pheno's", log="y")

#point(ThreshSvd, Corr_top, Corr_all)
#Results_1k_30 <- ThreshSvd_Corr
