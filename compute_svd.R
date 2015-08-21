# computes/approximates Beta using SVD
# see beta_setup.R for additional fxns/reading input
library(reshape2) # for melt


# TODO for full file, will have to window X and/or Y
snps.txt.file<-"GTEx_data/Lung10k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 

X<-importX(snps.txt.file)
Y<-importY(expr.txt.file)

# say there are n indiv, m snp, k probes. for real data, m >> k
n<-nrow(X); m<-ncol(X); k<-ncol(Y)

# test whether top percent_true snps are in top percent_svd snps
percent_svd <- 0.05 # use percentile ranges 0-1, not 0-100%
percent_true <- 0.01


# say we use b components to approximate X. exact is b = min(n, m).
# and c components to approximate Y. exact is c = min(n,k)
b<-100
# seems like a lot of improvement from b=60 to b=70 - plot/quantitfy?
c<-min(n,k)

### instead of cross prod, use SVD to get Beta ###
svdX <- svd(X, nu=b, nv=b) 
Dx <- diag(svdX$d[1:b])

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
rm("Beta_long", "Pvals_long") # rm once they are in Results_svd

# sort snp-pheno pairs by pval
# commented out bc sort the combined results instead
if (FALSE){
Results_svd_sorted<-Results_svd[order(Results_svd$pval),]
rownames(Results_svd_sorted)<-NULL
}

### Compare to results of X^T*Y method ###

# attach "true" (X^T*Y) beta and pvals to df
# get percentiles
Results_all<-cbind(Results_svd, 
                   Results_xty[,3:4],
                   rank(Results_svd$pval_svd)/(m*k),
                   rank(Results_xty$pval_true)/(m*k)
                  )
colnames(Results_all)[7] <- "percentile_pval_svd"
colnames(Results_all)[8] <- "percentile_pval_true"


# get rank instead of percentiles
if (FALSE){
Results_all<-cbind(Results_svd, 
                   Results_xty[,3:4],
                   rank(Results_svd$pval_svd),
                   rank(Results_xty$pval_true)
)
colnames(Results_all)[7] <- "rank_pval_svd"
colnames(Results_all)[8] <- "rank_pval_true"

topTrue <- Results_all[which(Results_all$rank_pval_true <= percent_true*m*k),]

}

if(FALSE){
write.table(Results_all, file = "Results_svdX70_8-17-15.Rmatrix")
write.table(topTrue, file = "topTrue_svdX70_8-17-15.Rmatrix")

}

# Pick top percent_true of snp-pheno pairs by true pval
# Are they ranked highly by svd pval as well?
topTrue <- Results_all[which(Results_all$percentile_pval_true <= percent_true),]
topTrue <- topTrue[order(topTrue$pval_true),]
topTrue <- topTrue[order(topTrue$pval_svd, decreasing=T),]


# (TODO later - sort and pick by svd threshold, see how many cross true thresh)
# (TODO get theoretical thresh)

# TODO use svd of Y or not?
if (FALSE){
svdY <- svd(Y, nu=c, nv=c)
Dy <- diag(svdY$d[1:c])
Y2 <- svdY$u %*% Dy %*% t(svdY$v)  # Y = U D V'

u1 <- svdX$u
v1 <- Dx %*% t(svdX$v)
u2 <- svdY$u
v2 <- Dy %*% t(svdY$v)
}