# computes/approximates Beta using SVD
# see beta_setup.R for additional fxns/reading input
library(reshape2) # for melt

snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 

X<-importX(snps.txt.file)
Y<-importY(expr.txt.file)

# say there are n indiv, m snp, k probes. for real data, m >> k
n<-nrow(X); m<-ncol(X); k<-ncol(Y)

# say we use b components to approximate X. exact is b = min(n, m).
# and c components to approximate Y. exact is c = min(n,k)
b<-100
c<-min(n,k)

### instead of cross prod, use SVD to get Beta ###
svdX <- svd(X, nu=b, nv=b) 
Dx <- diag(svdX$d[1:b])

X2 <- svdX$u %*% Dx %*% t(svdX$v)  # X = U D V'
rownames(X2) <- rownames(X); colnames(X2) <- colnames(X2)
Beta_svd <- crossprod(X2, Y)/n
SigmaHat_svd<-findSigmaHats(X2,Y)
Pvals_svd<-findPvals(Beta_svd, SigmaHat_svd, n)

# TODO add rownames to pvals_svd

# TODO sort and pick top x% of snp-pheno pairs by pval
# (todo later - sort and pick by svd threshold, see how many cross true thresh)

# make dataframe of snp, pheno, pval
Pvals_long<-melt(Pvals_svd, id.vars=1:2)
Pvals_long2<-melt(Pvals_svd)
# TODO use svd of Y or not?
svdY <- svd(Y, nu=c, nv=c)
Dy <- diag(svdY$d[1:c])
Y2 <- svdY$u %*% Dy %*% t(svdY$v)  # Y = U D V'

u1 <- svdX$u
v1 <- Dx %*% t(svdX$v)
u2 <- svdY$u
v2 <- Dy %*% t(svdY$v)