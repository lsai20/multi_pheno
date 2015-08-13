# computes/approximates Beta using SVD
# see beta_setup.R for additional fxns/reading input

snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung30.expr.txt" 

X<-importX(snps.txt.file)
Y<-importY(expr.txt.file)

# say there are n indiv, m snp, k probes. for real data, m >> k
n<-nrow(X); m<-ncol(X); k<-ncol(Y)

# say we use b components to approximate X (exact is b = min(n))
# and c components to approximate Y
b<-10
c<-0

### instead of cross prod, to get Beta, use SVD ###
svdX <- svd(X, nu=w, nv=w) #  X = U D V'
Dx <- diag(svdX$d[1:w])
dim(svdX$u); dim(Dx); dim(t(svdX$v))

X2 <- svdX$u %*% Dx %*% t(svdX$v)
Beta_svd <- crossprod(X2, Y)/n
SigmaHat_svd<-findSigmaHats(X2,Y)
Pvals_svd<-findPvals(Beta_svd, SigmaHat, n)


# TODO use svd of Y or not?
svdY <- svd(Y, nu=7, nv=9)
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
