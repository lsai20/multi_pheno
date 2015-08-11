# say there are n indiv, m snp, k probes. for real data, m >> k



n<-nrow(X); m<-ncol(X); k<-ncol(Y)

# X^T = m x n, Y = n x k
### instead of cross prod, to get Beta, use SVD ###
#Beta <- crossprod(X,Y)/nrow(X)  # m x k matrix of betas
svdX <- svd(X, nu=110, nv=110)
#dim(svdX$u); dim(diag(svdX$d)); dim( t(svdX$v))
#Dx <- diag(svdX$d)
Dx <- t(svdX$u) %*% X %*% svdX$v # D = U' X V
X2 <- svdX$u %*% Dx %*% t(svdX$v) #  X = U D V'
rownames(X2)<-rownames(X2); colnames(X2)<-colnames(X)
Beta_svd <- crossprod(X2, Y)/n
Beta_true <- crossprod(X,Y)/n

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
