

if (FALSE){
# svd example
X <- replicate(10000, rnorm(2000)) 
s <- svd(X)
D <- diag(s$d)
X2 <- s$u %*% D %*% t(s$v) #  X = U D V'
D2 <- t(s$u) %*% X %*% s$v #  D = U' X V
}


# Compare timing of svd and X^T*Y
n <- 1000; m <- 5000; k <- 2000
X <- replicate(n, rnorm(m))
Y <- replicate(n, rnorm(k))

# SVD
t0 <- proc.time()
s <- svd(X)
t1 <- proc.time()
svdTime <- t1-t0


# X^T * Y (cross prod)
t0 <- proc.time()
xty <- crossprod(X,Y)
t1 <- proc.time()
xtyTime <- t1-t0


# X^T * Y (tranpose X, then times Y)
t0 <- proc.time()
xty2 <- t(X) %*% Y
t1 <- proc.time()
xty2Time <- t1-t0
