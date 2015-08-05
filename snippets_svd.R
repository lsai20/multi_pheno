# uses some functions from compute_XT_Y.R (import manually for now)

if (FALSE){
  # svd example
  X <- replicate(10000, rnorm(2000)) 
  s <- svd(X)
  D <- diag(s$d)
  X2 <- s$u %*% D %*% t(s$v) #  X = U D V'
  D2 <- t(s$u) %*% X %*% s$v #  D = U' X V
}


# Compare timing of svd and X^T*Y
# Use a very rectangular matrix, typical for GE data
n <- 10; m <- 5000; k <- 500
X <- replicate(m, rnorm(n))  # n x m
Y <- replicate(k, rnorm(n))  # n x k

# SVD, default
t0 <- proc.time()
s <- svd(X)
t1 <- proc.time()
svdTime <- t1-t0

#SVD, 1 component
t0 <- proc.time()
s1 <- svd(X, nu=1, nv=1)
t1 <- proc.time()
svd1Time <- t1-t0

# SVD, full
t0 <- proc.time()
s <- svd(X, nu=1000, nv=5000)
t1 <- proc.time()
svdFullTime <- t1-t0

# SVD, default, transpose
Xt <- t(X)
t0 <- proc.time()
s <- svd(Xt)
t1 <- proc.time()
svdTransTime <- t1-t0

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


print(svdTime)
print(svd1Time)
print(svdFullTime)
print(svdTransTime)
print(xtyTime)
print(xty2Time)
