# say there are n indiv, m snp, k probes

# G, data frame of 0/1/2 genotypes, n x m
# X, matrix of normalized genotypes, n x m
# Y, matrix of normalized phenotypes, n x k


# first 1k snps are all from chr 1
# chr 1 has 2k genes, chr 2 has 3k so hopefully result within first 5k 

snps.txt.file<-"GTEx_data/Lung1k.snps.txt" 
expr.txt.file<-"GTEx_data/Lung5k.expr.txt" 

Gt<-read.table(snps.txt.file, header=TRUE, sep="", row.names=1)
Yt.df<-read.table(expr.txt.file, header=TRUE, sep="", row.names=1)

G<-t(Gt)
X<-scale(G)  # standardize genotypes to mean 0, var 1
Y<-t(data.matrix(Yt.df))

# X^T = m x n, Y = n x k

dim(X)
dim(Y)

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