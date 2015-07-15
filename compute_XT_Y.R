# say there are n indiv, m snp, k probes

# G, data frame of 0/1/2 genotypes, m x n
# X, matrix of normalized genotypes, m x n
# Y, matrix of normalized phenotypes, k x n

snps.txt.file<-"Lung.snps.txt"
expr.txt.file<-"Lung.expr.txt"
G<-read.table(snps.txt.file, header=TRUE, sep=" ", row.names=1)
Y.df <-read.table(expr.txt.file, header=TRUE, sep=" ", row.names=1)

# TODO matrixify/normalize G
# TODO convert Y.df to matrix, keep probe labels

if (false) {
# A = U S V^T
# U, columns are eigenvectors of A*A^T
# S, Sigma, diagonal matrix of eigenvalues
# V, columns are eigenvectors of A^T*A
}