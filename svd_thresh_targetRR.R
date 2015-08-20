# for varying rank, find 99% RR thresh


RR_target = 0.99    # target recall rate
thresh_true = 10^-6 # desired threshold
rank_X = 100        # rank of svdX (exact is 119)

# TODO 
# given svd X of particular rank,
# given pval_svd and pval_true for each snp-pheno pair,
# find 99% thresh (i.e. thresh for pval_svd st. pval_svd < thresh contains 99% of true)

# eventually for loop for rank_X ranging 1 to 119
#want 99% of true pvals to have p_svd within thresh_svd
#any 99% of true pvals ok. (i.e if 100 true vals, okay to get any 99, not top 99)

# TODO
# test on small data set - 1k snp, 30 pheno

topTrue <- Results_all[which(Results_all$pval_true <= thresh_true),]
topTrue <- topTrue[order(topTrue$pval_svd),]

cor(Results_all$pval_svd, Results_all$pval_true)
cor(topTrue$pval_svd, topTrue$pval_true)

