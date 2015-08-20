# for varying rank, find 99% RR thresh


RR_target = 0.99    # target recall rate
thresh_true = 10^-6 # desired threshold
rank_X = 100        # rank of svdX (exact is 119)

# TODO 
# given svd X of particular rank,
# given pval_svd and pval_true for each snp-pheno pair,
# find 99% thresh (i.e. thresh for pval_svd st. pval_svd < thresh contains 99% of true)

# eventually for loop for rank_X ranging 1 to 119

 