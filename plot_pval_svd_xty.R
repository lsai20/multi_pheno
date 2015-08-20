# plotting

# scatter plot, #000000 with 0x33 = 3/16 opacity
plot(topTrue$pval_true, topTrue$pval_svd, col="#00000033")

library(ggplot2)
library(hexbin)

# hex bin plots
ggplot(Results_all,aes(x=pval_true,y=pval_svd)) + stat_binhex()
plot(hexbin(Results_all$pval_true, Results_all$pval_svd))
plot(hexbin(topTrue$pval_true, topTrue$pval_svd))

# scatter plots for top 1% pval
ggplot(topTrue,aes(x=pval_true,y=pval_svd)) + geom_point(alpha = 0.2) +
  theme_light() + ggtitle("SVD vs true pval for top 1% true p-value") + coord_fixed(ratio=1)

ggplot(topTrue,aes(x=percentile_pval_true,y=percentile_pval_svd)) +
  geom_point(alpha = 0.2) + theme_light() + 
  ggtitle("SVD vs true percentile for top 1% true p-value") + coord_fixed(ratio=1)

# smoothed scatter for all snp-pheno pairs, base R
smoothScatter(topTrue$pval_true, topTrue$pval_svd, nbin=500, nrpoints=500)

# translucent scatter for all snp-pheno pairs, ggplot
ggplot(Results_all,aes(x=pval_true,y=pval_svd)) + geom_point(alpha = 0.1) + theme_light() + ggtitle("SVD vs true pval for all snp-pheno pairs")

# smooth scatters with legend
#http://stackoverflow.com/questions/14271584/
# r-legend-for-color-density-scatterplot-
# produced-using-smoothscatter
library(fields)
smoothScatterLegend <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)  


oldpar <- par()
par(mar = c(5,4,4,5.5) + .1)

smoothScatter(topTrue$pval_true, topTrue$pval_svd, nrpoints=500, 
              postPlotHook=smoothScatterLegend)
title("svdX100 smoothScatter( ... nrpoints=500, 
      postPlotHook=smoothScatterLegend)")
dev.off()

smoothScatter(Results_all$pval_true, Results_all$pval_svd, nrpoints=500, 
              postPlotHook=smoothScatterLegend)
title("svdX119 smoothScatter(..., nrpoints=500, 
      postPlotHook=smoothScatterLegend)")

#dev.off()
