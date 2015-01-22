
senstab = function(x, filt=force) {
 parms = t(sapply(x, "[[", 2))
 tabs = lapply(x, function(z) getTab(z[[1]]))
 counts = sapply(tabs,function(x) approx(x[,2], x[,3], .05)$y)
 counts01 = sapply(tabs,function(x) approx(x[,2], x[,3], .01)$y)
 an = as.numeric
 ans = data.frame(MAF=an(parms[,1]), radius=an(parms[,2]), `FDR 5%`=counts, `FDR 1%`=counts01, check.names=FALSE)
 ans = filt(ans)
 class(ans) = c("senstab", "data.frame")
 ans
}

plot.senstab = function(x, ...) {
mm = melt(x, id.vars=c("MAF", "radius"))
vind = which(names(mm) == "variable")
names(mm)[vind] = "criterion"
ggplot(mm, aes(x=MAF, y=value, shape=criterion)) + geom_point() +
  facet_grid(~radius) + 
  ylab("genes with significant eQTL at given FDR after filtering") + 
  xlab("lower bound on MAF") + ggtitle("cis radius (bp)") +
  theme(axis.text.x=element_text(angle=45))
}
