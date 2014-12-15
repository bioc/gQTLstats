#plot.sensobj.thresh = function(sensobj, thresh=.05,
#   returnDF=FALSE, filter = force ) {
#  require(gQTLstats)  # for getTab etc
#  ncalls = sapply(sensobj, 
#     function(x) approx(getTab(x[[1]])$fdr, 
#       getTab(x[[1]])$ncalls, thresh)$y)
#  mafs = sapply(sensobj, function(x) as.numeric(x$parms[1]))
#  dists = sapply(sensobj, function(x) as.numeric(x$parms[2]))
#  dfth = data.frame(calls=ncalls, maf=mafs, dist=dists)
#  if (returnDF) return(dfth)
#  ggplot(data=dfth, aes(y=call, x=maf)) + facet_grid(.~dist) +
#     ggtitle("cis radius (bp)") + ylab(paste0(
#         "# probes with FDR < ", thresh)) + 
#         xlab("lower bound on MAF")
#}

table.sensobj.thresh = function(sensobj, thresh=.05) {
  require(gQTLstats)  # for getTab etc
  ncalls = sapply(sensobj, 
     function(x) approx(getTab(x[[1]])$fdr, 
       getTab(x[[1]])$ncalls, thresh)$y)
  mafs = sapply(sensobj, function(x) as.numeric(x$parms[1]))
  dists = sapply(sensobj, function(x) as.numeric(x$parms[2]))
  dfth = data.frame(calls=ncalls, maf=mafs, dist=dists, thresh=thresh)
  class(dfth) = c("table.sensobj", "data.frame")
  dfth
}

plot.table.sensobj = function(x, y, ...) {
  x = data.frame(unclass(x))  # :)
  nthresh = length(ut <- unique(x$thresh))
  if (nthresh == 1) {
      gf = facet_grid(.~dist)
      yl = ylab(paste0("# probes with FDR < ", ut)) 
      }
  else {
      gf = facet_grid(thresh~dist)
      yl = ylab(paste0("# probes with FDR < given threshold"))
      }
  ggplot(data=x, aes(x=maf, y=calls)) + geom_point() + gf +
     ggtitle("cis radius (bp)") + yl +
         xlab("lower bound on MAF") +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


