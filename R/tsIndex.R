
tsIndex.reg = function(tsin, ind) {
#
# obtains the ranges defining genomic tiles in a registry
# comprising a TransStore; metadata on ranges includes
# an environment sn2job giving the map from SNP ids to
# jobs (components of registry)
#
  reg = tsin@allRegistries[[ind]]
  nj = length(findDone(reg))
  smap = new.env(hash=TRUE)
  rngs = sapply(1:nj, function(x) {
       r = loadResult(reg,x)$obs;
       if (is.null(r)) return(NULL);
       ns = names(r)
       Biobase::multiassign(ns, rep(x, length(ns)), smap)
       ans = range(r)
       mcols(ans)$job = x
       ans
       })
  bad = which(sapply(rngs,is.null))
  if (length(bad)>0) ans = do.call(c, rngs[-which(sapply(rngs,is.null))])
    else ans = do.call(c, rngs)
  metadata(ans)$sn2job = smap
  ans
}

