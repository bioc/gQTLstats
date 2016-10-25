# a TransStore is a list of registries
# each registry is a collection of job outputs
# to traverse the store (possibly in parallel)
# we want efficient direct access to job outputs
# and verification that traversal is complete

.transStoreAccessor = function(ts) {
 rs = allregs(ts)
 js = lapply(rs, findDone) # vector of job names for each
 function(i,j) loadResult( rs[[i]], js[[i]][j] )
}

tsByRank = function(ts, rank=1) {
 N = nregs(ts)
 acc = .transStoreAccessor(ts)
 z = foreach(i = 1:N) %dopar% {
   nj = length(findDone(allregs(ts)[[i]]))
   lapply(1:nj, function(j) { statsByRank(acc(i, j), rank)})
   }
 do.call(c, unlist(z))
}

allregs = function(tstore) tstore@allRegistries
nregs = function(tstore) length(allregs(tstore))
nthreg = function(n, tstore) {
  stopifnot(n >= 1)
  stopifnot(n <= nregs(tstore))
  allregs(tstore)[[n]]
}

#statsByRank = function(tstore, rank) {
#  probe1 = loadResults(allregs(tstore)[[1]])[[1]]
#  nregs = length(allregs(tstore))
#  therank = nrow(probe1$dist)
#  stopifnot(rank <= therank)
#  inigr = probe1$obs
#  ini = mcols(inigr)[,1:4]
#  feats = inigr$elnames[,rank]
#  scores = inigr$scorebuf[,rank]
#  permscores = probe1$perm$scorebuf[,rank]
#  ini = cbind(ini, DataFrame(feats,scores,permscores))
#  mcols(inigr) = ini
#  inigr
#}


statsByRank = function(job, rank=1) {
  inigr = job$obs
  if (class(inigr)!="GRanges"){
    return(NULL) # we can have empty tiles
    }
  stopifnot(class(inigr)=="GRanges")
  ini = mcols(inigr)[,1:4]
  feats = inigr$elnames[,rank]
  scores = inigr$scorebuf[,rank]
  permscores = job$perm$scorebuf[,rank]
  ini = cbind(ini, DataFrame(feats,scores,permscores))
  mcols(inigr) = ini
  inigr
}
