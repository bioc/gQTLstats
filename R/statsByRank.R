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

tsByRank = function(ts, rank=1, filt=force) {
 N = nregs(ts)
 acc = .transStoreAccessor(ts)
 z = foreach(i = 1:N) %dopar% {
   nj = length(findDone(allregs(ts)[[i]]))
   lapply(1:nj, function(j) { statsByRank(acc(i, j), rank, filt=filt)})
   }
 do.call(c, unlist(z))
}

tsByRankAccum = function(ts, maxrank=3, filt=force) {
 N = nregs(ts)
 acc = .transStoreAccessor(ts)
 accum = vector("list", maxrank)
 for (cur in 1:maxrank) {
  accum[[cur]] = {
   tmp = foreach(i = 1:N) %dopar% {
    nj = length(findDone(allregs(ts)[[i]]))
    lapply(1:nj, function(j) { statsByRank(acc(i, j), cur, filt=filt)})
    }
   do.call(c, unlist(tmp))
   }
  } # now all ranks available
 ans = accum[[1]] # basis 
 tmpmat_scores = matrix(NA, nr=length(ans), nc=maxrank)
 tmpmat_permscores = matrix(NA, nr=length(ans), nc=maxrank)
 tmpmat_obsdists = matrix(NA, nr=length(ans), nc=maxrank)
 tmpmat_feats = matrix(NA_character_, nr=length(ans), nc=maxrank)
 for (i in 1:maxrank) {
   tmpmat_permscores[,i] = accum[[i]]$permscores
   tmpmat_scores[,i] = accum[[i]]$scores
   tmpmat_feats[,i] = accum[[i]]$feats
   tmpmat_obsdists[,i] = accum[[i]]$obsdist
   }
 ans$allpermscores = tmpmat_permscores
 ans$allscores = tmpmat_scores
 ans$alldists = tmpmat_obsdists
 ans$allfeats = tmpmat_feats
 ans
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


statsByRank = function(job, rank=1, filt=force) {
  jj = job$obs
  if (length(jj)==0) return(NULL)
  mcols(jj)$dist = t(job$dist)  # will choose column later
  inigr = filt(jj)
  pgr = filt(job$perm)  # should be same records
  if (class(inigr)!="GRanges"){
    return(NULL) # we can have empty tiles
    }
  stopifnot(class(inigr)=="GRanges")
  ini = mcols(inigr)[,1:4]
  feats = inigr$elnames[,rank]
  scores = inigr$scorebuf[,rank]
  permscores = pgr$scorebuf[,rank]
  obsdist = inigr$dist[,rank]  # dists were transposed
  #permdist = inigr$pdist[rank,]
  ini = cbind(ini, DataFrame(feats,scores,permscores,obsdist))#, permdist))
  mcols(inigr) = ini
  inigr
}
