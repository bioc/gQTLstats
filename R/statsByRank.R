# Note: functions with name suffix _sing were implemented for the
# case of a single permutation for plug in FDR
# a TransStore is a list of registries
# each registry is a collection of job outputs
# to traverse the store (possibly in parallel)
# we want efficient direct access to job outputs
# and verification that traversal is complete

.transStoreAccessor_sing = function(ts) {
 rs = allregs(ts)
 js = lapply(rs, findDone) # vector of job names for each
 function(i,j) loadResult( rs[[i]], js[[i]][j] )
}

tsByRank_sing = function(tsin, rank=1, mcol2keep=c("REF", "ALT", "snp", "MAF", "z.HWE"), filt=force) {
 N = nregs(tsin)
 acc = .transStoreAccessor_sing(tsin)
 z = foreach(i = 1:N) %dopar% {
   nj = length(findDone(allregs(tsin)[[i]]))
   lapply(1:nj, function(j) { statsByRank_sing(acc(i, j), rank=rank, filt=filt, mcol2keep=mcol2keep)})
   }
 do.call(c, unlist(z))
}

tsByRankAccum_sing = function(tsin, maxrank=3, mcol2keep=c("REF", "ALT", "snp", "MAF", "z.HWE"), filt=force) {
 N = nregs(tsin)
 acc = .transStoreAccessor_sing(tsin)
 accum = vector("list", maxrank)
 for (cur in 1:maxrank) {
  accum[[cur]] = {
   tmp = foreach(i = 1:N) %dopar% {
    nj = length(findDone(allregs(tsin)[[i]]))
    lapply(1:nj, function(j) { statsByRank_sing(acc(i, j), rank=cur, filt=filt,
         mcol2keep=mcol2keep)})
    }
   do.call(c, unlist(tmp))
   }
  } # now all ranks available
 ans = accum[[1]] # basis 
 tmpmat_scores = matrix(NA, nrow=length(ans), ncol=maxrank)
 tmpmat_permscores = matrix(NA, nrow=length(ans), ncol=maxrank)
 tmpmat_permdists = matrix(NA, nrow=length(ans), ncol=maxrank)
 tmpmat_obsdists = matrix(NA, nrow=length(ans), ncol=maxrank)
 tmpmat_feats = matrix(NA_character_, nrow=length(ans), ncol=maxrank)
 for (i in 1:maxrank) {
   tmpmat_permscores[,i] = accum[[i]]$permscores
   tmpmat_permdists[,i] = accum[[i]]$permdist
   tmpmat_scores[,i] = accum[[i]]$scores
   tmpmat_feats[,i] = accum[[i]]$feats
   tmpmat_obsdists[,i] = accum[[i]]$obsdist
   }
 ans$allpermscores = tmpmat_permscores
 ans$allpermdists = tmpmat_permdists
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


statsByRank_sing = function(job, rank=1, filt=force,
    mcol2keep=c("REF", "ALT", "snp", "MAF", "z.HWE")) {
  jj = job$obs
  if (length(jj)==0) return(NULL)
  mcols(jj)$dist = t(job$dist)  # will choose column later
  inigr = filt(jj)
  pp = job$perm
  mcols(pp)$pdist = t(job$pdist)
  pgr = filt(pp)  # should be same records -- can only filter on SNP features at this stage
  if (class(inigr)!="GRanges"){
    return(NULL) # we can have empty tiles
    }
  stopifnot(class(inigr)=="GRanges")
  ini = mcols(inigr)[,mcol2keep]
  feats = inigr$elnames[,rank]
  scores = inigr$scorebuf[,rank]
  permscores = pgr$scorebuf[,rank]
  obsdist = inigr$dist[,rank]  # dists were transposed
  permdist = pgr$pdist[,rank]
  ini = cbind(ini, DataFrame(feats,scores,permscores,obsdist, permdist))
  mcols(inigr) = ini
  inigr
}


statsByRank = function (job, rank = 1, filt = force, mcol2keep = c("REF", "ALT", 
    "snp", "MAF", "z.HWE")) 
{
    jj = job$obs
    if (length(jj) == 0) 
        return(NULL)
    mcols(jj)$dist = t(job$dist)
    inigr = filt(jj)
    if (class(inigr) != "GRanges") {
        return(NULL)
    }
    ppp = job$perms
    nperm = length(ppp)
    for (i in 1:nperm) {
        pp = ppp[[i]]
        mcols(pp)$pdist = t(job$pdists[[i]])
        if (i==1) pgr = filt(pp)
        else mcols(pgr) = cbind(mcols(pgr), mcols(filt(pp)))
    }
    stopifnot(class(inigr) == "GRanges")
    ini = mcols(inigr)[, mcol2keep]
    feats = inigr$elnames[, rank]
    scores = inigr$scorebuf[, rank]
    obsdist = inigr$dist[, rank]
    permscores = matrix(NA, nrow=length(inigr), ncol=nperm)
    colnames(permscores) = paste0("permscore.", 1:nperm)
    permdists = matrix(NA, nrow=length(inigr), ncol=nperm)
    colnames(permdists) = paste0("pdist.", 1:nperm)
    psinds = grep("scorebuf", names(mcols(pgr)))
    pdinds = grep("pdist", names(mcols(pgr)))
    for (i in 1:nperm) {
       permscores[,i] = mcols(pgr)[, psinds[i]][,rank]
       permdists[,i] = mcols(pgr)[, pdinds[i]][,rank]
       }
    ini = cbind(ini, DataFrame(feats, scores, permscores, obsdist, 
        permdists))
    mcols(inigr) = ini
    inigr
}

 findDone2 = function(x) {
     d = findDone(x)
     if (length(x$availJobs)>0) return(intersect(d, x$availJobs))
     d
   }

.transStoreAccessor2 = function (ts) 
{
    rs = ts@allRegistries
    js = lapply(rs, findDone2)
    function(i, j) loadResult(rs[[i]], js[[i]][j])
}

tsByRankAccum = function (tsin, maxrank = 3, mcol2keep = c("REF", "ALT", "snp", 
    "MAF", "z.HWE"), filt = force) 
{
    N = nregs(tsin)
    acc = .transStoreAccessor2(tsin)
    accum = vector("list", maxrank)
    for (cur in 1:maxrank) {
        accum[[cur]] = {
            tmp = foreach(i = 1:N) %dopar% {
                nj = length(findDone2(allregs(tsin)[[i]]))
                lapply(1:nj, function(j) {
                  statsByRank(acc(i, j), rank = cur, filt = filt, 
                    mcol2keep = mcol2keep)
                })
            }
            do.call(c, unlist(tmp))
        }
    }
    ini = accum[[1]][, mcol2keep]
    alldists = matrix(NA, nrow=length(ini), ncol=maxrank)
    allscores = matrix(NA, nrow=length(ini), ncol=maxrank)
    allfeats = matrix(NA_character_, nrow=length(ini), ncol=maxrank)
    for (i in 1:maxrank) {
      allscores[,i] = accum[[i]]$scores
      alldists[,i] = accum[[i]]$obsdist
      allfeats[,i] = accum[[i]]$feats
      }
    nperms = length(pinds <- grep("permscore", names(mcols(accum[[1]]))))
    pdinds = grep("pdist", names(mcols(accum[[1]])))
    permscoresByRank = vector("list", maxrank)
    permdistsByRank = vector("list", maxrank)
    for (i in 1:maxrank) {
       permscoresByRank[[i]] = mcols(accum[[i]])[, pinds]
       permdistsByRank[[i]] = mcols(accum[[i]])[, pdinds]
    }
    for (i in c("allscores", "alldists", "allfeats"))
      mcols(ini)[[i]] = get(i)
    permscoreNames = paste0("permscoresByRank", 1:maxrank)
    permdistNames = paste0("pdistsByRank", 1:maxrank)
    for (i in 1:maxrank) {
      mcols(ini)[[permscoreNames[i]]] = permscoresByRank[[i]]
      mcols(ini)[[permdistNames[i]]] = permdistsByRank[[i]]
    }
    ini
}

