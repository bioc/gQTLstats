#library(gQTLstats)
# commented code can go into unit test
#if (!exists("result")) data("1-result")
setClass("TransChunk", representation(
  obs="GRanges", perm="GRanges", obsDist="matrix",
     permDist="matrix", bufsize="numeric"))
setMethod("show", "TransChunk", function(object) {
 cat("TransChunk instance for", length(object@obs), "SNPs.\nHead of observed results:\n")
 show(head(object@obs))
})
TransChunk = function(x) {
 nx = names(x)
 stopifnot(all(c("obs", "perm", "dist", "pdist")==nx))
 new("TransChunk", obs=x$obs, perm=x$perm,
    obsDist=x$dist, permDist=x$pdist, bufsize=nrow(x$dist))
}
#t1 = TransChunk(result)
ScoreDF = function(tc, use="obs") {
 bufsize = tc@bufsize
 inUse = slot(tc, use)
 distToUse = ifelse(use=="obs", "obsDist", "permDist")
 nsnp = length(names(inUse))
 rank = rep(1:bufsize, nsnp)
 snpid = rep(names(inUse), each=bufsize)
 loc = rep(start(inUse), each=bufsize)
 chr = rep(as.character(seqnames(inUse)), each=bufsize)
 MAF = rep(mcols(inUse)$MAF, each=bufsize)
 scores = as.numeric(t(inUse$scorebuf))
 probes = as.character(t(inUse$elnames))
 DataFrame(snpid=snpid, chr=chr, loc=loc, probes=probes, rank=rank, 
   scores=scores,
   MAF=MAF, dist=as.numeric(slot(tc, distToUse)))
}
filteredDFwPerm = function(tc, filter=force) {
  odf = filter(ScoreDF(tc, "obs"))
  permdf = filter(ScoreDF(tc, "perm"))
  metadata(odf)$filteredPermScores = permdf$scores
  odf
}
#ft1 = filteredDFwPerm(t1, filter=function(x) x[x$MAF > 0.05 & x$dist >= Inf & x$rank==1,])
#  


transTable = function( reg, filter = force, jobs ) {
#
# this code will use the BatchJobs registry reg to
# collect information on SNP-feature associations, filtered
# according to the function in filter (which operates on a DataFrame
# and returns one)
#
# a DataFrame is returned with plug-in FDR according to the
# filtering scheme
#
  if (missing(jobs)) jobs = findDone(reg)
  lens = foreach(i = jobs, .combine=c) %dopar% length(loadResult(reg, i))
  bad = which(lens==0)
  if (length(bad)>0) jobs = jobs[-bad]
  objs = foreach(i = jobs) %dopar% filteredDFwPerm( TransChunk(loadResult(reg, i)), 
                                     filter=filter )
  allperms = unlist(lapply(objs, function(x) metadata(x)[["filteredPermScores"]]))
  fullDF = do.call(rbind, objs)
  nperm = length(allperms)
  nobs = nrow(fullDF)
  needrepl = (nperm < nobs)
  if (nperm != nobs) allperms = allperms[ sample(1:nperm, size=nobs, replace = needrepl) ]
  fullDF$pifdr = pifdr(fullDF$scores, allperms)
  fullDF$permscores = allperms
  metadata(fullDF)$filterfun = filter
  fullDF
}


setClass("TransStore", representation(allRegistries="list",
  numSubmitted="numeric", numDone="numeric", nloci="numeric",
   jobinfos="list"))
setGeneric("describe", function(object)standardGeneric("describe"))
setMethod("describe", "TransStore", function(object) {
  ids = sapply(object@allRegistries, function(z) z$id)
  fracs = paste("(", paste(object@numDone, object@numSubmitted, sep="/"), ")", sep="")
  maxmem = sapply(object@jobinfos, function(x) max(x$memory, na.rm=TRUE))
  maxtime = sapply(object@jobinfos, function(x) max(x$time.running, na.rm=TRUE))
  DataFrame(regids=ids, fracDone=fracs, nloci=object@nloci,
     maxmem=maxmem, maxtime=maxtime)
  })

TransStore = function(regs, paths=NULL) {
   requireNamespace("doParallel")
   if (!is.null(paths)) regs = lapply(paths, loadRegistry)
   findJobs2 = function(x) {
     if (length(x$availJobs)>0) return(x$availJobs)
     findJobs(x)
   }
   findDone2 = function(x) {
     d = findDone(x)
     if (length(x$availJobs)>0) return(intersect(d, x$availJobs))
     d
   }
   getJobInfo2 = function(x) {
     d = findDone2(x)
     getJobInfo(x)[as.character(d),]
   }
   subs = sapply(regs, function(x) length(findJobs2(x)))
   dones = sapply(regs, function(x) length(findDone2(x)))
   todrop = which(dones==0)
   if (length(todrop)>0) {
      message(paste(length(todrop), "registries have no completed jobs, dropped"))
      regs = regs[-todrop]
      subs = subs[-todrop]
      dones = dones[-todrop]
      }
   lens = foreach(i=1:length(regs)) %dopar%
      sum(sapply(findDone2(regs[[i]]), function(x)length(loadResult(regs[[i]],x)[[1]])))  # first element is locus-centric GRanges
   jis = lapply(regs, function(x) getJobInfo2(x))
   new("TransStore", allRegistries=regs, nloci=unlist(lens),
      numSubmitted=subs, numDone=dones, jobinfos=jis)
}
setMethod("show", "TransStore", function(object) {
  cat("TransStore instance with ", length(object@allRegistries),
 " registries.  Description:\n", sep="")
  show(describe(object))
  })

