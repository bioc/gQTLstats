

vecToQuantiles = function( vec,
     probs=c(seq(0,.999,.001), 1-(c(1e-4,1e-5,1e-6,1e-7))), 
     ids=NULL, ..., checkField=FALSE, filter=force ) {
# emulates store-based interface
#
  quantile(filter(vec), probs )
}

#vecToHist = function( vec, breaks, filter=force)
#   hist(vec, breaks=breaks)
#
#vecsToFDR = function(vec, permvec, xprobs = c(seq(0, 0.999, 0.001), 1 - (c(1e-04,
#    1e-05, 1e-06, 1e-07))), xfield = "chisq",
#    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
#       grep("permScore", names(mcols(x)))])),
#       filter=force, .id4coln=1, ids=NULL) {
# theCall = match.call()
# ntests = length(vec)
# message("counting #NA...")
# nna = sum(is.na(vec))
# ntests = ntests - nna
# message("obtaining assoc quantiles...")
# xq = vecToQuantiles(vec, probs=xprobs, filter=filter, ids=ids) # nxq
# message("computing perm_assoc histogram....")
# yh.orig = vecToHist(permvec, getter=getter, breaks=c(0,xq,1e10), filter=filter, ids=ids)$counts # nxq+2
# yh = numeric(length(yh.orig)-1)
# yh[1] = yh.orig[1] + yh.orig[2]  # fuse first 2 elements for left boundary
# yh[2:length(yh)] = yh.orig[-c(1,2)] # transfer
# oy = nperm*ntests - cumsum(yh) # how many perm scores exceed cuts in assoc score
# ncalls = ntests*(1-xprobs)
# trimToUnit = function(x) pmax(0, pmin(1, x))
# fdr = trimToUnit(oy/(nperm*ncalls))
# ans = data.frame(assoc=xq, fdr=fdr, ncalls=ncalls, avg.false=oy/nperm)
# new("FDRsupp", tab=ans, theCall=theCall, sessinfo=sessionInfo(),
#    filterUsed=filter)
#}
#
#

storeToQuantiles = function( store, field, 
     probs=c(seq(0,.999,.001), 1-(c(1e-4,1e-5,1e-6,1e-7))), 
     ids=NULL, ..., checkField=FALSE, filter=force ) {
#
# uses ff for interim storage of data for summarization
#
  src = storeToFf( store, field, ids=ids, checkField=checkField, filter=filter )
  quantile(src, probs )
}

storeToHist = function( store, getter =
   function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[,grep("permScore",
                         names(mcols(x)))])),
            breaks, ids=NULL, filter=force ) {
   if (missing(breaks)) stop("breaks must be supplied and must cover range of data")
   if (is.null(ids)) ids=store@validJobs
##BP   tmp = bplapply(ids, function(x) {
   tmp = foreach(x=ids) %dopar% { #bplapply(ids, function(x) {
      getter(loadAndFilterResult(reg=store@reg, id=x, filter=filter))
##BP      })
   }
##BP   tmp = bplapply(tmp, function(x) try(hist(x , breaks=breaks, plot=FALSE )))
   tmp = foreach(x=tmp) %dopar% {
         try(hist(x , breaks=breaks, plot=FALSE ))
   }
     # can fail with NAs in x
   st = sapply(tmp, class)
   if (any(st == "try-error")) {
      warning("had an error in histogram computation, returning all attempts")
      return(tmp)
      }
   ans = tmp[[1]]
   if (length(tmp)==1) return(ans)
   for (i in 2:length(tmp)) ans$counts = ans$counts + tmp[[i]]$counts
   ans
}
  

storeToFDR = function(store, xprobs = c(seq(0, 0.999, 0.001), 1 - (c(1e-04,
    1e-05, 1e-06, 1e-07))), xfield = "chisq",
    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
       grep("permScore", names(mcols(x)))])),
       filter=force, .id4coln=1, ids=NULL) {
 theCall = match.call()
 r1 = storeApply(store, force, ids=.id4coln, flatten1=TRUE)[[1]]
 nperm.potential = length(grep("permScore", names(mcols(r1))))
 pgotten = getter(r1)  # need to explicitly determine number of perms selected
           # as it need not be the full number in store
 nperm = length(pgotten)/length(r1)
 stopifnot(nperm>0)
 message("counting tests...")
 ntests = sum(unlist(storeApply( store, function(x) length(filter(x)), ids=ids )) , na.rm=TRUE)
 message("counting #NA...")
 nna = sum(unlist(storeApply( store, function(x)sum(is.na(mcols(filter(x))[,xfield])), ids=ids )) , na.rm=TRUE)
 ntests = ntests - nna
 message("obtaining assoc quantiles...")
 xq = storeToQuantiles(store, field=xfield, probs=xprobs, filter=filter, ids=ids) # nxq
 message("computing perm_assoc histogram....")
 yh.orig = storeToHist(store, getter=getter, breaks=c(0,xq,1e10), filter=filter, ids=ids)$counts # nxq+2
 yh = numeric(length(yh.orig)-1)
 yh[1] = yh.orig[1] + yh.orig[2]  # fuse first 2 elements for left boundary
 yh[2:length(yh)] = yh.orig[-c(1,2)] # transfer
 oy = nperm*ntests - cumsum(yh) # how many perm scores exceed cuts in assoc score
 ncalls = ntests*(1-xprobs)
 trimToUnit = function(x) pmax(0, pmin(1, x))
 fdr = trimToUnit(oy/(nperm*ncalls))
 ans = data.frame(assoc=xq, fdr=fdr, ncalls=ncalls, avg.false=oy/nperm)
 new("FDRsupp", tab=ans, FDRfunc=NULL, FDRmodel=NULL, theCall=theCall, sessinfo=sessionInfo(),
    filterUsed=filter)
}

addFDRfunc = function(supp, f) {
  stopifnot(is(supp, "FDRsupp"))
  supp@FDRfunc = f
  supp
}
addFDRmodel = function(supp, m) {
  stopifnot(is(supp, "FDRsupp"))
  supp@FDRmodel = m
  supp
}
 

maxByProbeOLD = function (job, res, filter)
{
    as(filter(res), "data.frame") %>% dplyr::select(snp, probeid,
           chisq, permScore_1, permScore_2, permScore_3,
           seqnames, start) %>% group_by(probeid) %>%
      summarize(chisq=max(chisq), permScore_1=max(permScore_1),
            permScore_2=max(permScore_2), permScore_3=max(permScore_3))
}

maxByFeature = function (job, res, resfilter, feature="SNP")
{
#
# a problem is that a given SNP might reside in multiple jobs?
#
# currently assumes a field 'chisq' holds the quantity of interest
#
# modifications 1 Oct 2015 -- address variable numbers of permutations
#
#
# two helper functions that naively assist in the creation of
# unpredictable select predicates.  we precompute to cater for
# up to 10 permutations.  probably will not need more
#
# the idea here is that we get the largest observed per-gene score
# statistic, and the largest per-gene statistic for each permutation
#
gennl = function(in1 =
    list(~snp, ~probeid, ~MAF, ~chisq, ~mindist, ~seqnames, 
             ~start, ~permScore_1, ~permScore_2, ~permScore_3, 
             ~permScore_4, ~permScore_5, ~permScore_6, ~permScore_7,
             ~permScore_8, ~permScore_9, ~permScore_10), permn) {
 if (permn > 10) stop("doesn't handle > 10 perms ... simple fix available")
 if (permn < 1) stop("doesn't handle < 1 perms")
 in1[1:(7+permn)]
}

genmaxl = function(n) {
 l1 = list(permScore_1=~max(permScore_1), permScore_2=~max(permScore_2), 
           permScore_3=~max(permScore_3), permScore_4=~max(permScore_4),
           permScore_5=~max(permScore_5), permScore_6=~max(permScore_6),
           permScore_7=~max(permScore_7), permScore_8=~max(permScore_8),
           permScore_9=~max(permScore_9), permScore_10=~max(permScore_10))[1:n]
 addit = list(chisq=~max(chisq), seqnames=~nth(seqnames,1), start=~nth(start,1),
               MAF=~nth(MAF,1), mindist=~nth(mindist,1))
 c(l1, addit)
}

    permnms = grep("permScore", names(mcols(res)), value=TRUE)
    nperm = length(grep("permScore", names(mcols(res))))

ndots = gennl(permn=nperm)
maxdots = genmaxl(n=nperm)

    pregroup = as(resfilter(res), "data.frame") %>% dplyr::select_(
            .dots=ndots) 
    if (feature=="SNP") {
        postg = pregroup %>% group_by(snp) %>% arrange(desc(chisq),snp) 
        maxdots = c(maxdots, list(probeid=~nth(probeid,1)))
        }
    else if (feature=="probeid") {
        postg = pregroup %>% group_by(probeid) %>% arrange(desc(chisq),probeid)
        maxdots = c(maxdots, list(snp=~nth(snp,1)))
        }
    postg %>% dplyr::summarise_(.dots=maxdots)
}

maxBySNP = function(job, res, resfilter)
  maxByFeature(job, res, resfilter, feature="SNP")
maxByProbe = function(job, res, resfilter)
  maxByFeature(job, res, resfilter, feature="probeid")

enumerateByFDR = function (store, fdrsupp, threshold = 0.05,
    filter=force, ids=NULL, trimToUnit=TRUE) 
{
    if (!is.function(getFDRfunc(fdrsupp))) {
       message("error: fdrsupp has no FDR interpolating function.")
       message("please add one using setFDRfunc().")
       stop("cannot find FDRfunc element in fdrsupp.")
       }
    trimmer = force
    if (trimToUnit) trimmer = function(x) pmax(0,pmin(1,x))
    selector = function(x) {
        x$estFDR = getFDRfunc(fdrsupp)(x$chisq)
        x$estFDR = trimmer(x$estFDR)
        x[which(x$estFDR <= threshold)]
    }
    ans = storeApply(store, function(x) filter(selector(x)), ids=ids)
#
# it is possible that an empty range will result, if no pair
# has sufficiently low FDR
#
#    ul = ans  # revision to storeApply that flattens should remove need for this
#    zl = sapply(ul,length)
#    if (any(zl==0)) ul=ul[-which(zl==0)]
    ans = unlist(GRangesList(ans))
    metadata(ans)$enumCall = match.call()
    metadata(ans)$enumSess = sessionInfo()
    metadata(ans)$fdrCall = fdrsupp@theCall
    ans
}

dfrToFDR = function(dfr, xprobs = c(seq(0, 0.999, 0.001), 1 - (c(1e-04,
    1e-05, 1e-06, 1e-07))), xfield = "chisq",
    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
       grep("permScore", names(mcols(x)))])), filter=force) {
#
# when in-memory is OK
#
# data.frame method
 theCall = match.call()
 dfr = filter(dfr)
 stopifnot("chisq" %in% names(dfr))
 nperm = length(grep("permScore", names(dfr)))
 message("counting tests...")
 ntests = sum(!is.na(dfr$chisq))
 message("counting #NA...")
 nna = sum(is.na(dfr$chisq))
 ntests = ntests - nna
 message("obtaining assoc quantiles...")
 xq = quantile(dfr$chisq, probs=xprobs) # storeToQuantiles(store, field=xfield, probs=xprobs, filter=filter) # nxq
 message("computing perm_assoc histogram....")
 yh.orig = hist( c(dfr$permScore_1, dfr$permScore_2, dfr$permScore_3), breaks=c(0,xq,1e10), plot=FALSE)$counts   #storeToHist(store, getter=getter, breaks=c(0,xq,1e10), filter=filter)$counts # nxq+2
 yh = numeric(length(yh.orig)-1)
 yh[1] = yh.orig[1] + yh.orig[2]  # fuse first 2 elements for left boundary
 yh[2:length(yh)] = yh.orig[-c(1,2)] # transfer
 oy = nperm*ntests - cumsum(yh) # how many perm scores exceed cuts in assoc score
 ncalls = ntests*(1-xprobs)
 trimToUnit = function(x) pmax(0, pmin(1, x))
 fdr = trimToUnit(oy/(nperm*ncalls))
 ans = data.frame(assoc=xq, fdr=fdr, ncalls=ncalls, avg.false=oy/nperm)
 new("FDRsupp", tab=ans, theCall=theCall, sessinfo=sessionInfo(),
    filterUsed=filter)
}

# we use dplyr here for group_by performance
#
storeToFDRByProbe = function( store, xprobs = seq(0, 0.99, 0.01),
       xfield = "chisq", filter=force, ids=NULL) {
#    maxbp = unlist(storeApply(store, function(x) {
#       maxByProbe(1, filter(x), force) }, ids=ids), recursive=FALSE)
    maxbp = storeApply(store, function(x) {
       maxByProbe(1, filter(x), force) }, ids=ids)
suppressWarnings({  # expected unequal factor levels
    maxbp = bind_rows(maxbp)  # note filter was already applied
})
    dfrToFDR(maxbp, xprobs = xprobs, xfield = xfield, filter=filter)
}

storeToMaxAssocBySNP = function( store, chr=1, xprobs = seq(0, 0.999, 0.001),
    xfield = "chisq", nperm=3, resfilter=force) {
    stopifnot(is.atomic(chr) & length(chr)==1)
    if (length(store@rangeMap)==0) stop("please use a store with valid rangeMap component.")
    rmap = store@rangeMap
    ac = as.character
    jobsToDo = as.integer(rmap[ which(ac(seqnames(rmap))==ac(chr)) ]$jobid)
    stopifnot(length(jobsToDo)>0)
    #maxbsByJob = unlist(storeApply(store, function(x) {  # new storeApply defaults
    #   maxBySNP(1, resfilter(x), force) }, ids=jobsToDo ), recursive=FALSE)
    maxbsByJob = storeApply(store, function(x) {
       maxBySNP(1, resfilter(x), force) }, ids=jobsToDo )
    message("binding all jobs")  # SNPs may occupy different jobs so another agg
    aggr = do.call(rbind, maxbsByJob)
    message("done.")
    message("group by and max")
    # after collecting, use
    #dfrToFDR( aggr, xprobs = xprobs, xfield=xfield )
    ans <- aggr %>% group_by(snp) %>% arrange(desc(chisq),snp) %>%
               summarize( chisq = max(chisq),
                 permScore_1 = max(permScore_1),
                 permScore_2 = max(permScore_2),
                 permScore_3 = max(permScore_3),
                 start=nth(start,1), seqnames=nth(seqnames,1),
                 MAF = max(MAF), probeid=nth(probeid,1), mindist=nth(mindist,1))
    GRanges(seqnames=ans$seqnames, IRanges(ans$start, width=1),  # to allow standard 'store' ops
             chisq=ans$chisq, permScore_1=ans$permScore_1,
             permScore_2=ans$permScore_2,
             permScore_3=ans$permScore_3, snp=ans$snp, MAF=ans$MAF,
             probeid=ans$probeid, mindist=ans$mindist)
}
