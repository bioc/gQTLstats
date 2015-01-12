

storeToQuantiles = function( store, field, 
     probs=c(seq(0,.999,.001), 1-(c(1e-4,1e-6,1e-6,1e-7))), 
     ids=NULL, ..., checkField=FALSE, filter=force ) {
#
# uses ff for interim storage of data for summarization
#
  src = storeToFf( store, field, ids=ids, checkField=checkField, filter=filter )
  quantile(src, probs )
}

storeToHist = function( store, getter =
   function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[,c("permScore_1",
            "permScore_2", "permScore_3")])), breaks, ids=NULL,
            filter=force ) {
   if (missing(breaks)) stop("breaks must be supplied and must cover range of data")
   if (is.null(ids)) ids=store@validJobs
   tmp = bplapply(ids, function(x) {
      getter(loadAndFilterResult(reg=store@reg, id=x, filter=filter))
      })
   tmp = bplapply(tmp, function(x) try(hist(x , breaks=breaks, plot=FALSE )))
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
    1e-06, 1e-06, 1e-07))), xfield = "chisq",
    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
       c("permScore_1", "permScore_2", "permScore_3")])), nperm=3,
    filter=force) {
 theCall = match.call()
 message("counting tests...")
 ntests = sum(unlist(storeApply( store, function(x) length(filter(x)) )) , na.rm=TRUE)
 message("counting #NA...")
 nna = sum(unlist(storeApply( store, function(x)sum(is.na(filter(x)$chisq)) )) , na.rm=TRUE)
 ntests = ntests - nna
 message("obtaining assoc quantiles...")
 xq = storeToQuantiles(store, field=xfield, probs=xprobs, filter=filter) # nxq
 message("computing perm_assoc histogram....")
 yh.orig = storeToHist(store, getter=getter, breaks=c(0,xq,1e10), filter=filter)$counts # nxq+2
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
 
.maxByProbe = function (job, res, filter)
{
    cur = filter(res)
    sr = split(cur, cur$probeid)
    tmpp = vector("list", length(sr))
    for (i in 1:length(sr)) {
        tmp = sr[[i]][which.max(sr[[i]]$chisq)]
        mcols(tmp)[, "permScore_1"] = max(sr[[i]]$permScore_1)
        mcols(tmp)[, "permScore_2"] = max(sr[[i]]$permScore_2)
        mcols(tmp)[, "permScore_3"] = max(sr[[i]]$permScore_3)
        tmpp[[i]] = tmp
    }
    unlist(GRangesList(tmpp))
}

maxByProbe = function (job, res, filter)
{
    as(filter(res), "data.frame") %>% dplyr::select(snp, probeid,
           chisq, permScore_1, permScore_2, permScore_3,
           seqnames, start) %>% group_by(probeid) %>%
      summarize(chisq=max(chisq), permScore_1=max(permScore_1),
            permScore_2=max(permScore_2), permScore_3=max(permScore_3))
}

maxBySNP = function (job, res, filter)
{
#
# you probably need to use ffdf for this...
#
    as(filter(res), "data.frame") %>% dplyr::select(snp, probeid, MAF,
           chisq, permScore_1, permScore_2, permScore_3,
           seqnames, start) %>% group_by(snp) %>%
      summarize(chisq=max(chisq), permScore_1=max(permScore_1),
            permScore_2=max(permScore_2), permScore_3=max(permScore_3),
            seqnames=nth(seqnames,1), start=nth(start,1), MAF=nth(MAF,1))  # last 3 const within SNP
}

storeToFDRByProbe = function( store, xprobs = c(seq(0, 0.999, 0.001), 1 - (c(1e-04,
    1e-05, 1e-06, 1e-07))), xfield = "chisq",
    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
       c("permScore_1", "permScore_2", "permScore_3")])), nperm=3, filter=force) {
    curbb = options()$BBmisc.ProgressBar.style
    options(BBmisc.ProgressBar.style="off")
    tf = tempfile()
    suppressMessages({
    newreg = makeRegistry("byp", file.dir=tf)
    batchMapResults( store@reg, newreg, ids=store@validJobs,
       fun=function(job, res, filter) .maxByProbe(job, res, filter), 
       more.args=list(filter=filter) )
    submitJobs( newreg )#, store@validJobs )
    waitForJobs( newreg )
    })
    if (length(findNotDone(newreg))>0) {
      warning("could not finish maxByProbe; returning recipient registry")
      return(newreg)
      }
    newstore = ciseStore( newreg, FALSE, FALSE )
    on.exit({unlink(tf, recursive=TRUE);
             options(BBmisc.ProgressBar.style=curbb)
             })
    ans = storeToFDR( newstore, xprobs = xprobs, xfield="chisq", getter=getter,
        nperm=nperm )  # data already filtered
    ans@filterUsed = filter
    ans@theCall = match.call()
    ans
}

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
    ul = unlist(ans, recursive=FALSE)
    zl = sapply(ul,length)
    if (any(zl==0)) ul=ul[-which(zl==0)]
    ans = unlist(GRangesList(ul))
    metadata(ans)$enumCall = match.call()
    metadata(ans)$enumSess = sessionInfo()
    metadata(ans)$fdrCall = fdrsupp@theCall
    ans
}

dfrToFDR = function(dfr, xprobs = c(seq(0, 0.999, 0.001), 1 - (c(1e-04,
    1e-06, 1e-06, 1e-07))), xfield = "chisq",
    getter = function(x) as.numeric(S4Vectors::as.matrix(mcols(x)[, 
       c("permScore_1", "permScore_2", "permScore_3")])), nperm=3,
    filter=force) {
#
# when in-memory is OK
#
# data.frame method
 theCall = match.call()
 dfr = filter(dfr)
 stopifnot("chisq" %in% names(dfr))
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

storeToFDRByProbe = function( store, xprobs = seq(0, 0.99, 0.01),
       xfield = "chisq", nperm=3, filter=force) {
    maxbp = unlist(storeApply(store, function(x) {
       maxByProbe(1, filter(x), force) }), recursive=FALSE)
    maxbp = rbind_all(maxbp)  # note filter was already applied
    dfrToFDR(maxbp, xprobs = xprobs, xfield = xfield, filter=force)
}

storeToMaxAssocBySNP = function( store, chr=1, xprobs = seq(0, 0.999, 0.001),
    xfield = "chisq", nperm=3, filter=force) {
    stopifnot(is.atomic(chr) & length(chr)==1)
    if (length(store@rangeMap)==0) stop("please use a store with valid rangeMap component.")
    rmap = store@rangeMap
    ac = as.character
    jobsToDo = as.integer(rmap[ which(ac(seqnames(rmap))==ac(chr)) ]$jobid)
    stopifnot(length(jobsToDo)>0)
    maxbsByJob = unlist(storeApply(store, function(x) {
       maxBySNP(1, filter(x), force) }, ids=jobsToDo ), recursive=FALSE)
 #   maxbp = rbind_all(maxbsByJob)  # note filter was already applied
 #  here you need to group by once more and max by snp
 #   maxbp  #dfrToFDR(maxbp, xprobs = xprobs, xfield = xfield, filter=force)
 #   aggr = maxbsByJob[[1]]
 #   for (i in 2:length(maxbsByJob)) {
 #       cat(i)
 #       aggr = rbind( aggr, maxbsByJob[[i]] )
 #       }
    message("binding all jobs")
    aggr = do.call(rbind, maxbsByJob)
    message("group by and max")
    # after collecting, use
    #dfrToFDR( aggr, xprobs = xprobs, xfield=xfield )
    ans <- aggr %>% group_by(snp) %>%
               summarize( chisq = max(chisq),
                 permScore_1 = max(permScore_1),
                 permScore_2 = max(permScore_2),
                 permScore_3 = max(permScore_3),
                 start=nth(start,1), seqnames=nth(seqnames,1),
                 MAF = max(MAF))
    GRanges(seqnames=ans$seqnames, IRanges(ans$start, width=1),  # to allow standard 'store' ops
             chisq=ans$chisq, permScore_1=ans$permScore_1,
             permScore_2=ans$permScore_2,
             permScore_3=ans$permScore_3, snp=ans$snp, MAF=ans$MAF)
}
