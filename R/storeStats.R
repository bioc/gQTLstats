

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
 fdr = oy/(nperm*ncalls)
 ans = data.frame(assoc=xq, fdr=pmin(1, fdr), ncalls=ncalls, avg.false=oy/nperm)
 new("FDRsupp", tab=ans, theCall=theCall, sessinfo=sessionInfo())
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
 
maxByProbe = function (job, res, filter)
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
       fun=function(job, res, filter) maxByProbe(job, res, filter), 
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
    storeToFDR( newstore, xprobs = xprobs, xfield="chisq", getter=getter,
        nperm=nperm )  # data already filtered
}

enumerateByFDR = function (store, fdrsupp, threshold = 0.05,
    filter=force) 
{
    if (!is.function(getFDRfunc(fdrsupp))) {
       message("error: fdrsupp has no FDR interpolating function.")
       message("please add one using setFDRfunc().")
       stop("cannot find FDRfunc element in fdrsupp.")
       }
    selector = function(x) {
        x$estFDR = getFDRfunc(fdrsupp)(x$chisq)
        x[which(x$estFDR <= threshold)]
    }
    ans = storeApply(store, function(x) filter(selector(x)))
    ans = unlist(GRangesList(unlist(ans,recursive=FALSE)))
    metadata(ans)$enumCall = match.call()
    metadata(ans)$enumSess = sessionInfo()
    metadata(ans)$fdrCall = fdrsupp@theCall
    ans
}

