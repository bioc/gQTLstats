setOldClass("function")
setOldClass("gam") # may be too restrictive to assume GAM in use
.nullgam = ""
class(.nullgam) = "gam"
setClassUnion("gamORNULL", c("gam", "NULL"))
setClassUnion("functionORNULL", c("function", "NULL"))

setClass("FDRsupp", representation(tab="data.frame",
   FDRfunc="functionORNULL", FDRmodel="gamORNULL", theCall="call",
   sessinfo="ANY", filterUsed="function"))

setMethod("initialize", "FDRsupp", function(.Object, 
        tab=data.frame(), FDRfunc=NULL, FDRmodel=NULL, 
        sessinfo=sessionInfo(), theCall=call(" "),
        filterUsed=force, ...) {
 .Object@tab = tab
 .Object@FDRfunc = FDRfunc
 .Object@FDRmodel = FDRmodel
 .Object@sessinfo = sessinfo
 .Object@theCall = theCall
 .Object@filterUsed = filterUsed
 .Object
})

setGeneric("getTab", function(x)standardGeneric("getTab"))
setMethod("getTab", "FDRsupp", function(x) x@tab)
setGeneric("getFDRfunc", function(x)standardGeneric("getFDRfunc"))
setMethod("getFDRfunc", "FDRsupp", function(x) x@FDRfunc)
setGeneric("getFDRmodel", function(x)standardGeneric("getFDRmodel"))
setMethod("getFDRmodel", "FDRsupp", function(x) x@FDRmodel)

directPlot = function(FDRsupp) {
 stopifnot(is.function(FDRsupp@FDRfunc))
 tab = getTab(FDRsupp)
 plot( tab$fdr, getFDRfunc(FDRsupp)(tab$assoc),
   xlab="empirical FDR", ylab="modeled FDR")
 abline(0,1)
}

txsPlot.old = function(FDRsupp) {
# plot in transformed space
 stopifnot(is.function(FDRsupp@FDRfunc))
 tab = getTab(FDRsupp)
 plot(qlogis(fdr)~log(assoc), data=getTab(FDRsupp))
 lines(seq(-6,5,.1), predict(FDRsupp@FDRmodel, newdata=
   list(assoc=exp(seq(-6,5,.1)))))
}

txsPlot = function (FDRsupp, xmax=50) 
{
    stopifnot(is.function(FDRsupp@FDRfunc))
    tab = getTab(FDRsupp)
    plot(-log10(fdr+1e-6) ~ assoc, data = getTab(FDRsupp), xlim=c(0,xmax))
    lines(1:xmax, predict(FDRsupp@FDRmodel, newdata = list(assoc =1:xmax)))
        
}



setMethod("show", "FDRsupp", function(object) {
 cat("FDRsupp instance with ", nrow(object@tab), " rows.\n", sep="")
 print(head(object@tab,3))
 cat("...\n")
 print(tail(object@tab,3))
 if (is.function(object@FDRfunc)) cat("An interpolating function is available.\n")
 else cat("No interpolating function is available; use 'setFDRfunc'.\n")
}
)

#setFDRfunc = function(FDRsupp, fudge=1e-6, zthresh=30, ...) {
##
## may want to generalize from smoother choice lo()
##
# qlfdrmod = gam(qlogis(fdr+fudge)~lo(assoc,...), data=getTab(FDRsupp))
# FDRsupp = addFDRmodel(FDRsupp, qlfdrmod)
# addFDRfunc(FDRsupp, function(assoc)
#   ifelse(assoc < zthresh, plogis(predict(getFDRmodel(FDRsupp),
#     newdata=list(assoc=assoc, data=getTab(FDRsupp))))-fudge, 0))
#}

setFDRfunc = function (FDRsupp, fudge = 1e-06, zthresh = 30, maxch=30, ...)
{
  gt = getTab(FDRsupp)
  fdrmod = gam(-log10(fdr+fudge)~s(assoc,bs="tp"), data=gt,
    subset=assoc<(1.1*maxch))
  FDRsupp = addFDRmodel(FDRsupp, fdrmod)
  FDRsupp@FDRfunc = function(assoc)
      10^-predict(fdrmod, newdata=list(assoc=assoc)) - fudge
  FDRsupp
}

