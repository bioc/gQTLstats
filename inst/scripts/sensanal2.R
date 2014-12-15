
# Following moved to parms.R
# MAFS = c(.03, .04, .05, .075, .10, .125, .15)
# dists = c(5000, 7500, 10000, 15000, 20000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000)
# parms = expand.grid(MAFS, dists)
library(BatchJobs)  # for bigStore manip
library(gQTLstats)
load("../bigStore.rda")
library(BiocParallel)

# could use multilevel parallelism here
# because it is a somewhat large, fragile job, BatchJobs
# is a relevant tool for iteration.  but storeToFDRByProbe is
# already using bplapply.  so register 3 cores for it 
# and specify 15 cpu for BatchJobs in .BatchJobs.R

register(MulticoreParam(workers=3))

sens1 = makeRegistry("sens1", file.dir="sens1", 
    packages=c("BiocParallel", "gQTLstats", "dplyr"),
    src.files="parms.R")

sens4One = function(z) {
      library(BiocParallel)
      load("../bigStore.rda")
      register(MulticoreParam(workers=3))  
      ans = storeToFDRByProbe(bigStore, xprobs=seq(.01,.99,.01),
      filter=function(x) x[which(x$MAF >= parms[z,1] & 
          x$mindist <= parms[z,2])])
      ans = setFDRfunc(ans, span=.35)
      list(fdrsupp=ans, parms=parms[z,])
 }

batchMap(sens1, sens4One, 1:nrow(parms))
submitJobs(sens1)

#sens = lapply(1:nrow(parms), function(z) {
#  ans = storeToFDRByProbe(bigStore, xprobs=seq(.01,.99,.01),
#      filter=function(x) x[which(x$MAF >= parms[z,1] & 
#        x$mindist <= parms[z,2])])
#  ans = setFDRfunc(ans, span=.35)
#  ans
#})

