
library(gQTLstats)
library(BatchJobs)
myr = makeRegistry("myr", packages=c("gQTLstats", "GenomicRanges"))
doit = function(x){
   load("../bigStore.rda");  # this is a ciseStore for GEUVADIS
   storeToMaxAssocBySNP(bigStore, chr=x) }
batchMap(myr, doit, c(1,22:15,2:14))
submitJobs(myr)

