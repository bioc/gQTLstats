

library(geuvPack)
data(geuFPKM)
seqlevelsStyle(geuFPKM) = "NCBI"
library(GenomeInfoDb)
ok = which(seqnames(geuFPKM) %in% c(1:22, "X", "Y"))
geuFPKM = geuFPKM[ok,]
library(gQTLBase)
featlist = balancedFeatList( geuFPKM, max=10 )
lens = sapply(featlist,length)
featlist = featlist[ which(lens>0) ]

library(BatchJobs)
reg7 = makeRegistry("reg7_t250k_c500k",  # tile/cis
  packages=c("GenomicRanges", "gQTLstats", "geuvPack",
             "Rsamtools", "VariantAnnotation"), seed=1234)
myf = function(i) {
   if (!exists("geuFPKM")) data(geuFPKM)
   seqlevelsStyle(geuFPKM) = "NCBI"
   curse = geuFPKM[i,]
   load("gsvs.rda")
   svmat = gsvs$sv
   colnames(svmat) = paste0("SV", 1:ncol(svmat))
   colData(curse) = cbind(colData(curse), DataFrame(svmat))
   fmla = as.formula(paste("~", paste0(colnames(svmat), collapse="+")))
   curse = regressOut(curse, fmla)
   pn = gtpath( paste0("chr", as.character(seqnames(curse)[1])) )
   if (as.character(seqnames(curse)[1]) %in% c("4", "9")) {
            tf = TabixFile(pn, index=paste0(basename(pn), ".tbi"))
            }
   else tf = TabixFile(pn)
   cisAssoc( curse, tf, cisradius=500000 )
   }
batchMap(reg7, myf, featlist )
submitJobs(reg7, job.delay = function(n,i) runif(1,1,3))
  

