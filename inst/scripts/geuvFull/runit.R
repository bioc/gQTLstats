library(geuvPack)
data(geuFPKM)
seqlevelsStyle(geuFPKM) == "NCBI"
load("tg_250k.rda")
reqs = 1:length(tg_250k)

library(BatchJobs)
reg5 = makeRegistry("reg5_t250k_c500k",  # tile/cis
  packages=c("GenomicRanges", "gQTLstats", "geuvPack",
             "Rsamtools", "VariantAnnotation"), seed=1234)
myf = function(i) {
   load("tg_250k.rda")
   curti = tg_250k[i]
   if (!exists("geuFPKM")) data(geuFPKM)
   seqlevelsStyle(geuFPKM) = "NCBI"
   curse = subsetByOverlaps(geuFPKM, curti)
   if (nrow(curse)==0) return(NA)
   load("gsvs.rda")
   svmat = gsvs$sv
   colnames(svmat) = paste0("SV", 1:ncol(svmat))
   colData(curse) = cbind(colData(curse), DataFrame(svmat))
   fmla = as.formula(paste("~", paste0(colnames(svmat), collapse="+")))
   curse = regressOut(curse, fmla)
   pn = gtpath( as.numeric(seqnames(curti)) )
   if (as.character(seqnames(curti)) %in% c("4", "9")) {
            tf = TabixFile(pn, index=paste0(basename(pn), ".tbi"))
            }
   else tf = TabixFile(pn)
   cisAssoc( curse, tf, cisradius=500000 )
   }
batchMap(reg5, myf, reqs)
submitJobs(reg5)
  

