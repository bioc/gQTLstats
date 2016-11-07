bindHWE.clo = function(samps, tivcf, genome="hg19") function( job, res, ... ) {
  if (is.null(res)) return(NULL)
  rng = range(res$obs)
  vcfdata = readVcf(tivcf, 
     ScanVcfParam(which=rng, samples=samps, info=NA), genome=genome)
  vcfdata = vcfdata[ names(res$obs), ]
  gt = genotypeToSnpMatrix(vcfdata)[[1]]
  sm = col.summary(gt)
  mcols(res$obs)$z.HWE = sm[, "z.HWE"]
  mcols(res$perm)$z.HWE = sm[, "z.HWE"]
  res
}

bhwedemo = function() {
tf17 = TabixFile(dir("~/Research/VCF/", patt=".*17.*gz$", full=TRUE))
#bindHWE = bindHWE.clo( colnames(geuFPKM), tf17 )
#newreg = makeRegistry("updated", packages=c("gQTLstats", "GenomicRanges"))
#batchMapResults( g17@allRegistries[[1]], newreg, bindHWE )
#submitJobs( newreg )

newreg = loadRegistry("updated-files")

n17 = TransStore(list(newreg))
tt17 = tsByRankAccum(n17, 6, filt=function(x) x[x$MAF>.1 & abs(x$z.HWE)<1 ])
}

