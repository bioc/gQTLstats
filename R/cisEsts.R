
cisEsts = function( summex, vcf.tf, rhs=~1, nperm=3, cisradius=50000, 
    genome="hg19", assayind=1, lbmaf=1e-6, lbgtf = 1e-6, dropUnivHet=TRUE,
    infoFields = c("LDAF", "SVTYPE"),
    simpleSNV=TRUE) {
#
# LDAF is concession to bug in readVcf where specifying only SVTYPE
# leads to error on 1KG VCF data
#
 #
 # take all features from RangedSummarizedExperiment
 # harmonize samples between summex and vcf.tf (TabixFile)
 # obtain genotypes of variants cis to features in summex -- using only SNVs!!
 # compute associations between stx(features) and vtx(genotypes)
 # store in a GRanges ordered by variant address with relevant metadata
 #
 # obtain sample IDs and harmonize genotypes and molec phenotype assay data
 #
 thecall = match.call()
 usn = unique(seqnames(summex))
 if(length(usn)>1) stop("current implementation insists that length(unique(seqnames(summex)))==1 as VCF assumed chr-specific")
 sampidsInSumm = colnames(summex)
 sampidsInVCF = sampsInVCF(vcf.tf)
 vh = attr(sampidsInVCF, "vh")
 oksamp = intersect(sampidsInSumm, sampidsInVCF)
 stopifnot(length(oksamp)>0)
 summex = summex[, oksamp]
 #
 # harmonize annotation for seqnames -- could use style methods here
 #
 sn = force(as.character(seqnames(summex)))
 stopifnot(length(ctouse <- unique(sn))==1)
 #
 # generate cis search space for assay probes
 #
 cisr = trim(rowRanges(summex)+cisradius)
 requestSize = length(cisr)
 seqlevels(cisr) = force(seqlevels(cisr)) # must use VCF-oriented seqlevels
 #
 # first pass at genotype data retrieval
 #
 vp = ScanVcfParam(fixed="ALT", info=infoFields, geno="GT", 
     samples=oksamp, which=cisr)  # which will sort variants into groups defined by probes
 vp2 = ScanVcfParam(samples=oksamp[1], which=cisr)  #  lightweight 
 vdata = readVcf(vcf.tf, genome=genome, param=vp) # compressed
 liteGT = readGT(vcf.tf, param=vp2)  # for verification of readVcf
 #
 # retain only SNVs with MAF > lbmaf
 # 12/26/2014 -- the exclusion of non-SNVs has become complex
 # we introduce attempt to capture SVTYPE field above.  this may need
 # to be moved up to interface
#
# if we can avoid complex SNV handling
#
 if (!simpleSNV) {
 svinfo = info(vdata)$SVTYPE
 if (length(svinfo)>0) {
   ok = which(is.na(svinfo))
   vdata = vdata[ok,]
#
# but the example extract has an ALT entry of <DEL> for which SVTYPE is NA
#
   ael = elementNROWS(alt(vdata))
   vdata = vdata[ which(ael==1), ]
   stopifnot(length(alt(vdata)) == length(unlist(alt(vdata))))
   todrop = which(!(unlist(alt(vdata)) %in% c("A", "C", "T", "G")))
   if (length(todrop)>0) vdata = vdata[-todrop,]
   tmpalt = try( DNAStringSetList(alt(vdata)) )
   if (inherits( tmpalt, "try-error" )) stop("attempt to reclass ALT fails after SV exclusion")
   alt(vdata) = tmpalt
   }
 }
 nonSNV = which(!isSNV(vdata))
 if (length(nonSNV)>0) {
    vdata = vdata[-nonSNV,]
    liteGT = liteGT[-nonSNV,,drop=FALSE]
    }
 gtdata = genotypeToSnpMatrix(vdata)
 uhetinds = NULL
 if (dropUnivHet) {
    message("checking for universal heterozygous loci for exclusion (as dropUnivHet == TRUE) ...")
    gtchar = as(gtdata[[1]],"character")  # could be slow
    uhetinds = which(apply(gtchar,2, function(x) all(x %in% c("A/B", "NA"))))
    if ((nu <- length(uhetinds))>0) 
      warning(paste0("found ", nu, " universally heterozygous loci."))
    if (nu == ncol(gtdata[[1]])) {
        warning("all loci universally heterozygous, returning NULL")
        return(NULL)
        }
    message("done checking.")
    }
 csumm = col.summary(gtdata[[1]])
 inmafs = csumm[,"MAF"]
 ingtmat = csumm[,c("P.AA", "P.AB", "P.BB")]
 lowgt = apply(ingtmat,1,min,na.rm=TRUE)
 bad = union(which(inmafs<lbmaf | lowgt<lbgtf), uhetinds)
 if (length(bad)>0) {
   vdata = vdata[-bad,]
   liteGT = liteGT[-bad,,drop=FALSE]
   gtdata = genotypeToSnpMatrix(vdata)
   }
 varrd = rowRanges(vdata)  # would like to use this as the backbone of test result report
 #
 # use a list mapping probes to SNVs in cis to organize the testing
 #
 uo = unique(varrd$paramRangeID)
 nRequestsSatisfied = length(uo)
 if (requestSize != nRequestsSatisfied) warning("number unique values of paramRangeID returned differs from number requested")
 snpbyprobe = split(names(varrd), varrd$paramRangeID)[uo] # will reorder without uo
 probes2test = names(snpbyprobe)
 numdata = assays(summex)[[assayind]][probes2test,,drop=FALSE]
 #
 # force the formula to have form ex~[rhs]
 #
 infmla = as.formula(paste("ex", paste(as.character(rhs), collapse="")))
 tsts = vector("list", length(probes2test) )
 summs = col.summary(gtdata$genotypes)
 mafs = summs[,"MAF"]
 names(mafs) = rownames(summs)

 doEsts = TRUE
 if (doEsts) ests = vector("list", length(probes2test) )

 #
 # loop over cis map to collect tests
 #
 if (!exists(".Random.seed")) .xyzzy = runif(1)
 iniSeed = .Random.seed
 suppressWarnings({
 for (i in 1:length(probes2test)) {
   ex = numdata[ probes2test[i], ]
   tsts[[i]] = snp.rhs.tests( formula=infmla, 
           snp.data=gtdata$genotypes[, snpbyprobe[[ probes2test[i] ]] ], family="gaussian",
           data=data.frame(ex=ex, as(colData(summex), "data.frame") ), uncertain=TRUE )
    if (doEsts) {
           ests[[i]] = snp.rhs.estimates( formula=infmla, 
            snp.data=gtdata$genotypes[, snpbyprobe[[ probes2test[i] ]] ], family="gaussian",
            data=data.frame(ex=ex, as(colData(summex), "data.frame") ), uncertain=TRUE )
                 }
   }
  })
 #
 # test under permutation for plug-in FDR
 #
 perms = vector("list", nperm)
 for (j in 1:nperm) {
   perms[[j]] = vector("list", length(probes2test))
   suppressWarnings({
   for (i in 1:length(probes2test)) {
     ex = as.numeric(numdata[ probes2test[i], ])
     if (any(is.na(ex))) {
         print("perm")
         print(j)
         print(probes2test[i])
         stop("NA sneaks in")
         }
     perms[[j]][[i]] = snp.rhs.tests( formula=infmla, 
             snp.data=gtdata$genotypes[, snpbyprobe[[ probes2test[i] ]] ], family="gaussian",
             data=data.frame(ex=sample(ex), as(colData(summex), "data.frame") ), uncertain=TRUE )
     }
    })
   }
 names(tsts) = probes2test
 #
 # bind test results to varrd GRanges instance -- note that ALT may be either
 # CharacterList or DNAStringSetList and this may need attention on collection
 #
 chisqs = unlist(lapply(tsts, chi.squared))
 varrd$chisq = chisqs
 pnames = paste0("permScore_", 1:nperm)
 for (i in 1:nperm) {
   mcols(varrd)[ ,pnames[i] ] = unlist(lapply(perms[[i]], chi.squared))
   }
 varrd$snp = names(varrd)
 varrd$MAF = as.numeric(mafs[varrd$snp])
 varrd$probeid = as.character(varrd$paramRangeID)
  if (doEsts) {
      varrd$beta = 
        unlist(lapply(ests, function(x) sapply(x@.Data, function(x) x$beta)))
      varrd$se.beta = 
        unlist(lapply(ests, function(x) sapply(x@.Data, function(x) sqrt(x$Var.beta))))
      }

 metadata(varrd)$sessInfo = sessionInfo()
 metadata(varrd)$init.Random.seed = iniSeed
 metadata(varrd)$dimSummex = dim(summex)
 metadata(varrd)$rowRangesSummex = rowRanges(summex) # should be small
 metadata(varrd)$vcf.tf = vcf.tf # should be small
 metadata(varrd)$vcfHeader = vh # should be small
 metadata(varrd)$requestSize = requestSize
 metadata(varrd)$nRequestsSatisfied = nRequestsSatisfied
 metadata(varrd)$dimliteGT = dim(liteGT)
 metadata(varrd)$theCall = thecall
 rm(liteGT)
 names(varrd) = NULL
 
 snpl = start(varrd)
 gstart = abs(start(summex[varrd$probeid,]))
 gend = abs(end(summex[varrd$probeid,]))

 dists = pmin(abs(snpl-gstart), abs(snpl-gend))
 dists[ which((snpl >= gstart) & (snpl <= gend))] = 0
 varrd$mindist = dists

 varrd
}

