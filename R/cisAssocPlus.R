
sampsInVCF = function(tf) {
#
# probe into VCF file to determine sample names
#  perhaps you don't need the chr?  can determine from header?
#
 ss <- samples(vh <- scanVcfHeader(tf))
 attr(ss, "vh") <- vh
 ss
}

#snvsOnly = function(v) {
#
# confine VCF instance to loci with single nucleotide REF and ALT
#
# 12/25 use new utility
#  v[ width(ref(v)) == 1 & width(GenomicRanges::unlist(alt(v)))==1, ]
#  v[ which(isSNV(v)), ] -- seems unduly limited
#   ev = expand(v)
#   ev[ which(width(ref(ev))==1 & nchar(alt(ev))==1), ] 
#  }


cisAssoc = function( summex, vcf.tf, rhs=~1, nperm=3, cisradius=50000, 
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
   ael = elementLengths(alt(vdata))
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

# use queryVCF
.SnpMatrixCisToSummex = function(summex, vcf.tf,
   radius=50000L, genome="hg19") {
    usn = unique(seqnames(summex))
    if (length(usn) > 1) 
        stop("current implementation insists that length(unique(seqnames(summex)))==1 as VCF assumed chr-specific")
    sampidsInSumm = colnames(summex)
    sampidsInVCF = sampsInVCF(vcf.tf)
    oksamp = intersect(sampidsInSumm, sampidsInVCF)
    stopifnot(length(oksamp) > 0)
    summex = summex[, oksamp]
#    sn = snfilt(as.character(seqnames(summex)))
#    stopifnot(length(ctouse <- unique(sn)) == 1)
    cisr = rowRanges(summex) + radius
#    seqlevels(cisr) = snfilt(seqlevels(cisr))
    vp = ScanVcfParam(fixed = "ALT", info = NA, geno = "GT", 
        samples = oksamp, which = cisr)
    vdata = readVcf(vcf.tf, genome = genome, param = vp)
    rdd = rowRanges(vdata)
#    vdata = snvsOnly(vdata) -- dead function
    genotypeToSnpMatrix(vdata)
}

.eqBox = function( gene, snp, se, tf, radius=1e6, genome="hg19", ...) {
  stopifnot(gene %in% rownames(se))
  LL = SnpMatrixCisToSummex(se[gene,], tf, radius=radius, genome=genome)[[1]]
  stopifnot(snp %in% colnames(LL))
  okids = intersect(colnames(se), rownames(LL))
  stopifnot(length(okids)>0)
  ex = assay(se[ gene, okids])
  gt = as(LL[okids, snp], "character")
  boxplot(split(ex,gt), xlab=snp, ylab=gene, ...)
}

.eqDesc = function (gene, snp, se, tf, radius=1e6, genome="hg19", ...) 
{
    stopifnot(gene %in% rownames(se))
    LL = SnpMatrixCisToSummex(se[gene, ], tf, radius=radius)[[1]]
    stopifnot(snp %in% colnames(LL))
    okids = intersect(colnames(se), rownames(LL))
    stopifnot(length(okids) > 0)
    ex = assay(se[gene, okids])
    gt = as(LL[okids, snp], "character")
    sapply(split(ex, gt), length)
}

prep.cisAssocNB = function( summex, vcf.tf, geneind=1, snpind=1, rhs=~1, nperm=3, cisradius=50000, 
    genome="hg19", assayind="counts", lbmaf=1e-6, dropUnivHet=TRUE,
    infoFields = c("LDAF", "SVTYPE") ) {
#
#
# strategy is to a) select a gene for cis identification of SNP
#                b) select a SNP in this set for analysis
#                c) use all genes for variance estimation
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
 ng = length(rowRanges(summex))
 usn = unique(seqnames(summex))
 if(length(usn)>1) stop("current implementation insists that length(unique(seqnames(summex)))==1 as VCF assumed chr-specific")
 sampidsInSumm = colnames(summex)
 sampidsInVCF = sampsInVCF(vcf.tf)
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
 cisr = rowRanges(summex)[geneind]+cisradius
 seqlevels(cisr) = force(seqlevels(cisr)) # must use VCF-oriented seqlevels
 #
 # first pass at genotype data retrieval
 #
 vp = ScanVcfParam(fixed="ALT", info=infoFields, geno="GT", 
     samples=oksamp, which=cisr)  # which will sort variants into groups defined by probes
 vdata = readVcf(vcf.tf, genome=genome, param=vp) # compressed
 #
 # retain only SNVs with MAF > lbmaf
 # 12/26/2014 -- the exclusion of non-SNVs has become complex
 # we introduce attempt to capture SVTYPE field above.  this may need
 # to be moved up to interface
 svinfo = info(vdata)$SVTYPE
 if (length(svinfo)>0) {
   ok = which(is.na(svinfo))
   vdata = vdata[ok,]
#
# but the example extract has an ALT entry of <DEL> for which SVTYPE is NA
#
   ael = elementLengths(alt(vdata))
   vdata = vdata[ which(ael==1), ]
   stopifnot(length(alt(vdata)) == length(unlist(alt(vdata))))
   todrop = which(!(unlist(alt(vdata)) %in% c("A", "C", "T", "G")))
   if (length(todrop)>0) vdata = vdata[-todrop,]
   tmpalt = try( DNAStringSetList(alt(vdata)) )
   if (inherits( tmpalt, "try-error" )) stop("attempt to reclass ALT fails after SV exclusion")
   alt(vdata) = tmpalt
   }
 nonSNV = which(!isSNV(vdata))
 if (length(nonSNV)>0) vdata = vdata[-nonSNV,]
 gtdata = genotypeToSnpMatrix(vdata)
 uhetinds = NULL
 if (dropUnivHet) {
    message("checking for universal heterozygous loci for exclusion (as dropUnivHet == TRUE) ...")
    gtchar = as(gtdata[[1]],"character")  # could be slow
    uhetinds = which(apply(gtchar,2, function(x) all(x %in% c("A/B", "NA"))))
    if ((nu <- length(uhetinds))>0) 
      warning(paste0("found ", nu, " universally heterozygous loci."))
    message("done checking.")
    }
 csumm = col.summary(gtdata[[1]])
 inmafs = csumm[,"MAF"]
 bad = union(which(inmafs<lbmaf), uhetinds)
 if (length(bad)>0) {
   vdata = vdata[-bad,]
   gtdata = genotypeToSnpMatrix(vdata)
   }
 varrd = rowRanges(vdata)  # would like to use this as the backbone of test result report
# rad = rowRanges(vdata[snpind])+cisradius
# fo = findOverlaps( rowRanges(summex), rad )
# stopifnot(length(subjectHits(fo))>0)
# summex = summex[ subjectHits(fo), ]
 colData(summex)$snp = as(gtdata[[1]], "numeric")[,snpind]
 metadata(colData(summex))$snpgrange = varrd[snpind]
 #list(gtdata=gtdata[[1]][,snpind], varrd=varrd[snpind], summex=summex)
 summex
}

AllAssoc = function( summex, vcf.tf, variantRange, rhs=~1, nperm=3, 
    genome="hg19", assayind=1, lbmaf=1e-6, lbgtf = 1e-6, dropUnivHet=TRUE,
    infoFields = c("LDAF", "SVTYPE")) {
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
 #
 # first pass at genotype data retrieval
 #
 vp = ScanVcfParam(fixed="ALT", info=infoFields, geno="GT", 
     samples=oksamp, which=variantRange)
 vp2 = ScanVcfParam(samples=oksamp[1], which=variantRange)  #  lightweight 
 vdata = readVcf(vcf.tf, genome=genome, param=vp) # compressed
 liteGT = readGT(vcf.tf, param=vp2)  # for verification of readVcf
 #
 # retain only SNVs with MAF > lbmaf
 # 12/26/2014 -- the exclusion of non-SNVs has become complex
 # we introduce attempt to capture SVTYPE field above.  this may need
 # to be moved up to interface
 nonSNV = which(!isSNV(vdata))
 if (length(nonSNV)>0) {
    vdata = vdata[-nonSNV,]
    liteGT = liteGT[-nonSNV,,drop=FALSE]
    }
 if (length(vdata) == 0) {
    message("no SNP identified in range:")
    show(variantRange)
    message("returning NULL")
    return(NULL)
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
 probes2test = rownames(summex)
 numdata = assays(summex)[[assayind]][probes2test,,drop=FALSE]
 #
 # force the formula to have form ex~[rhs]
 #
 infmla = as.formula(paste("ex", paste(as.character(rhs), collapse="")))
 tsts = vector("list", length(probes2test) )
 summs = col.summary(gtdata$genotypes)
 mafs = summs[,"MAF"]
 names(mafs) = rownames(summs)
 #
 # loop over probe list map to collect tests
 #
 if (!exists(".Random.seed")) .xyzzy = runif(1)
 iniSeed = .Random.seed
 suppressWarnings({
 for (i in 1:length(probes2test)) {
   ex = numdata[ probes2test[i], ]
   tsts[[i]] = snp.rhs.tests( formula=infmla, 
           snp.data=gtdata$genotypes[, names(varrd) ], family="gaussian",
           data=data.frame(ex=ex, as(colData(summex), "data.frame") ), uncertain=TRUE )
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
     perms[[j]][[i]] = snp.rhs.tests( formula=infmla, 
             snp.data=gtdata$genotypes[, names(varrd) ], family="gaussian",
             data=data.frame(ex=sample(ex), as(colData(summex), "data.frame") ), uncertain=TRUE )
     }
    })
   }
 names(tsts) = probes2test
 #
 # bind test results to varrd GRanges instance -- note that ALT may be either
 # CharacterList or DNAStringSetList and this may need attention on collection
 #
 #chisqs = unlist(lapply(tsts, chi.squared))
## bind on probe-specific assoc scores for all SNP

 allc = sapply(tsts, chi.squared)  # matrix with SNPs as rows
 savec = colnames(allc)
 colnames(allc)=paste0(savec, "_obs")
 mcols(varrd) = cbind(mcols(varrd), DataFrame(allc))

 pnames = paste0("permScore_", 1:nperm)
 for (i in 1:nperm) {
   allps = sapply( perms[[i]], chi.squared )
   colnames(allps) = paste0(savec, "_permScore_", i)
   mcols(varrd) = cbind(mcols(varrd), DataFrame(allps))
   }
 varrd$snp = names(varrd)
 varrd$MAF = as.numeric(mafs[varrd$snp])
 varrd$probeid = as.character(varrd$paramRangeID)
 metadata(varrd)$sessInfo = sessionInfo()
 metadata(varrd)$init.Random.seed = iniSeed
 metadata(varrd)$dimSummex = dim(summex)
 metadata(varrd)$rowRangesSummex = rowRanges(summex) # should be small
 metadata(varrd)$vcf.tf = vcf.tf # should be small
 metadata(varrd)$vcfHeader = vh # should be small
 metadata(varrd)$dimliteGT = dim(liteGT)
 metadata(varrd)$theCall = thecall
 rm(liteGT)
# names(varrd) = NULL
 
# dropping distance calcs as they can be done post hoc with metadata
# snpl = start(varrd)
# gstart = abs(start(summex[varrd$probeid,]))
# gend = abs(end(summex[varrd$probeid,]))
#
# dists = pmin(abs(snpl-gstart), abs(snpl-gend))
# dists[ which((snpl >= gstart) & (snpl <= gend))] = 0
# varrd$mindist = dists

 varrd
}
