


cisCount = function( summex, vcf.tf, rhs=~1, cisradius=50000, 
    genome="hg19", assayind=1, lbmaf=1e-6, lbgtf = 1e-6, dropUnivHet=TRUE,
    infoFields = c("LDAF", "SVTYPE"),
    simpleSNV=TRUE ) {
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
 if (length(nonSNV)>0) vdata = vdata[-nonSNV,]
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
   gtdata = genotypeToSnpMatrix(vdata)
   }
 varrd = rowRanges(vdata)  # would like to use this as the backbone of test result report
 length(varrd)
}
