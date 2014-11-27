
queryVCF = function (gr, vcf.tf, samps, genome = "hg19", getSM=TRUE) 
{
#
# support anonymous extraction of VCF/SnpMatrix by range
#
    vp = ScanVcfParam(fixed = "ALT", info = NA, geno = "GT", 
        which = gr)
    if (!missing(samps)) {
      present = samples(scanVcfHeader(vcf.tf))
      if((dl <- length(setdiff(samps, present))) > 0) {
        message(paste0("NOTE: there were ", dl, " samples not found (of ",
            length(samps)," requested)."))
        }
      vcfSamples(vp) = (int <- intersect(samps, present))
      stopifnot(length(int)>0)
      }
    readout = readVcf(vcf.tf, genome = genome, param = vp)
    sm = NULL
    if (getSM) sm = genotypeToSnpMatrix(readout)
    list(readout=readout, sm=sm)
}

prepEqData = function (gene, se, tf, snpgr, genome = "hg19") {
    stopifnot(gene %in% rownames(se))
    esamps = colnames(se)
    gtstuff = queryVCF(gr=snpgr, vcf.tf=tf, samps=esamps,
               genome=genome, getSM=TRUE) # [[1]]:readout, [[2]]:(SnpMatrix,map)
    smat = gtstuff[[2]][[1]]
    stopifnot(prod(dim(smat))>0)
    okids = intersect(esamps, rownames(smat))
    ex = assay(se[gene, okids])
    gt = as.character(as(smat[okids,], "character"))
    list(ex=ex, gt=gt, coln=colnames(smat))
}

eqBox2 = function (gene, se, tf, snpgr, genome = "hg19", 
    ...) {
    ans = prepEqData( gene, se, tf, snpgr, genome )
    boxplot(split(ans$ex, ans$gt), xlab = ans$coln, ylab = gene, ...)
}

eqDesc2 = function (gene, se, tf, snpgr, genome = "hg19") {
    ans = prepEqData( gene, se, tf, snpgr, genome )
    sapply(split(ans$ex, ans$gt),length)
}

