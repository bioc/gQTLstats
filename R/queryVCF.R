
queryVCF = function (gr, vcf.tf, samps, genome = "hg19", getSM=TRUE, snvOnly=TRUE) 
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
    if (snvOnly) readout = readout[which(isSNV(readout)),]
    sm = NULL
    if (getSM) sm = genotypeToSnpMatrix(readout)
    list(readout=readout, sm=sm)
}

prepEqData = function (gene, se, tf, snpgr, genome = "hg19",
    forceRs=TRUE) {
    stopifnot(gene %in% rownames(se))
    esamps = colnames(se)
    gtstuff = queryVCF(gr=snpgr, vcf.tf=tf, samps=esamps,
               genome=genome, getSM=TRUE) # [[1]]:readout, [[2]]:(SnpMatrix,map)
    smat = gtstuff[[2]][[1]]
    if (ncol(smat)>1 & forceRs)
       smat = smat[, grep("rs", colnames(smat))] # assume DELLY or other ids need to be dropped, and there is only on dbSNP id, but if only one variant we are OK
    stopifnot(prod(dim(smat))>0)
    okids = intersect(esamps, rownames(smat))
    ex = assay(se[gene, okids])
    gt = as.character(as(smat[okids,], "character"))
    list(ex=ex, gt=gt, coln=colnames(smat), sampids=okids)
}

gQTLs = function(filtgr, se, tf, genome="hg19", forceRs=TRUE, chunksize=50) {
    # purpose is to generate a SummarizedExperiment from filtered ciseStore
    # filtgr is GRanges from filtering ciseStore, tf is chromosome-specific VCF so
    # be sure that filtgr is properly restricted
    thecall = match.call()
    stopifnot(all(c("probeid", "snp") %in% names(mcols(filtgr))))
    vh = scanVcfHeader(tf)
    vcfids = vcfSamples(vh)
    okids = intersect(vcfids, colnames(se))
    stopifnot(length(okids)>0)
    ex = assay(se[ filtgr$probeid, okids ])
#
# to improve robustness, in case readVcf does not return everything
# we expect (perhaps two snp names share an addr, we use snp:gene
# as rownames
#
    rownames(ex) = paste(filtgr$snp, ":", filtgr$probeid, sep="")
    chunks = chunk(1:length(filtgr), chunk.size=chunksize)
    snpd = foreach (ch = chunks) %dopar% {
          q1 = queryVCF(gr=filtgr[ch], vcf.tf=tf, samps=okids,
               genome=genome, getSM=TRUE)$sm[[1]] 
          t(as(q1[okids,], "character"))
          }
    snpd = do.call(rbind, snpd) #
#    rownames(snpd) = NULL
    # we may have a mismatch if vcf does not supply all snp ...
    # needs further checking
    if (nrow(snpd) != nrow(ex)) {
      message("some SNP requests seem unfilled ... multiple names for same addr?")
      exkeys = gsub(":.*", "", rownames(ex))
      oksnp = intersect(rownames(snpd), exkeys)
      ex = ex[ which(exkeys %in% oksnp), ] # will drop non-intersecting
      snpd = snpd[ which(rownames(snpd) %in% oksnp), ]  # likewise
      filtgr = filtgr[ which(filtgr$snp %in% oksnp) ] # likewise
      }
    stopifnot( nrow(snpd) == nrow(ex) )
    rownames(ex) = NULL
    se = SummarizedExperiment(List(calls=snpd, exprs=ex), colData=colData(se)[colnames(ex),])
    names(filtgr) = filtgr$snp
    rowRanges(se) = filtgr
    metadata(se)$thecall = thecall
    metadata(se)$tfpath = path(tf)
    se
}
  

eqBox2 = function (gene, se, tf, snpgr, genome = "hg19", forceRs = TRUE,
    ...) {
    ans = prepEqData( gene, se, tf, snpgr, genome, forceRs=forceRs )
    boxplot(split(ans$ex, ans$gt), xlab = ans$coln, ylab = gene, ...)
}

eqDesc2 = function (gene, se, tf, snpgr, genome = "hg19", forceRs=TRUE) {
    ans = prepEqData( gene=gene, se=se, tf=tf, snpgr=snpgr, genome=genome, forceRs=forceRs )
    sapply(split(ans$ex, ans$gt),length)
}

setClass("SnpToGeneQTL", representation(exprs="numeric",
   genotypes="ANY", gene="character", snpid="character",
   snpgr="GRanges", genome="character", boxdata="ANY",
   covars="DataFrame"))
setGeneric("SnpToGeneQTL", function(gene, snpid, snpgr, tf, se, genome, forceRs) 
   standardGeneric("SnpToGeneQTL"))

setMethod("SnpToGeneQTL", c(gene="character", snpid="character", snpgr="GRanges",
     tf="TabixFile", se="SummarizedExperiment"), function(gene, snpid, snpgr, tf, se,
         genome="hg19", forceRs=TRUE) {
    ans = prepEqData( gene, se, tf, snpgr, genome, forceRs=forceRs )
    ab = boxplot(split(ans$ex, ans$gt), plot=FALSE)
    new("SnpToGeneQTL", exprs=as.numeric(ans$ex), genotypes=ans$gt, gene=gene, snpgr=snpgr, genome=genome, boxdata=
       ab, covars=colData(se)[colnames(ans$ex),], snpid=snpid)
})

setMethod("SnpToGeneQTL", c(gene="missing", snpid="missing", snpgr="GRanges",
     tf="TabixFile", se="SummarizedExperiment"), function(gene, snpid, snpgr, tf, se,
         genome="hg19", forceRs=TRUE) {
    snpid = mcols(snpgr)$snp
    gene = mcols(snpgr)$probeid
    ans = prepEqData( gene, se, tf, snpgr, genome, forceRs=forceRs )
    ab = boxplot(split(ans$ex, ans$gt), plot=FALSE)
    new("SnpToGeneQTL", exprs=as.numeric(ans$ex), genotypes=ans$gt, gene=gene, snpgr=snpgr, genome=genome, boxdata=
       ab, covars=colData(se)[colnames(ans$ex),], snpid=snpid)
})




setMethod("show", "SnpToGeneQTL", function(object){
cat("gQTLstats SnpToGeneQTL instance:\n")
cat("SNP ", object@snpid, "; gene: ", object@gene, "\n", sep="")
cat("use boxswarm() for visualization\n")
})
setGeneric("boxswarm", function(x, colvar, colmap, inxlab, inylab, ...)
     standardGeneric("boxswarm"))
setMethod("boxswarm", "SnpToGeneQTL", function(x, colvar, colmap, inxlab, inylab, ...) {
 if (!missing(colvar)) {
     cov2use = x@covars[[colvar]]
     collist = split(cov2use, x@genotypes)
     ucov = unique(cov2use)
     if (missing(colmap)) {
        colmap = 1:length(ucov)
        names(colmap) = ucov
        }
     pwcols = lapply(collist, function(x) colmap[x])
     }
 else pwcols=NULL
 if (missing(inylab)) inylab=g1@gene
 if (missing(inxlab)) inxlab=g1@snpid
 beeswarm(split(x@exprs, x@genotypes), pwcol=pwcols, pch=19, xlab=inxlab, ylab=inylab, ...)
 bd = x@boxdata
 bd$out=NULL
 bd$group=NULL
 bxp(bd, add=TRUE, boxfill="transparent")
})

gQTLswarm = function (se, ind, covar = NULL, inpch=19, xlab, ylab, featTag="probeid", ...) 
{
    gt = assays(se)[[1]][ind, ]
    ex = assays(se)[[2]][ind, ]
    el = split(ex, gt)
    colp = NULL
    if (!is.null(covar)) 
        colp = split(as.numeric(factor(colData(se)[[covar]])), 
            gt)
    if (missing(xlab)) xlab=rowRanges(se)$snp[ind]
    if (missing(ylab)) ylab=mcols(rowRanges(se))[[featTag]][ind]
#    beeswarm(el, pwcol = colp, pch = inpch, xlab=xlab, ylab=ylab, ...)
    bp = boxplot(el, xlab=xlab, ylab=ylab, notch=TRUE, pch=" ", ...) #plot = FALSEbb)
    beeswarm(el, pwcol = colp, pch = inpch, add=TRUE)
    #bxp(bp, add = TRUE, bg = "transparent")
}

