eqBox3 = function (gene, se, tf, snpgr, geneAnno, genome = "hg19", forceRs = TRUE, 
    ...) 
{
#
# geneAnno is a named vector with rownames of se (and 'gene') as names,
# and symbols as values
#
    require(HardyWeinberg)
    ans = gQTLstats:::prepEqData(gene, se, tf, snpgr, genome, forceRs = forceRs)
    ggg = table(ans$gt)
    names(ggg) = gsub("\\/", "", names(ggg))
    #hwp = HWExact(ggg, verbose=FALSE)$pval
    hwp = HWChisq(ggg, verbose=FALSE)$pval
    bb = boxplot(split(ans$ex, ans$gt), plot = FALSE)
    ind = list(...)
    if ("xlab" %in% names(ind)) 
        thexlab = ind$xlab
    else thexlab = ans$coln
    if ("ylab" %in% names(ind)) 
        theylab = ind$ylab
    else theylab = gene
    thexlab = paste(thexlab, "(HWE p =", round(hwp,3),")")
    beeswarm(split(ans$ex, ans$gt), xlab = thexlab, ylab = theylab, 
        main = paste(geneAnno[gene], as.character(seqnames(se[gene, 
            ])), sep = ": "), ...)
    bxp(bb, add = TRUE, boxfill = "transparent")
}
