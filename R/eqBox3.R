eqBox3 = function (gene, se, tf, snpgr, geneAnno, genome = "hg19", forceRs = TRUE, 
    ...) 
{
#
# geneAnno is a named vector with rownames of se (and 'gene') as names,
# and symbols as values
#
    ans = prepEqData(gene, se, tf, snpgr, genome, forceRs = forceRs)
    ggg = table(ans$gt)
    names(ggg) = gsub("\\/", "", names(ggg))
    #hwp = HWExact(ggg, verbose=FALSE)$pval
    hwp = HardyWeinberg::HWChisq(ggg, verbose=FALSE)$pval
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
eqBox4 = function (gene, se, tf, snpgr, geneAnno, genome = "hg19", forceRs = TRUE, 
    ...) 
{
#
# geneAnno is a named vector with rownames of se (and 'gene') as names,
# and symbols as values
#
    ans = prepEqData(gene, se, tf, snpgr, genome, forceRs = forceRs)
    ggg = table(ans$gt)
    names(ggg) = gsub("\\/", "", names(ggg))
    #hwp = HWExact(ggg, verbose=FALSE)$pval
    hwp = HardyWeinberg::HWChisq(ggg, verbose=FALSE)$pval
    bb = boxplot(split(ans$ex, ans$gt), plot = FALSE)
    ind = list(...)
    if ("xlab" %in% names(ind)) 
        thexlab = ind$xlab
    else thexlab = ans$coln
    if ("ylab" %in% names(ind)) 
        theylab = ind$ylab
    else theylab = gene
    thexlab = paste(thexlab, "(HWE p =", round(hwp,3),")")
    dd = data.frame(ex=as.numeric(ans$ex), gt=factor(ans$gt),
          id=colnames(ans$ex))
    maintitle= paste(geneAnno[gene], as.character(seqnames(se[gene, 
            ])), sep = ": ")
    ggp = ggplot(dd, aes(x=gt, y=ex, text=id)) + geom_boxplot() +
           geom_quasirandom() + labs(x=thexlab, y=theylab) +
           ggtitle(maintitle)
    if (!is.null(ggp)) ggplotly(ggp)
#    beeswarm(split(ans$ex, ans$gt), xlab = thexlab, ylab = theylab, 
#        main = paste(geneAnno[gene], as.character(seqnames(se[gene, 
#            ])), sep = ": "), ...)
#    bxp(bb, add = TRUE, boxfill = "transparent")
}
