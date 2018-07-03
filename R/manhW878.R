
manhWngr = function( store, probeid="ENSG00000183814.10", 
  sym="LIN9", fdrsupp, namedGR, slstyle="NCBI", xlab.in=sym,
  ylab.in="-log10 FDR", applyFDRfilter=TRUE ) {
#
# manhattan plot for eQTL association FDRs, with coloring
# via a column in a named GRanges -- column 'name' will determine
# coloring scheme
#
 stopifnot( "name" %in% names(mcols(namedGR)) )
 stopifnot( is.character(mcols(namedGR)[,"name"]) )
 sel1 = extractByProbes(store, probeid)
 if (applyFDRfilter) sel1 = fdrsupp@filterUsed(sel1)
 sel1$ml10fdr = -log10( getFDRfunc(fdrsupp)(sel1$chisq)+1e-6 ) # evaluate FDR
 sel1$ml10fdr = as.numeric(sel1$ml10fdr)
 df = as(sel1, "data.frame")
 seqlevelsStyle(namedGR) = slstyle
 hkeep = namedGR[ which(namedGR %over% sel1) ]
 nn = nearest(sel1, hkeep)
 lab = hkeep$name[nn]
 df$lab = lab
 genem = as(gmod2(sym, collector=transcriptsBy), "data.frame")
 genem$lab="transcripts"
 genem$ml10fdr = -.25 # FAKE
 okn = intersect(names(df), names(genem))
 df2 = rbind( df[,okn], genem[,okn] )
 xlab2 = paste0( xlab.in, " (", genem$seqnames[1], ")")
 ggplot( data=df2, aes(y=ml10fdr, x=start, colour=lab)) + 
     geom_point() + xlab(xlab2) + ylab(ylab.in)
}


