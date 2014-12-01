
manhW878 = function( store, probeid="ENSG00000183814.10", 
  sym="LIN9", fdrsupp, slstyle="NCBI", xlab.in=sym,
  ylab.in="-log10 FDR", applyFDRfilter=TRUE ) {
 sel1 = extractByProbes(store, probeid)
 if (applyFDRfilter) sel1 = fdrsupp@filterUsed(sel1)
 sel1$ml10fdr = -log10( getFDRfunc(fdrsupp)(sel1$chisq)+1e-6 ) # evaluate FDR
 df = as(sel1, "data.frame")
 if (!exists("hmm878")) data(hmm878)  # in gQTLstats
 seqlevelsStyle(hmm878) = slstyle
 hkeep = hmm878[ which(hmm878 %over% sel1) ]
 nn = nearest(sel1, hkeep)
 lab = hkeep$name[nn]
 df$lab = lab
 genem = as(gmod2(sym, collector=transcriptsBy), "data.frame")
 genem$lab="transcripts"
 genem$ml10fdr = -.25 # FAKE
 okn = intersect(names(df), names(genem))
 df2 = rbind( df[,okn], genem[,okn] )
 ggplot( data=df2, aes(y=ml10fdr, x=start, colour=lab)) + 
     geom_point() + xlab(xlab.in) + ylab(ylab.in)
}


