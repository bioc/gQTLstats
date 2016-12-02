distToGene = function(buf, annogr) { if (is.null(buf)) return(buf);
  sapply(1:length(buf), function(i) {
  ac = as.character
  ans = rep(NA, length( buf$elnames[i,] ) )
  ans[ ac(seqnames(buf[i])) != ac(seqnames(annogr[ buf$elnames[i,] ])) ] = Inf  # off chrom
  ans[ overlapsAny(annogr[ buf$elnames[i,] ] , buf[i]) ] = 0
  tss = resize(annogr[ buf$elnames[i,] ] , 1 )
  undone = which(is.na(ans))
  if (length(undone)>0) ans[undone] = start(buf[i]) - start(tss[undone])
  ans
})}

collapse_multiPerm = function( se, fblockList, tf, varrng, nperms, bufsize=10 ) {
#
# se = SummarizedExperiment
# fblock = list of feature selections (e.g., chopped up gene list)
# tf = TabixFile to vcf
# varrng = GRanges defining variants under analysis
#
# start iteration over blocks of features
#
# call AllAssoc, computes all pairwise tests
#
 tt0 = try(AllAssoc(se[fblockList[[1]],], tf, variantRange=varrng, nperm=nperms))
 if (is.null(tt0)) {  # AllAssoc will return NULL if no SNP in range or _all_ are Het
   message("variantRange does not yield testable SNP, returning NULL")
   return(NULL)
   }
#
# continue iteration, filtering strongest associations seen so far into buffer 
#
 for (j in 2:length(fblockList)) { # test remaining chunks of transcriptome, retaining best so far
   tt = try(AllAssoc(se[fblockList[[j]],], tf, variantRange=varrng, nperm=nperms))
   tt0 = gQTLstats:::collapseToBuf( tt0, tt, bufsize=bufsize )
   permOut = vector("list", nperms)
   for (k in 1:nperms) {
     if (j == 2) permOut[[k]] = tt0 # must get started, then continue permuted buffer
     permOut[[k]] = gQTLstats:::collapseToBuf( permOut[[k]], tt , bufsize=bufsize, frag=paste0("_permScore_", k))
     gc()
     }
   }
  list(obs=tt0, perms=permOut)
}
