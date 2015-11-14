transAssoc = function(variantGR, exSE, vcfgen, bufsize=10, nperm=3,
     exChLen=2*bufsize, ...) {
  stopifnot(genome(variantGR)[1] == genome(exSE)[1])
  stopifnot(seqlevelsStyle(variantGR)[1] == seqlevelsStyle(exSE)[1])
  sn = as.character(seqnames(variantGR))
  stopifnot(length(unique(sn))==1)
  varseq = as.character(seqnames(variantGR)[1])
  gtref = vcfgen(varseq)
  chinds = chunk(1:length(exSE), chunk.size=exChLen)
  stopifnot(length(chinds)>=3)
  aa = AllAssoc(exSE[chinds[[1]],], gtref, variantGR, ...)
  ab = AllAssoc(exSE[chinds[[2]],], gtref, variantGR, ...)
  collec = collapseToBuf(aa, ab, bufsize=bufsize)
  collec_perm = collapseToBuf(aa, ab, bufsize=nperm*bufsize, frag="_permScore")
  for (i in 3:length(chinds)) {
     tmp = AllAssoc(exSE[chinds[[i]],], gtref, variantGR, ...)
     collec = collapseToBuf(collec, tmp, bufsize=bufsize)
     collec_perm = collapseToBuf(collec, tmp, bufsize=nperm*bufsize, frag="_permScore")
  }
  mcols(collec)$perm_scorebuf = mcols(collec_perm)$scorebuf
  collec
}
