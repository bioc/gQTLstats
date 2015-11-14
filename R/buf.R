.collapseToBuf.raw = function(x, y, bufsize, frag) {
  stopifnot(length(y) == length(x))
  nsnp = length(x)
  obsinds = function(x) grep(frag, colnames(mcols(x)))
  obsmat = function(x) data.matrix(as.data.frame(mcols(x)[, obsinds(x)]))
  fnames = colnames(mcols(x))[obsinds(x)]
  fnames.y = colnames(mcols(y))[obsinds(y)]
  fnames = c(fnames, fnames.y)
  fnames = gsub(paste0(frag, ".*"), "", fnames)
  XO = obsmat(x)
  YO = obsmat(y)
  XYO = cbind(XO, YO)
  o = t(apply(XYO,1,order,decreasing=TRUE)) # nrow = nsnp
  obufinds = o[, 1:bufsize]
  dat = sapply(1:nsnp, function(z) XYO[z, obufinds[z,]])
  topn = sapply(1:nsnp, function(z) fnames[obufinds[z,]])
  cl = x
  mcols(cl) = mcols(cl)[,c("REF", "ALT", "snp", "MAF")]
  mcols(cl)$scorebuf = t(dat)
  mcols(cl)$elnames = t(topn)
  cl
}
  
collapseToBuf = function(x, y, bufsize=5, frag="_obs") {
  nx = names(mcols(x))
  ny = names(mcols(y))
  if ((!("scorebuf" %in% nx)) & !("scorebuf" %in% ny))
     return(.collapseToBuf.raw(x, y, bufsize, frag))
  if (("scorebuf" %in% nx) & (!("scorebuf" %in% ny)))
     return(.updateBuf(x, y, bufsize, frag))
  if ((!("scorebuf" %in% nx)) & ("scorebuf" %in% ny))
     return(.updateBuf(y, x, bufsize, frag))
  stop("don't see how to start or update the inputs")
}

.updateBuf = function(x, y, bufsize, frag) {
# x is GRanges with scorebuf and elnames
# y is GRanges with [probename]_obs and [probename]_permScore[n]
  nsnp = length(x)
  stopifnot(length(y) == nsnp)
  bufmat = mcols(x)[,"scorebuf"]
  obsinds = function(x) grep(frag, colnames(mcols(x)))
  obsmat = function(x) data.matrix(as.data.frame(mcols(x)[, obsinds(x)]))
  ymat = obsmat(y)
  ynames = colnames(ymat)
  ynames = gsub(paste0(frag, ".*"), "", ynames)
  ynmat = t(matrix(rep(ynames, nrow(ymat)), nc=nrow(ymat)))
  XYO = cbind(bufmat, ymat)
  alln = cbind(mcols(x)[,"elnames"], ynmat)
  o = t(apply(XYO,1,order,decreasing=TRUE)) # nrow = nsnp
  obufinds = o[, 1:bufsize]
  dat = sapply(1:nsnp, function(z) XYO[z, obufinds[z,]])
  topn = sapply(1:nsnp, function(z) alln[z, obufinds[z,]])
  cl = x
  mcols(cl) = mcols(cl)[,c("REF", "ALT", "snp", "MAF")] # restart
  mcols(cl)$scorebuf = t(dat)
  mcols(cl)$elnames = t(topn)
  cl
}
