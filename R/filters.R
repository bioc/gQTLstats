setGeneric("clipPCs", 
 function(x, inds2drop, center=TRUE) standardGeneric("clipPCs"))
 
setMethod("clipPCs", 
  c("SummarizedExperiment", "numeric", "logical"), function(x, inds2drop, center=TRUE){
   .clipPCs.SE(se=x, inds2drop, center)
})

setMethod("clipPCs", 
  c("SummarizedExperiment", "numeric", "missing"), function(x, inds2drop, center=TRUE){
   .clipPCs.SE(se=x, inds2drop, TRUE)
})


reconstruct = function(ex, inds2drop, center=TRUE) {
    ex = scale(ex, center=center, scale = FALSE)
    ss = svd(ex)
    d = ss$d
    d[inds2drop] = 0
    t(ss$u %*% diag(d) %*% t(ss$v))
}

.clipPCs.SE = function(se, inds2drop, center=TRUE) {
     assn = names(assays(se))
     message(paste("clipping PCs", 
          paste0(selectSome(inds2drop),collapse=","), "from", assn[1], collapse=""))
     ex = t(assays(se)[[1]])
     recon = reconstruct(ex, inds2drop, center)
     assays(se)[[1]] = recon
     exptData(se)$PCsClipped = inds2drop
     se
     }

regressOut = function(sms, x, ...) {
 stopifnot(is(x, "SummarizedExperiment")) 
 mm = model.matrix(rhs, data=colData(x))
 message("using assay() to extract 'expression' matrix from SummarizedExperiment")
 ex = assay(x)
 f = lmFit(ex, mm, ...)
 r = ex - (f$coef %*% t(f$design))
 assay(x) = r
 x
}

