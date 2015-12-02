setGeneric("clipPCs", 
 function(x, inds2drop, center=TRUE) standardGeneric("clipPCs"))
 
#setMethod("clipPCs", 
#  c("SummarizedExperiment", "numeric", "logical"), function(x, inds2drop, center=TRUE){
#   .clipPCs.SE(se=x, inds2drop, center)
#})

setMethod("clipPCs", 
  c("RangedSummarizedExperiment", "numeric", "logical"), function(x, inds2drop, center=TRUE){
   .clipPCs.SE(se=x, inds2drop, center)
})

#setMethod("clipPCs", 
#  c("SummarizedExperiment", "numeric", "missing"), function(x, inds2drop, center=TRUE){
#   .clipPCs.SE(se=x, inds2drop, TRUE)
#})

setMethod("clipPCs", 
  c("RangedSummarizedExperiment", "numeric", "missing"), function(x, inds2drop, center=TRUE){
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
     cn = colnames(se)
     rn = rownames(se)
     assn = names(assays(se))
     message(paste("clipping PCs", 
          paste0(selectSome(inds2drop),collapse=","), "from", assn[1], collapse=""))
     ex = t(assays(se)[[1]])
     recon = reconstruct(ex, inds2drop, center)
     assays(se)[[1]] = recon   # BAD ... assay(x) seems to do better
     metadata(se)$PCsClipped = inds2drop
     colnames(se)=cn
     rownames(se)=rn
     se
     }

regressOut = function(x, rhs, ...) {
 stopifnot(is(x, "SummarizedExperiment") ||
           is(x, "RangedSummarizedExperiment")) 
 stopifnot(is(rhs, "formula")) 
 rn = rownames(x)
 cn = colnames(x)
 mm = model.matrix(rhs, data=colData(x))
 message("using assay() to extract 'expression' matrix from RangedSummarizedExperiment")
 ex = assay(x)
 f = lmFit(ex, mm, ...)
 r = ex - (f$coef %*% t(f$design))
 assay(x) = r
 rownames(x) = rn
 colnames(x) = cn
 x
}

chantmpl = function() "/proj/rerefs/reref00/1000Genomes/Phase3_v5/ALL/ALL.%%N%%.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"

