.validVcfStack = function(object) {
 pa = .paths(object)
 npa = names(pa)
 if(!is.character(npa) || length(npa)!=length(pa)) return("object@paths must be character with one element per path")
 TRUE
}

setClass("VcfStack", representation(paths="character",
   seqinfo="Seqinfo",
   colData="DataFrame"),
   validity=.validVcfStack  # verify 1-1 mapping from paths to seqinfo
)

setClass("RangedVcfStack", contains="VcfStack",
  representation=representation(rowRanges="GRanges", sampleNames="character"))

setGeneric("bindRanges", function(vs, gr, ...) standardGeneric("bindRanges"))
setMethod("bindRanges", c("VcfStack", "GRanges"), function(vs, gr, ...) {
 new("RangedVcfStack", vs, rowRanges=gr)
})
setMethod("bindRanges", c("RangedVcfStack", "GRanges"), function(vs, gr, ...) {
 vs@rowRanges = gr
 vs
})
setGeneric("bindSampleNames", function(vs, sn, ...) standardGeneric("bindSampleNames"))

setMethod("bindSampleNames", c("RangedVcfStack", "character"), function(vs, sn, ...) {
 vs@sampleNames = sn
 vs
})

.paths = function(ob) ob@paths

setMethod("[", c("VcfStack", "GenomicRanges", "character", "missing"),
   function(x, i, j, drop) {
    querseq = as.character(seqnames(i))
    stopifnot(length(unique(querseq))==1) # is this good enough?
    path2use = .paths(x)[querseq]
    param = ScanVcfParam(which=i)
    vcfSamples(param) = j
    readVcf(path2use, param=param, genome=genome(i)[1])
   })

setMethod("[", c("VcfStack", "GenomicRanges", "missing", "missing"),
   function(x, i, j, drop) {
    querseq = as.character(seqnames(i))
    stopifnot(length(unique(querseq))==1) # is this good enough?
    path2use = .paths(x)[querseq]
    param = ScanVcfParam(which=i)
    readVcf(path2use, param=param, genome=genome(i)[1])
   })



setMethod("show", "VcfStack", function(object) {
 cat("VcfStack instance with", length(object@paths), "paths.\n")
 cat("Genome build recorded as ", as.character(genome(object@seqinfo)[1]), ".\n",
     sep="")
 cat("use '[ [GRanges], [sampleids] ]' to extract VariantAnnotation VCF.\n")
 show(object@seqinfo)
})

VcfStack = function(paths, seqinfo, colData=DataFrame(),
   set.seqlstyle="NCBI") {
 pt = paths
 sn = names(paths)
 si = seqinfo[sn]
 seqlevelsStyle(si) = set.seqlstyle
 names(paths) = seqnames(si)
 tmp = new("VcfStack", paths=paths, colData=colData, seqinfo=si)
 tmp
}
 
getVCFPath = function(vs, chrtok) {
  stopifnot(is.atomic(chrtok), length(chrtok)==1)
  .paths(vs)[chrtok]
  }


setMethod("samples", "RangedVcfStack", function(object) {
  if (length(sn <- object@sampleNames)>0) return(sn)
  samples(scanVcfHeader(object@paths[[1]]))
})
setMethod("subsetSample", c("RangedVcfStack"), function(x, j, ...) {
  if (length(sn <- x@sampleNames)>0) {
      x@sampleNames = intersect(j, x@sampleNames)
      return(x)
      }
  sn = samples(scanVcfHeader(x@paths[[1]]))
  x@sampleNames=intersect(j, sn)
  x
})
setMethod("[", c("RangedVcfStack", "GRanges", "ANY", "missing"),
   function(x, i, j, drop) {
     x@rowRanges=i
     if(!missing(j)) x@sampleNames=j
     x
   })  # just updates slots endomorphically
setMethod("[", c("RangedVcfStack", "missing", "missing", "missing"),
   function(x, i, j, drop) {
     x
})  # above needed for biocMultiAssay validity method checker which runs x[]
#
setMethod("features", "RangedVcfStack", function(x) {
  x@rowRanges
})
setMethod("assay", c("RangedVcfStack", "missing"), function(x, i) {
  sn = x@sampleNames
  vp = ScanVcfParam(which=x@rowRanges)
  if (length(sn)>0) rd = x[ x@rowRanges, sn ]
  else rd = x[ x@rowRanges, ]
  mat = as( genotypeToSnpMatrix(rd)$genotypes, "numeric" )
  t(mat)
})
