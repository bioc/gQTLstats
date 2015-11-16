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

.paths = function(ob) ob@paths

setMethod("[", c("VcfStack", "GenomicRanges", "character", "missing"),
   function(x, i, j, drop) {
    querseq = as.character(seqnames(i))
    stopifnot(length(unique(querseq))==1)
    path2use = .paths(x)[querseq]
    param = ScanVcfParam(which=i)
    vcfSamples(param) = j
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
 
