\name{tqbrowser}
\alias{tqbrowser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
general browsing facility for trans-gQTL
}
\description{
Provide a general browsing facility for trans-gQTL.
}
\usage{
tqbrowser(mae, felname, gelname, tiling, tsbra, 
  annovec, band.init = "6q12", ermaset, gwascat, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mae}{
Instance of \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
}
  \item{felname}{
character naming the element of \code{mae} holding assay features
}
  \item{gelname}{
character naming the element of \code{mae} holding a 
\code{\link[GenomicFiles]{VcfStack-class}} instance for genotypes
}
  \item{tiling}{
a tiling of the genome used to partition large genotype resource
}
  \item{tsbra}{
an instance of the output of \code{\link{tsByRankAccum}} that
collects association statistics and metadata on general searches
for genotype-feature association
}
  \item{annovec}{
a named character vector mapping between identifiers used to identify features in 
\code{experiments(mae)[[felname]]} and tokens to be used
in display -- the names of \code{annovec} are the rownames to be
translated to the associated value in the display.
}
  \item{band.init}{
an initial tile selection
}
  \item{ermaset}{
instance of \code{\link[erma]{ErmaSet-class}}
}
  \item{gwascat}{
instance of \code{\link[gwascat]{gwaswloc-class}}
}
  \item{\dots}{
not currently used
}
}
\details{
starts a shiny app
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (interactive()) {
oa = options()$example.ask
options(example.ask=FALSE)
#
# this example assumes you have a working internet connection
# it will collect genotype information from a S3 bucket
# where 1000 genomes VCF resides
#
# obtain infrastructure
#
# message("note: as of Dec 17 2016 this function will trigger transient errors... ignore them") # solved with req()
message("loading packages...")
packs = c("VariantAnnotation", "GenomicFiles", "ggvis", "plotly",
  "Rsamtools", "GenomeInfoDb", "geuvPack", "shiny", "ggplot2",
  "MultiAssayExperiment", "AnnotationHub", "ldblock", "erma")
suppressPackageStartupMessages({
r = sapply(packs, require, character.only=TRUE)
})
stopifnot(all(r))

# use S3 bucket to get genotypes, create VcfStack wrapper
#
message("create VcfStack...")
tf17 = ldblock::s3_1kg("17")
vcst = VcfStack(c("17"=path(tf17)), index=FALSE)
seqlevelsStyle(vcst) = "NCBI"

# obtain expression data for GEUVADIS samples
#
message("obtain expression data...")
if (!exists("geuFPKM")) data(geuFPKM)
data(gen2sym)
seqlevelsStyle(geuFPKM) = "NCBI"
#
# bind to MAE
#
el = ExperimentList(list(geu=geuFPKM, vcf=vcst))
message("create MultiAssayExperiment...")
suppressWarnings({ # samples don't line up between expression and genotype, we know this
mae = MultiAssayExperiment(el, colData=colData(el[[1]]))
})
#
# obtain and clean up cytoband representation
# cyto37n created as follows:
#ah = AnnotationHub()
#cyto37 = ah[["AH5012"]]
#seqlevelsStyle(cyto37) = "NCBI"
#cyto37 = as(cyto37, "GRanges")
#sn = as.character(seqnames(cyto37))
#mcols(cyto37)$name = paste0(sn, mcols(cyto37)$name)
#names(cyto37) = mcols(cyto37)$name
#seqlengths(cyto37)["MT"] = 16569
message("obtain cytoband index...")
data(cyto37n)
data(tbgaOrmdl3)  # saved output of tsByRankAccum, giving association scores
#
message("obtain gwas catalog...")
library(gwascat)
data(ebicat37)
# obtain chromatin state calls from erma
message("obtain chromatin state calls...")
erset = makeErmaSet()
#
# target and invoke browser
#
okba = c("17q12", "17q21.1", "17q21.2")
on.exit(options(example.ask=oa))
print(tqbrowser( mae, "geu", "vcf", cyto37n[okba], 
   tbgaOrmdl3, gen2sym, band.init="17q12", ermaset=erset, gwascat=ebicat37 ))
} # end interactivity check
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
