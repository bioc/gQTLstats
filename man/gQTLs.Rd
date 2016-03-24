\name{gQTLs}
\alias{gQTLs}
\alias{gQTLswarm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use SummarizedExperiment to manage a collection of gQTL results of interest
}
\description{
use SummarizedExperiment to manage a collection of gQTL results of interest
}
\usage{
gQTLs(filtgr, se, tf, genome = "hg19", forceRs = TRUE, chunksize = 50)
gQTLswarm(se, ind, covar = NULL, inpch = 19, xlab, ylab, featTag="probeid", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filtgr}{
a GRanges instance typically obtained by filtering a ciseStore instance
}
  \item{se}{
SummarizedExperiment with individual level expression
and sample-level data from which filtgr statistics were derived
}
  \item{tf}{
TabixFile for VCF on which filtgr statistics are based
}
  \item{genome}{
tag for \code{\link[VariantAnnotation]{readVcf}}
}
  \item{forceRs}{
if TRUE insist that snp ids include 'rs'
}
  \item{chunksize}{
VCF processing proceeds via foreach in chunks of size chunksize
}
  \item{se}{
a SummarizedExperiment instance, output of \code{gQTLs}
}
  \item{ind}{
index into rows of se to be used for visualization, must be length 1
}
  \item{covar}{
a character string indicating a variable in colData(se) to be used
to color the points
}
  \item{inpch}{
pch setting for dots in swarm
}
  \item{xlab}{
xlabel for beeswarm plot, defaults to snp id as recovered from rowRanges(se)$snp
}
  \item{ylab}{
ylabel for beeswarm plot, defaults to probe id as recovered from rowRanges(se)$probeid
}
  \item{featTag}{
element of mcols(rowRanges(se)) used to find ylab text, defaults to 'probeid',
'symbol' is often preferred
}
  \item{\dots}{
passed to \code{\link[beeswarm]{beeswarm}}
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a SummarizedExperiment instance with two assays, the
first is genotype the second is expression
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
very preliminary
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(Rsamtools)
tf = TabixFile(system.file("vcf/litv.vcf.gz", package="gQTLstats"))
data(sigInlit) # 33 loci with significant cis eQTL on a specific filtering
library(geuvPack)
data(geuFPKM)
require(doParallel)
registerDoSEQ()
gdem = gQTLs(sigInlit, geuFPKM, tf, genome = "hg19")
gQTLswarm(gdem, 1, "popcode")
}
\keyword{ manip }
