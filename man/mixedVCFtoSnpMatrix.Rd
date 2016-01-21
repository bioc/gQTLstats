\name{mixedVCFtoSnpMatrix}
\alias{mixedVCFtoSnpMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
amalgamate called genotypes and imputed allelic dosages in VCF to
SnpMatrix representation
}
\description{
amalgamate called genotypes and imputed allelic dosages in VCF to
SnpMatrix representation
}
\usage{
mixedVCFtoSnpMatrix(vcf, preferGT = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{
object inheriting from \code{\link[VariantAnnotation]{CollapsedVCF-class}}
}
  \item{preferGT}{
logical.  VCF allows loci for samples to be reported in various
formats, and a given locus can have a call tagged GT and a
genotype probability or likelihood representation tagged GP or GL.
\code{\link[VariantAnnotation]{genotypesToSnpMatrix}} has an
\code{uncertain} parameter that, if TRUE, will transform GP or GL
content to allelic dose.  Note that only the "first" dosage type appearing in the header will be transformed.  Thus if GP is first in the header but
a given locus is tagged only with GL, the genotype for
thus locus will be recorded as NA.
}
}
\details{
emulates output from \code{\link[VariantAnnotation]{genotypesToSnpMatrix}}
}
\value{
list with elements \code{genotypes} and \code{map}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[VariantAnnotation]{genotypesToSnpMatrix}}
}
\examples{
fn = system.file("vcf/polytypeSNV.vcf", package="gQTLstats")
require("VariantAnnotation")
require("snpStats")
vv = readVcf(fn, genome="hg19") # only 4th SNP will have dosage coding
mixedVCFtoSnpMatrix(vv)$genotypes@.Data
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
