\name{eqBox2}
\alias{eqBox2}
\alias{eqDesc2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualization of expression or other assay measure
against genotypes extracted from VCF
}
\description{
visualization of expression or other assay measure
against genotypes extracted from VCF
}
\usage{
eqBox2(gene, se, tf, snpgr, radius = 1e+06, genome = "hg19", ...)
eqDesc2(gene, se, tf, snpgr, radius = 1e+06, genome = "hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
%%     ~~Describe \code{gene} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{tf}{
%%     ~~Describe \code{tf} here~~
}
  \item{snpgr}{
%%     ~~Describe \code{snpgr} here~~
}
  \item{radius}{
%%     ~~Describe \code{radius} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(Rsamtools)
require(GenomicRanges)
mygr = GRanges("1", IRanges(54683925, width=1))
gene = "ENSG00000231581.1"
library(geuvPack)
data(geuFPKM)
#tf = gtpath(1)
tf = TabixFile(system.file("vcf/small_1.vcf.gz", package="gQTLstats"))
eqBox2(gene, se=geuFPKM, tf, mygr )
eqDesc2(gene, se=geuFPKM, tf, mygr )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }
