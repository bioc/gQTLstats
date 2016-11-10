\name{eqBox2}
\alias{eqBox2}
\alias{eqBox3}
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
eqBox2(gene, se, tf, snpgr, genome = "hg19", forceRs=TRUE, ...)
eqDesc2(gene, se, tf, snpgr, genome = "hg19", forceRs=TRUE)
eqBox3(gene, se, tf, snpgr, geneAnno, genome = "hg19", forceRs = TRUE,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
an element of rownames(se) from which a vector of assay
values will be created
}
  \item{se}{
a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object
}
  \item{tf}{
instance of class \code{\link[Rsamtools]{TabixFile-class}}, defining
paths to a tabix-indexed VCF and index file
}
  \item{snpgr}{
instance of \code{\link[GenomicRanges]{GRanges-class}} identifying
the SNP to be visualized
}
  \item{genome}{
tag identifying reference genome
}
  \item{forceRs}{ In the 1000 genomes VCF, there are sometimes
variants identified with DELLY that are grabbed by readVcf on
an SNV address.  Set forceRs to TRUE to retain only variants
with 'rs' in the name.  Has no effect if readVcf extracts only
a single variant.
}
  \item{geneAnno}{
named vector, \code{geneAnno[gene]} will be used to
annotate display
}
  \item{\dots}{
extra arguments passed to beeswarm
}
}
\details{
In 1.5.4, altered to supply beeswarm data visualization in addition
to boxplot.  Use additional option corral="gutter" to reduce 
horizontal sprawl in large samples.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(Rsamtools)
require(SummarizedExperiment)
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
