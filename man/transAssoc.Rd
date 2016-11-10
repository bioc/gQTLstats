\name{transAssoc}
\alias{transAssoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute 'trans' SNP-feature associations by wrapping AllAssoc
}
\description{
compute 'trans' SNP-feature associations by wrapping AllAssoc,
retaining only the strongest associations (and similarly
filtered association scores computed under permutation)
}
\usage{
transAssoc(variantGR, exSE, vcfgen, bufsize = 10, nperm = 3, exChLen = 2 * bufsize, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variantGR}{
GRanges instance establishing scope of variants to test
}
  \item{exSE}{
SummarizedExperiment instance, all of whose features will
be tested for association with all SNP
}
  \item{vcfgen}{
a function returning a path to a tabix-indexed VCF file
from which SNP genotypes will be extracted
}
  \item{bufsize}{
Size of 'buffer' used to retain largest feature association scores
encountered during the search.  The scores and the names of
associated genes are retained in 'scorebuf' and 'elnames' components
of output GRanges
}
  \item{nperm}{
number of permutations of features against genotypes to be
performed for realizing null distribution of association scores
}
  \item{exChLen}{
size of chunks of exSE to be tested through calls to AllAssoc; this
is intended to allow control of RAM usage
}
  \item{\dots}{
arguments passed to AllAssoc
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a GRanges with mcols including
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{  # requires acccess to 1KG S3
library(geuvPack)
data(geuFPKM)
seqlevelsStyle(geuFPKM) = "NCBI"
mysr = GRanges("20", IRanges(33000055, 33020055))
genome(mysr) = "hg19"
tt = transAssoc(mysr, geuFPKM[1:16,],
    bufsize=3, exChLen=4, vcfgen=function(x)gtpath(paste0("chr", x)) )
colnames(mcols(tt))
table(as.character(mcols(tt)$elnames))
}
}
\keyword{ models }
