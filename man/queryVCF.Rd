\name{queryVCF}
\alias{queryVCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
obtain SnpMatrix from VCF genotypes
}
\description{
obtain SnpMatrix from VCF genotypes
}
\usage{
queryVCF(gr, vcf.tf, samps, genome = "hg19", getSM = TRUE,
   snvOnly=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
GRanges instance; SNPs lying within will be processed
}
  \item{vcf.tf}{
TabixFile instance pointing to VCF
}
  \item{samps}{
samples to be retained
}
  \item{genome}{
tag identifying build
}
  \item{getSM}{
logical; if FALSE, \code{\link[VariantAnnotation]{genotypeToSnpMatrix}}
will not be run and only the output of readVcf is returned.
}
  \item{snvOnly}{logical, if TRUE, will confine results to SNV}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list of length two 
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{readout }{output of readVcf}
  \item{sm }{output of genotypeToSnpMatrix run on the read result}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(Rsamtools)
tf20 = TabixFile(system.file("vcf/c20exch.vcf.gz", package="GGtools"))
require(geuvPack)
data(geuFPKM)
lgeu = geuFPKM[ which(seqnames(geuFPKM)=="chr20"), 
    which(geuFPKM$popcode=="CEU") ]
seqlevelsStyle(lgeu) = "NCBI"
rng = rowRanges(lgeu)[232] # CPNE1
myq = queryVCF( rng, tf20, samps=colnames(lgeu), genome="hg19" )
myq
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
