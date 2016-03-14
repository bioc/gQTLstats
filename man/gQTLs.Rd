\name{gQTLs}
\alias{gQTLs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use SummarizedExperiment to manage a collection of gQTL results of interest
}
\description{
use SummarizedExperiment to manage a collection of gQTL results of interest
}
\usage{
gQTLs(filtgr, se, tf, genome = "hg19", forceRs = TRUE, chunksize = 50)
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
  # coming
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
