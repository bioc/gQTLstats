\name{tsByRank}
\alias{tsByRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
harvest contents of a TransStore by rank in associations of features to SNP
}
\description{
Harvest contents of a TransStore by rank in associations of features to SNP.
}
\usage{
tsByRank(ts, rank = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ts}{
An instance of \code{\link{TransStore-class}}
}
  \item{rank}{
The rank of association score to retrieve.  Each variant has been
tested for association with each genomic feature (e.g.,
gene in a typical expression QTL study), but only the top ranking
associations are recorded for each variant.  If rank=k, for each variant, this
function retrieves the feature exhibiting the kth largest association 
recorded over all features.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \code{\link[GenomicRanges]{GRanges}} instance.
}
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
 if (require(geuvStore2) && require(doParallel)) {
  registerDoSEQ()
  r17 = g17transRegistry()
  r18 = g18transRegistry()
  g1718 = TransStore(list(r17, r18))
  tsByRank(g1718, 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
