\name{enumerateByFDR}
\alias{enumerateByFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filter a ciseStore instance using an FDR threshold
}
\description{
filter a ciseStore instance using an FDR threshold
}
\usage{
enumerateByFDR(store, fdrsupp, threshold = 0.05, filter=force,
   ids=NULL, trimToUnit=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link[gQTLBase]{ciseStore-class}}
}
  \item{fdrsupp}{
instance of \code{\link{FDRsupp-class}}
}
  \item{threshold}{
upper bound on FDR to be included
}
  \item{filter}{
The FDR can be computed for any association score.  To
return only records satisfying a given filter, supply
the filter function here.  It may be desirable to 
carry a filter function from the storeToFDR stage, and
this may be considered in future versions.
}
  \item{ids}{if NULL, process all results in store, otherwise
limit attention to jobs with id values in \code{ids}}
  \item{trimToUnit}{plug-in FDR estimates can sometimes
lie outside [0,1] owing to sparsity or defects of extrapolation;
if this parameter is TRUE, estimated FDR values outside [0,1] are moved to
the nearest boundary}
}
\details{
uses \code{\link{storeApply}}, which will
use BiocParallel infrastructure when available
}
\value{
A GRanges instance with store contents to which \code{estFDR}
is appended for each range.  The estFDR quantity is
predicted using the GAM model held in the FDRsupp instance.
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
require(geuvStore)
require(gQTLBase)
re = partialRegistry()
st = ciseStore(re, partialIds(), FALSE, FALSE)
data(filtFDR)
filtEnum = enumerateByFDR( st, filtFDR,
  filter=function(x)x[which(x$mindist <= 500000 & x$MAF >= 0.05)] )
names(metadata(filtEnum))
filtEnum[order(filtEnum$chisq, decreasing=TRUE)[1:2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
