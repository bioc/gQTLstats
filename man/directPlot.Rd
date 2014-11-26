\name{directPlot}
\alias{directPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualize relationship between empirical and
modeled FDR based on analysis of a gQTL store
}
\description{
visualize relationship between empirical and
modeled FDR based on analysis of a gQTL store
}
\usage{
directPlot(FDRsupp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FDRsupp}{
instance of \code{\link[gQTLstats]{FDRsupp-class}}
}
}
\details{
This plot is used to show the degree
of fit between a smooth model
relating modeled FDR to empirical FDR,
and the empirical FDR themselves.
It should be used in conjunction with \code{\link{txsPlot}}.

It is possible for an implausible squiggly model
to yield perfect agreement for all empirical FDR
estimates.  See the example.
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
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (FDRsupp) 
{
    stopifnot(is.function(FDRsupp@FDRfunc))
    tab = getTab(FDRsupp)
    plot(tab$fdr, getFDRfunc(FDRsupp)(tab$assoc), xlab = "empirical FDR", 
        ylab = "modeled FDR")
    abline(0, 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
