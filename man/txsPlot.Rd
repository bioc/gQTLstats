\name{txsPlot}
\alias{txsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualize transformed FDR against transformed association statistics
}
\description{
visualize transformed FDR against transformed association statistics
}
\usage{
txsPlot(FDRsupp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FDRsupp}{
%%     ~~Describe \code{FDRsupp} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (FDRsupp) 
{
    stopifnot(is.function(FDRsupp@FDRfunc))
    tab = getTab(FDRsupp)
    plot(qlogis(fdr) ~ log(assoc), data = getTab(FDRsupp))
    lines(seq(-6, 5, 0.1), predict(FDRsupp@FDRmodel, newdata = list(assoc = exp(seq(-6, 
        5, 0.1)))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
