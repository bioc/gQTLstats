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
txsPlot(FDRsupp, xmax=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FDRsupp}{
an instance of \code{\link{FDRsupp-class}}
}
  \item{xmax}{upper bound on xlim for display}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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
data(filtFDR)
txsPlot(filtFDR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }
