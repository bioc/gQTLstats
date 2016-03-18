\name{setFDRfunc}
\alias{setFDRfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimate and store function relating association scores
to approximate plug-in FDR
}
\description{
estimate and store function relating association scores
to approximate plug-in FDR
}
\usage{
setFDRfunc(FDRsupp, fudge = 1e-06, zthresh = 30, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FDRsupp}{
instance of \code{\link{FDRsupp-class}}
}
  \item{fudge}{
if FDR is zero, a log or logistic transform will fail; we add the small
positive number \code{fudge} to avoid this
}
  \item{zthresh}{
for association scores greater than this value, a hard value
of FDR 0 is assigned
}
  \item{\dots}{
arguments passed to \code{\link[mgcv]{s}} for the smooth
model relating association score to FDR at selected quantiles
of the association score distribution
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
returns an updated \code{\link{FDRsupp-class}} instance
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
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(filtFDR)
filtFDR2 = setFDRfunc(filtFDR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
