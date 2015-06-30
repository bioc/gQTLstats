\name{filtFDR}
\alias{filtFDR}
\alias{rawFDR}
\docType{data}
\title{
illustration of FDRsupp class
}
\description{
illustration of FDRsupp class
}
\usage{data("filtFDR")}
\format{
  A FDRsupp object.
}
\details{
filtFDR was constructed on geuvStore contents, filtering to
MAF at least five percent and radius at most 500kbp.
rawFDR uses the entire geuvStore contents, with 1Mbp radius and
1 percent MAF lower bound
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(filtFDR)
filtFDR
}
\keyword{datasets}
