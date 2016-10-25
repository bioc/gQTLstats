\name{TransStore}
\alias{TransStore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Instance constructor for managing trans gQTL results
}
\description{
Instance constructor for managing trans gQTL results
}
\usage{
TransStore(regs, paths = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regs}{
a list of \code{\link[BatchJobs]{Registry}} instances, typically
one per (variant-oriented) chromosome
}
  \item{paths}{
if desired, paths to folders for which \code{\link[BatchJobs]{loadRegistry}} succeeds, used instead of \code{regs}
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
instance of \code{\link{TransStore-class}}
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
%%  ~~further notes~~
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
  g1718
  }
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
