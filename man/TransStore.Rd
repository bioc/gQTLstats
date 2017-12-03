\name{TransStore}
\alias{TransStore}
\alias{tsIndex.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Instance constructor for managing trans gQTL results
}
\description{
Instance constructor for managing trans gQTL results
}
\usage{
TransStore(regs, paths = NULL)
tsIndex.reg(tsin, ind)
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
\item{tsin}{a TransStore instance}
\item{ind}{index of registry to index}
}
\details{
\code{tsIndex.reg} is experimental, producing a hash mapping snps
to registry job identifiers, to support rapid store-level retrieval
of locus-specific findings.
}
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
\dontrun{  # dec 2017 changes to BatchJobs lead to errors
if (require(geuvStore2) && require(doParallel)) {
  registerDoSEQ()
  r17 = g17transRegistry()
  r18 = g18transRegistry()
  g1718 = TransStore(list(r17, r18))
  g1718
  }
 } # end dontrun
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
