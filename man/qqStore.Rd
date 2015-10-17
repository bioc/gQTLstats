\name{qqStore}
\alias{qqStore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a binned QQplot for a sharded store
}
\description{
create a binned QQplot for a sharded store with association and permutation statistics
}
\usage{
qqStore(st, ids = NULL, 
  .probs = c(0, seq(0.6, 0.8, 0.2), 0.9, 0.95, 0.99, 0.999, 0.9999, 1), 
  xlim.in = c(0.2, 75), lowfac = 0.5, xlab = "Permutation distribution", 
  ylab = "Distribution of score statistic", countpos = 50, 
  plot.it = TRUE, doab = TRUE, scoreField = "chisq", 
  permField = "permScore_1", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{st}{
instance of \code{\link[gQTLBase]{ciseStore-class}}
}
  \item{ids}{
optional job id vector; if NULL, all jobs used
}
  \item{.probs}{
vector of probabilities for use with quantile evaluation,
as provided in ffbase, using \code{\link[gQTLstats]{storeToQuantiles}}
}
  \item{xlim.in}{
xlim setting for QQplot
}
  \item{lowfac}{
we use a log-log plot, and the first quantile (as prescribed
in \code{.probs}) is often close to
zero; we reassign it to lowfac*(second quantile)
}
  \item{xlab}{
label
}
  \item{ylab}{
label
}
  \item{countpos}{
where on the x axis will we stack the information on bin counts
}
  \item{plot.it}{
logical, if FALSE, a list is returned with elements on
quantile values and bin counts
}
  \item{doab}{
logical prescribing drawing of line of identity
}
  \item{scoreField}{
tag in store naming the statistic, typically 'chisq', can also
be 'tstat' for GTEx
}
  \item{permField}{
tag in store naming the field holding statistics on
realizations from permutation distribution
}
  \item{\dots}{
passed to storeToQuantiles
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
invisibly returns list with
 elements qx, qy, counts, fracs
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
\dontrun{
library(geuvStore)
library(gQTLBase)
gs = makeGeuvStore()
qqStore(gs, ids=partialIds()[1:20])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
