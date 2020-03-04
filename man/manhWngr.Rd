\name{manhWngr}
\alias{manhWngr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
manhattan plot with named GRanges
}
\description{
manhattan plot with named GRanges
}
\usage{
manhWngr(store, probeid = "ENSG00000183814.10", sym = "LIN9", fdrsupp, namedGR, slstyle = "NCBI", xlab.in = sym, ylab.in = "-log10 FDR", applyFDRfilter = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link[gQTLBase]{ciseStore-class}}
}
  \item{probeid}{
name of feature identifier to use for cis association
}
  \item{sym}{
symbol for feature identifier
}
  \item{fdrsupp}{
instance of \code{\link{FDRsupp-class}}
}
  \item{namedGR}{
\code{\link[GenomicRanges]{GRanges}} instance with 'name' in mcols element
}
  \item{slstyle}{
seqlevelsStyle
}
  \item{xlab.in}{
x axis label
}
  \item{ylab.in}{
y axis label
}
  \item{applyFDRfilter}{
if TRUE, use the filter defined in the \code{filterUsed}
element of the object supplied as \code{fdrsupp} 
on the output
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(geuvStore2)
require(gQTLBase)
store = makeGeuvStore2()
data(hmm878)
data(filtFDR)
manhWngr(store, fdrsupp=filtFDR, namedGR=hmm878)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
