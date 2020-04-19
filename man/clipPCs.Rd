\name{clipPCs}
\alias{clipPCs}
\alias{clipPCs,SummarizedExperiment,numeric,logical-method}
\alias{clipPCs,RangedSummarizedExperiment,numeric,logical-method}
\alias{clipPCs,SummarizedExperiment,numeric,missing-method}
\alias{clipPCs,RangedSummarizedExperiment,numeric,missing-method}
\alias{regressOut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transformations of expression data in smlSet instances
}
\description{
transformations of expression data in smlSet instances
or assay data in RangedSummarizedExperiment
}
\usage{
clipPCs(x, inds2drop, center = TRUE)

regressOut(x, rhs, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a \code{\link{RangedSummarizedExperiment}} object
}
  \item{inds2drop}{
Vector of PCs to be eliminated by setting the associated diagonal elements in the
SVD to zero before recomposing the matrix of expression values.  If the value 0
is present in inds2drop, the smlSet is returned unchanged, with a message.
}
  \item{center}{
logical, passed to \code{\link[stats]{prcomp}}
}
  \item{rhs}{formula fragment (no dependent variable) used to form residuals
in a reexpression of the expression matrix; variable bindings
found in pData of an ExpressionSet or colData of a RangedSummarizedExperiment}
  \item{\dots}{arguments passed to \code{\link[limma]{lmFit}}}
}
\details{
\code{clipPCs} is an operation on the n x p transposed matrix X of expression data.
The singular value decomposition X = UDV^t is formed, the diagonal elements of D
corresponding to \code{inds2drop} are set to zero yielding the diagonal matrix E, and
then Y = UEV^t is computed and transposed to replace the expression data.

\code{regressOut} obtains residuals after genewise regression of expression on
the design matrix specified by the \code{rhs}; \code{\link[limma]{lmFit}} is used to
compute coefficients, linear predictions and residuals.
}
\value{
a \code{\link{RangedSummarizedExperiment}} object
}

\references{
The use of PCA-based adjustments to remove mass extraneous effects from
expression matrices has been criticized in work of Oliver Stegle and Jeffrey Leek,
who offer Bayesian PEER and SVA respectively as alternative solutions.
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if(require(geuvPack)){
 data(geuFPKM)
 cg = clipPCs(geuFPKM, 1:10)
 ro = regressOut(cg, ~popcode)
 ro
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
