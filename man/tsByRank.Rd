\name{tsByRank}
\alias{tsByRank}
\alias{tsByRankAccum}
\alias{tsByRankAccum_sing}
\alias{tsByRank_sing}
\alias{getTransRegistries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
harvest contents of a TransStore by rank in associations of features to SNP
}
\description{
Harvest contents of a TransStore by rank in associations of features to SNP.
}
\usage{
tsByRankAccum(tsin, maxrank = 3, mcol2keep=c("REF", "ALT", "snp", "MAF", "z.HWE"), filt=force)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tsin}{
An instance of \code{\link{TransStore-class}}
}
\item{maxrank}{
The maximum rank of association scores to retrieve,
cumulatively.  Each variant has been
tested for association with each genomic feature (e.g.,
gene in a typical expression QTL study), but only the top ranking
associations are recorded for each variant.  If maxrank=k, 
for each variant, this
function retrieves the features exhibiting the kth largest association 
recorded over all features, along with all k-1 larger association scores.
}
\item{mcol2keep}{a character vector of metadata columns to
retain}
\item{filt}{a function accepting a GRanges and returning
a GRanges.  The mcols of the GRanges to be processed
will have elements
\code{c(mcol2keep, "scorebuf", "elnames", "dist")},
where the latter two are matrices with number of columns
equal to the \code{bufsize} of the \code{\link{transAssoc}} call
that generated \code{ts}. Only SNP-specific elements can be used
to define the filter.}
}
\details{
tsByRankAccum_sing and other functions with suffix _sing
were developed for the case of a single permutation

getTransRegistries() accesses objects packaged
for demonstration purposes

}
\value{
A \code{\link[GenomicRanges]{GRanges}} instance.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{ # BatchJObs issues 12/2017
 if (require(doParallel)) {
  registerDoSEQ()
 lit = TransStore(getTransRegistries()) # very limited slice
 tbga = tsByRankAccum(lit, maxrank=5)
 plot(ecdf(as.numeric(data.matrix(tbga$permscoresByRank1))), ylim=c(.99,1),
    main="eCDF of permutation dist. of association, by variant rank")
 exr = paste0("permscoresByRank", 2:5)
 for (i in 1:4)
   lines(ecdf(as.numeric(data.matrix(mcols(tbga)[[exr[i]]]))), col=i+1)
 legend(200, .994, lty=1, col=1:5, legend=paste("rank", 1:5))
 plot(ecdf(as.numeric(data.matrix(tbga$permscoresByRank1[,1]))), ylim=c(.99,1),
    main="between-permutation variation")
 lines(ecdf(as.numeric(data.matrix(tbga$permscoresByRank1[,2]))),col=2)
 lines(ecdf(as.numeric(data.matrix(tbga$permscoresByRank1[,3]))),col=3)
 lines(ecdf(as.numeric(data.matrix(tbga$permscoresByRank5[,1]))),col=4)
 lines(ecdf(as.numeric(data.matrix(tbga$permscoresByRank5[,2]))),col=5)
 lines(ecdf(as.numeric(data.matrix(tbga$permscoresByRank5[,3]))),col=6)
 legend(200, .994, col=1:6, lty=1, legend=c("rank 1 (perm 1)", "(2)", "(3)",
  "rank 5 (perm 1)", "(2)", "(3)"))
# head(tbga,2)
 # consider the following filtering utility
# tbfilt = function(tbg, seqnames.="17", minMAF=.1, minabsodist = 1e7,
#    nrec=1000) {
#   tbg = tbg[ which(as.character(seqnames(tbg)) \%in\% seqnames.) ]
#   tbg = tbg[ which(tbg$MAF > minMAF) ]
#   tbg[ order(tbg$scores, decreasing=TRUE) ][1:nrec]
# }
}
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
