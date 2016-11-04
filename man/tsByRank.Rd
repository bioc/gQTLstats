\name{tsByRank}
\alias{tsByRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
harvest contents of a TransStore by rank in associations of features to SNP
}
\description{
Harvest contents of a TransStore by rank in associations of features to SNP.
}
\usage{
tsByRank(tsin, rank = 1, filt=force)
tsByRankAccum(tsin, maxrank = 3, filt=force)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tsin}{
An instance of \code{\link{TransStore-class}}
}
  \item{rank}{
The rank of association score to retrieve.  Each variant has been
tested for association with each genomic feature (e.g.,
gene in a typical expression QTL study), but only the top ranking
associations are recorded for each variant.  If rank=k, for each variant, this
function retrieves the feature exhibiting the kth largest association 
recorded over all features.
}
\item{maxrank}{all features and scores will be retrieved
up to \code{maxrank}, see example}
\item{filt}{a function accepting a GRanges and returning
a GRanges.  The mcols of the GRanges to be processed
will have elements
"REF", "ALT", "snp", "MAF", "scorebuf", "elnames", "dist",
where the latter two are matrices with number of columns
equal to the \code{bufsize} of the \code{\link{transAssoc}} call
that generated \code{ts}. Only SNP-specific elements can be used
to define the filter.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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
 if (require(geuvStore2) && require(doParallel)) {
  registerDoSEQ()
  r17 = g17transRegistry()
  r18 = g18transRegistry()
  g1718 = TransStore(list(r17, r18))
  tbg = tsByRank(g1718, 1)

 # consider the following filtering utility
 tbfilt = function(tbg, seqnames.="17", minMAF=.1, minabsodist = 1e7,
    nrec=1000) {
   tbg = tbg[ which(as.character(seqnames(tbg)) \%in\% seqnames.) ]
   tbg = tbg[ which(tbg$MAF > minMAF) ] 
   tbg[ order(tbg$scores, decreasing=TRUE) ][1:nrec]
 }
 #
 g17 = TransStore(list(r17))
 tbga = tbfilt(tsByRankAccum(g17, 3)) # 1000 records
 head(tbga,2)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
