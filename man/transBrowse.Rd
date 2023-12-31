\name{transBrowse}
\alias{transBrowse}
\alias{transBrowse2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
shiny app to exhibit genotype:genomic feature distributions
}
\description{
exhibit genotype:genomic feature distributions with a shiny app
}
\usage{
transBrowse(tbg, anno, tivcf, se, title = "trans eQTL")
transBrowse2(tbga, annovec, tivcf, se, title = "trans eQTL", maxrank=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbg}{
filtered output of \code{\link{tsByRankAccum}}, see example
}
  \item{tbga}{
filtered output of \code{\link{tsByRankAccum}}, see example
}
  \item{anno}{
%#GRanges instance with content analogous to that of 
%#\code{gencodeV12} data element of \code{geuvPack}
a vector with 'feature symbols' (e.g.,, gene symbols) as values and
'feature names' (elements of rownames of \code{se}, e.g., ENSEMBLE gene ids)
as names
}
  \item{annovec}{
%#GRanges instance with content analogous to that of 
%#\code{gencodeV12} data element of \code{geuvPack}
a vector with 'feature symbols' (e.g.,, gene symbols) as values and
'feature names' (elements of rownames of \code{se}, e.g., ENSEMBLE gene ids)
as names
}
  \item{tivcf}{
reference to Tabix-indexed VCF 
}
  \item{se}{
SummarizedExperiment instance with rowname coincident with
\code{anno} and \code{tbg[["allfeats"]]}
}
  \item{title}{
optional string for title panel
}
  \item{maxrank}{
transBrowse2 works with the tsByRankAccum function that collects
scores down to a specified rank.  This parameter specifies the boundary.
}
}
\details{
This function is under development.  The intention is to allow
convenient visualization of off-chromosome genotype-feature
relationships.  AllAssoc collects association scores SNP-wise, and
saves the largest "K" scores obtained, along with feature identity
and location metadata.  The largest score obtained for a given SNP
is the rank 1 association, the next largest is rank 2, and so on.
}
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
\dontrun{
 # consider the following filtering utility
 tbfilt = function(tbg, seqnames.="17", minMAF=.1, minabsodist = 1e7,
    nrec=1000) {
   tbg = tbg[ which(as.character(seqnames(tbg)) \%in\% seqnames.) ]
   tbg = tbg[ which(tbg$MAF > minMAF & abs(tbg$obsdist) > minabsodist) ]
   tbg[ order(tbg$scores, decreasing=TRUE) ][1:nrec]
 }
 #
   registerDoSEQ()
   library(geuvStore2)
   r17 = g17transRegistry()
   g17 = TransStore(list(r17))
   tbg = tbfilt(tsByRankAccum(g17, 3, mcol2keep=c("REF", "snp", "MAF"))) # 1000 records
   tf17 = ldblock::s3_1kg("17")  # uses S3 bucket
   require(geuvPack)
   require(shiny)
   if (!exists("geuFPKM")) data(geuFPKM)
   if (!exists("gencodeV12")) data(gencodeV12)
   data(gen2sym)
   transBrowse2( tbg, gen2sym, tf17, geuFPKM, title="trans GEUV chr17")
 } # end dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
