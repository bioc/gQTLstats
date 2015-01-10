\name{cisAssoc}
\alias{cisAssoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test for variant-expression associations in cis, using VCF
}
\description{
test for variant-expression associations in cis, using VCF
and SummarizedExperiment representations
}
\usage{
cisAssoc(summex, vcf.tf, rhs = ~1, 
   nperm = 3, cisradius = 50000, 
   genome = "hg19", assayind = 1, lbmaf = 1e-06,
   dropUnivHet = TRUE, infoFields = c("LDAF", "SVTYPE"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summex}{
instance of \code{\link[GenomicRanges]{SummarizedExperiment-class}}
}
  \item{vcf.tf}{
instance of \code{\link[Rsamtools]{TabixFile}}, referring to a
tabix-indexed, bgzipped VCF file
}
  \item{rhs}{
formula `right hand side' for adjustments to be made
as \code{\link[snpStats]{snp.rhs.tests}} is run on each
expression vector
}
  \item{nperm}{
number of permutations to be used for plug-in FDR computation
}
  \item{cisradius}{
distance in bp around each gene body to be searched for SNP
association 
}
  \item{genome}{
tag suitable for use in GenomeInfoDb structures
}
  \item{assayind}{
index of \code{assays(summex)} to use for expression data retrieval
}
  \item{lbmaf}{
lower bound on MAF of SNP to use
}
  \item{dropUnivHet}{
logical, if TRUE, will check for columns of SnpMatrix instance
that possess no values other than "NA" and "A/B".
See \url{http://www.biostars.org/p/117155/#117270}
}
  \item{infoFields}{character -- VCF fields to retain in vcfInfo() part of query}
}
\details{
\code{\link[snpStats]{snp.rhs.tests}} is the 
workhorse for statistical modeling.  VCF content is
transformed to the byte-code (which allows for uncertain imputation)
and used in fast testing.
}
\value{
a \code{\link[GenomicRanges]{GRanges-class}} instance
with mcols including chisq, permScore...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjd@channing.harvard.edu>
}
\note{
seqlevelsStyle for summex and vcf.tf content must agree
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
 require(GenomeInfoDb)
 require(geuvPack)
 require(Rsamtools)
 data(geuFPKM)
 lgeu = geuFPKM[ which(seqnames(geuFPKM)=="chr20"), ]
 seqlevelsStyle(lgeu) = "NCBI"
 tf20 = TabixFile(system.file("vcf/c20exch.vcf.gz", package="GGtools"))
 if (require(VariantAnnotation)) scanVcfHeader(tf20)
 lgeue = clipPCs(lgeu[,which(lgeu$popcode=="CEU")], 1:2)
 litc = cisAssoc(lgeue[c(162,201),], tf20, nperm=2, lbmaf=.05, cisradius=50000)
 summary(litc$chisq)
 \dontrun{
 litc$pifdr = gQTLstats:::pifdr(litc$chisq, c(litc$permScore_1, litc$permScore_2))
 litc[which(litc$pifdr < .01)]
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
