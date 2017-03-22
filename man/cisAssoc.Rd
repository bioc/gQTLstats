\name{cisAssoc}
\alias{cisAssoc}
\alias{cisEsts}
\alias{AllAssoc}
\alias{cisCount}
\alias{collapse_multiPerm}
\alias{distToGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test for variant-expression associations in cis or generally, using 
VCF
}
\description{
test for variant-expression associations in cis or generally, using VCF
and RangedSummarizedExperiment representations
}
\usage{
cisAssoc(summex, vcf.tf, rhs = ~1, nperm = 3, cisradius =
                 50000, genome = "hg19", assayind = 1, lbmaf = 1e-06,
                 lbgtf = 1e-06, dropUnivHet = TRUE, infoFields =
                 c("LDAF", "SVTYPE"), simpleSNV = TRUE)
cisEsts(summex, vcf.tf, rhs = ~1, nperm = 3, cisradius =
                 50000, genome = "hg19", assayind = 1, lbmaf = 1e-06,
                 lbgtf = 1e-06, dropUnivHet = TRUE, infoFields =
                 c("LDAF", "SVTYPE"), simpleSNV = TRUE)
cisCount(summex, vcf.tf, rhs = ~1, cisradius =
                 50000, genome = "hg19", assayind = 1, lbmaf = 1e-06,
                 lbgtf = 1e-06, dropUnivHet = TRUE, infoFields =
                 c("LDAF", "SVTYPE"), simpleSNV = TRUE)
AllAssoc(summex, vcf.tf, variantRange, rhs = ~1, nperm = 3, 
    genome = "hg19", assayind = 1, lbmaf = 1e-06, lbgtf = 1e-06, 
    dropUnivHet = TRUE, infoFields = c("LDAF", "SVTYPE")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summex}{
a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object
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
lower bound on MAF of SNP to retain for analysis, computed using
\code{\link[snpStats]{col.summary}}
}
  \item{lbgtf}{
lower bound on genotype frequency of SNP to retain for analysis
}
  \item{dropUnivHet}{
logical, if TRUE, will check for columns of SnpMatrix instance
that possess no values other than "NA" and "A/B".
See \url{http://www.biostars.org/p/117155/#117270}
}
  \item{infoFields}{character -- VCF fields to retain in vcfInfo() part of query}
  \item{simpleSNV}{logical -- will use simple computation of
\code{\link[VariantAnnotation]{isSNV}} to filter variants for
analysis to SNV }
  \item{variantRange}{GRanges instance that defines the
scope of the VCF to be used for testing against all features on summex}
}
\details{
\code{\link[snpStats]{snp.rhs.tests}} is the 
workhorse for statistical modeling.  VCF content is
transformed to the byte-code (which allows for uncertain imputation)
and used in fast testing.

distToGene is a helper function that should be replaced with
something from the Bioconductor annotation subsystem
}
\value{
cisAssoc: a \code{\link[GenomicRanges]{GRanges-class}} instance
with mcols including chisq, permScore...

cisCount: enumerate locations in VCF that would be tested
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
#
# obtain geuvadis expression measures as FPKM
#
 data(geuFPKM)
#
# confine the chromosome 20
#
 lgeu = geuFPKM[ which(seqnames(geuFPKM)=="chr20"), ]
 seqlevelsStyle(lgeu) = "NCBI"
#
# acquire subset of genotypes on chr20
#
 tf20 = TabixFile(system.file("vcf/c20exch.vcf.gz", package="gQTLstats"))
 if (require(VariantAnnotation)) scanVcfHeader(tf20)
#
# perform a general technical confounder correction, and confine
# attention to CEU samples
#
 lgeue = clipPCs(lgeu[,which(lgeu$popcode=="CEU")], 1:2)
#
# obtain all score test statistics for SNP:gene pairs at radius 50k
#
 set.seed(1234)
 litc = cisAssoc(lgeue[c(162,201),], tf20, nperm=2, lbmaf=.05, cisradius=50000)
#
# obtain all estimates for SNP:gene pairs at radius 50k
#
 set.seed(1234)
 lite = cisEsts(lgeue[c(162,201),], tf20, nperm=2, lbmaf=.05, cisradius=50000)
 summary(litc$chisq)
 mysr = range(litc)
#
# compute the plug-in FDR
#
 litc$pifdr = gQTLstats:::pifdr(litc$chisq, c(litc$permScore_1, litc$permScore_2))
 litc[which(litc$pifdr < .01)]
#
# trans association testing.  leave to the user the question of
# whether a test is actually cis
#
 lita = AllAssoc(geuFPKM[1:10,], tf20, mysr)
 lita3 = AllAssoc(geuFPKM[11:20,], tf20, mysr)
 #lita5 = AllAssoc(geuFPKM[21:30,], tf20, mysr)
#
# This retains the top 5 (default) associations per SNP
#
 n1 = gQTLstats:::collapseToBuf(lita, lita3)
 #n1 = collapseToBuf(n1, lita5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
