\name{hmm878}
\alias{hmm878}
\docType{data}
\title{
labeled GRanges with ChromHMM chromatin states for GM12878
}
\description{
labeled GRanges with ChromHMM chromatin states for GM12878
}
\usage{data(hmm878)}
\format{
  The format is: \cr
Formal class 'GRanges' [package "GenomicRanges"] with 6 slots \cr
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. ..@ values         : Factor w/ 23 levels "chr1","chr2",..: 1 2 3 4 5 6 7 8 9 10 ... \cr
  .. .. ..@ lengths        : int [1:23] 54467 46499 37617 25155 30071 34846 29420 24506 24123 27263 ... \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. ..@ start          : int [1:571339] 10001 10601 11138 11738 11938 12138 14538 20338 22138 22938 ... \cr
  .. .. ..@ width          : int [1:571339] 600 537 600 200 200 2400 5800 1800 800 4000 ... \cr
  .. .. ..@ NAMES          : NULL \cr
  .. .. ..@ elementType    : chr "integer" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots \cr
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3 \cr
  .. .. ..@ lengths        : int 571339 \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. ..@ rownames       : NULL \cr
  .. .. ..@ nrows          : int 571339 \cr
  .. .. ..@ listData       :List of 4 \cr
  .. .. .. ..$ name   : chr [1:571339] "15_Repetitive/CNV" "13_Heterochrom/lo" "8_Insulator" "11_Weak_Txn" ... \cr
  .. .. .. ..$ score  : num [1:571339] 0 0 0 0 0 0 0 0 0 0 ... \cr
  .. .. .. ..$ itemRgb: chr [1:571339] "#F5F5F5" "#F5F5F5" "#0ABEFE" "#99FF66" ... \cr
  .. .. .. ..$ thick  :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. .. ..@ start          : int [1:571339] 10001 10601 11138 11738 11938 12138 14538 20338 22138 22938 ... \cr
  .. .. .. .. .. ..@ width          : int [1:571339] 600 537 600 200 200 2400 5800 1800 800 4000 ... \cr
  .. .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ elementType    : chr "ANY" \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ metadata       : list() \cr
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots \cr
  .. .. ..@ seqnames   : chr [1:23] "chr1" "chr2" "chr3" "chr4" ... \cr
  .. .. ..@ seqlengths : int [1:23] 249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 ... \cr
  .. .. ..@ is_circular: logi [1:23] FALSE FALSE FALSE FALSE FALSE FALSE ... \cr
  .. .. ..@ genome     : chr [1:23] "hg19" "hg19" "hg19" "hg19" ... \cr
  ..@ metadata       :List of 1 \cr
  .. ..$ url: chr "http://genome.ucsc.edu/cgi-bin/hgFileUi?g=wgEncodeBroadHmm&db=hg19" \cr
}
\details{
acquired using rtracklayer import from the bed file given at metadata(hmm878)[["url"]]
}
\source{
see details
}
\references{
Ernst J, Kellis M. Discovery and characterization of chromatin states for systematic annotation of the human genome. Nat Biotechnol. 2010 Aug;28(8):817-25.

Ernst J, Kheradpour P, Mikkelsen TS, Shoresh N, Ward LD, Epstein CB, Zhang X, Wang L, Issner R, Coyne M et al. Mapping and analysis of chromatin state dynamics in nine human cell types. Nature. 2011 May 5;473(7345):43-9.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(hmm878)
table(hmm878$name)
}
\keyword{datasets}
