\name{filtFDR}
\alias{filtFDR}
\alias{rawFDR}
\docType{data}
\title{
illustration of FDRsupp class
}
\description{
illustration of FDRsupp class
}
\usage{data("filtFDR")}
\format{
  The format is:
Formal class 'FDRsupp' [package "gQTLstats"] with 5 slots \cr
  ..@ tab     :'data.frame':	27 obs. of  4 variables: \cr
  .. ..$ assoc    : num [1:27] 0.00403 0.01632 0.03726 0.06758 0.10704 ... \cr
  .. ..$ fdr      : num [1:27] 0.946 0.942 0.937 0.932 0.925 ... \cr
  .. ..$ ncalls   : num [1:27] 1567862 1485343 1402824 1320305 1237786 ... \cr
  .. ..$ avg.false: num [1:27] 1483742 1399124 1314093 1230747 1145307 ... \cr
  ..@ FDRfunc :function (assoc)   \cr
  ..@ FDRmodel:List of 34 \cr
  .. ..$ smooth.frame       :'data.frame':	27 obs. of  1 variable: \cr
  .. .. ..$ lo(assoc, ...): smooth [1:27, 1] 0.00403 0.01632 0.03726 0.06758 0.10704 ... \cr
  .. .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. .. ..$ : NULL \cr
  .. .. .. .. ..$ : chr "assoc" \cr
  .. .. .. ..- attr(*, "degree")= num 1 \cr
  .. .. .. ..- attr(*, "span")= num 0.25 \cr
  .. .. .. ..- attr(*, "ncols")= num 1 \cr
  .. .. .. ..- attr(*, "call")= language gam.lo(data[["lo(assoc, ...)"]], z, w, span = 0.25, degree = 1, ncols = 1) \cr
  .. .. .. ..- attr(*, "class")= chr [1:2] "smooth" "matrix" \cr
  .. ..$ coefficients       : Named num [1:2] -0.1472 -0.0713 \cr
  .. .. ..- attr(*, "names")= chr [1:2] "(Intercept)" "lo(assoc, ...)" \cr
  .. ..$ residuals          : Named num [1:27] 0.016 -0.0246 -0.0445 -0.0186 -0.016 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. ..$ fitted.values      : Named num [1:27] 2.85 2.81 2.74 2.64 2.53 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. ..$ effects            : Named num [1:27] 14.6928 -27.7869 -0.045 -0.0191 -0.0166 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "(Intercept)" "lo(assoc, ...)" "" "" ... \cr
  .. ..$ weights            : Named num [1:27] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. ..$ rank               : int 2 \cr
  .. ..$ assign             : int [1:2] 0 1 \cr
  .. ..$ qr                 :List of 5 \cr
  .. .. ..$ qr   : num [1:27, 1:2] -5.196 0.192 0.192 0.192 0.192 ... \cr
  .. .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. .. ..$ : chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. .. .. .. ..$ : chr [1:2] "(Intercept)" "lo(assoc, ...)" \cr
  .. .. .. ..- attr(*, "assign")= int [1:2] 0 1 \cr
  .. .. ..$ qraux: num [1:2] 1.19 1.08 \cr
  .. .. ..$ rank : int 2 \cr
  .. .. ..$ pivot: int [1:2] 1 2 \cr
  .. .. ..$ tol  : num 1e-07 \cr
  .. .. ..- attr(*, "class")= chr "qr" \cr
  .. ..$ smooth             : num [1:27, 1] 3 2.96 2.89 2.79 2.69 ... \cr
  .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. ..$ : chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. .. .. ..$ : chr "lo(assoc, ...)" \cr
  .. ..$ nl.df              : Named num 10.6 \cr
  .. .. ..- attr(*, "names")= chr "lo(assoc, ...)" \cr
  .. ..$ df.residual        : num 14.4 \cr
  .. ..$ var                : num [1:27, 1] 0.454 0.291 0.183 0.278 0.315 ... \cr
  .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. ..$ : chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. .. .. ..$ : chr "lo(assoc, ...)" \cr
  .. ..$ additive.predictors: Named num [1:27] 2.85 2.81 2.74 2.64 2.53 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. ..$ R                  : num [1:2, 1:2] -5.2 0 -195.4 389.8 \cr
  .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. ..$ : chr [1:2] "(Intercept)" "lo(assoc, ...)" \cr
  .. .. .. ..$ : chr [1:2] "(Intercept)" "lo(assoc, ...)" \cr
  .. ..$ rank               : int 2 \cr
  .. ..$ family             :List of 11 \cr
  .. .. ..$ family    : chr "gaussian" \cr
  .. .. ..$ link      : chr "identity" \cr
  .. .. ..$ linkfun   :function (mu)   \cr
  .. .. ..$ linkinv   :function (eta)   \cr
  .. .. ..$ variance  :function (mu)   \cr
  .. .. ..$ dev.resids:function (y, mu, wt)   \cr
  .. .. ..$ aic       :function (y, n, mu, wt, dev)   \cr
  .. .. ..$ mu.eta    :function (eta)   \cr
  .. .. ..$ initialize:  expression({     n <- rep.int(1, nobs)     if (is.null(etastart) && is.null(start) && is.null(mustart) && ((family$link ==          "inverse" && any(y == 0)) || (family$link == "log" && any(y <=          0))))          stop("cannot find valid starting values: please specify some")     mustart <- y }) \cr
  .. .. ..$ validmu   :function (mu)   \cr
  .. .. ..$ valideta  :function (eta)   \cr
  .. .. ..- attr(*, "class")= chr "family" \cr
  .. ..$ deviance           : num 0.158 \cr
  .. ..$ aic                : num -35.1 \cr
  .. ..$ null.deviance      : num 1236 \cr
  .. ..$ iter               : int 2 \cr
  .. ..$ prior.weights      : num [1:27] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ y                  : Named num [1:27] 2.87 2.79 2.7 2.62 2.52 ... \cr
  .. .. ..- attr(*, "names")= chr [1:27] "5%" "10%" "15%" "20%" ... \cr
  .. ..$ df.null            : int 26 \cr
  .. ..$ nl.chisq           : Named num 464 \cr
  .. .. ..- attr(*, "names")= chr "lo(assoc, ...)" \cr
  .. ..$ call               : language gam(formula = qlogis(fdr + fudge) ~ lo(assoc, ...), data = getTab(FDRsupp)) \cr
  .. ..$ formula            :Class 'formula' length 3 qlogis(fdr + fudge) ~ lo(assoc, ...) \cr
  .. .. .. ..- attr(*, ".Environment")=<environment: 0x7fc2a2b75350>  \cr
  .. ..$ terms              :Classes 'terms', 'formula' length 3 qlogis(fdr + fudge) ~ lo(assoc, ...) \cr
  .. .. .. ..- attr(*, "variables")= language list(qlogis(fdr + fudge), lo(assoc, ...)) \cr
  .. .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1 \cr
  .. .. .. .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. .. .. .. ..$ : chr [1:2] "qlogis(fdr + fudge)" "lo(assoc, ...)" \cr
  .. .. .. .. .. ..$ : chr "lo(assoc, ...)" \cr
  .. .. .. ..- attr(*, "term.labels")= chr "lo(assoc, ...)" \cr
  .. .. .. ..- attr(*, "specials")=Dotted pair list of 3 \cr
  .. .. .. .. ..$ s     : NULL \cr
  .. .. .. .. ..$ lo    : int 2 \cr
  .. .. .. .. ..$ random: NULL \cr
  .. .. .. ..- attr(*, "order")= int 1 \cr
  .. .. .. ..- attr(*, "intercept")= int 1 \cr
  .. .. .. ..- attr(*, "response")= int 1 \cr
  .. .. .. ..- attr(*, ".Environment")=<environment: 0x7fc2a2b75350>  \cr
  .. ..$ data               :'data.frame':	27 obs. of  4 variables: \cr
  .. .. ..$ assoc    : num [1:27] 0.00403 0.01632 0.03726 0.06758 0.10704 ... \cr
  .. .. ..$ fdr      : num [1:27] 0.946 0.942 0.937 0.932 0.925 ... \cr
  .. .. ..$ ncalls   : num [1:27] 1567862 1485343 1402824 1320305 1237786 ... \cr
  .. .. ..$ avg.false: num [1:27] 1483742 1399124 1314093 1230747 1145307 ... \cr
  .. ..$ offset             : NULL \cr
  .. ..$ control            :List of 5 \cr
  .. .. ..$ epsilon   : num 1e-07 \cr
  .. .. ..$ maxit     : num 30 \cr
  .. .. ..$ bf.epsilon: num 1e-07 \cr
  .. .. ..$ bf.maxit  : num 30 \cr
  .. .. ..$ trace     : logi FALSE \cr
  .. ..$ method             : chr "glm.fit" \cr
  .. ..$ contrasts          : NULL \cr
  .. ..$ xlevels            : Named list() \cr
  .. ..- attr(*, "class")= chr [1:3] "gam" "glm" "lm" \cr
  ..@ theCall : language storeToFDR(store = prst, xprobs = c(seq(0.05, 0.95, 0.05), 0.975, 0.99,      0.995, 0.9975, 0.999, 0.9995, 0.9999, 0.99999), filter = function(x) x[which(x$MAF >=  ... \cr
  ..@ sessinfo:List of 6 \cr
  .. ..$ R.version :List of 14 \cr
  .. .. ..$ platform      : chr "x86_64-apple-darwin13.4.0" \cr
  .. .. ..$ arch          : chr "x86_64" \cr
  .. .. ..$ os            : chr "darwin13.4.0" \cr
  .. .. ..$ system        : chr "x86_64, darwin13.4.0" \cr
  .. .. ..$ status        : chr "Under development (unstable)" \cr
  .. .. ..$ major         : chr "3" \cr
  .. .. ..$ minor         : chr "2.0" \cr
  .. .. ..$ year          : chr "2014" \cr
  .. .. ..$ month         : chr "10" \cr
  .. .. ..$ day           : chr "14" \cr
  .. .. ..$ svn rev       : chr "66761" \cr
  .. .. ..$ language      : chr "R" \cr
  .. .. ..$ version.string: chr "R Under development (unstable) (2014-10-14 r66761)" \cr
  .. .. ..$ nickname      : chr "Unsuffered Consequences" \cr
  .. ..$ platform  : chr "x86_64-apple-darwin13.4.0 (64-bit)" \cr
  .. ..$ locale    : chr "en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8" \cr
  .. ..$ basePkgs  : chr [1:10] "stats4" "parallel" "stats" "graphics" ... \cr
  .. ..$ otherPkgs :List of 17 \cr
  .. .. ..$ gQTLstats    :List of 13 \cr
  .. .. .. ..$ Package        : chr "gQTLstats" \cr
  .. .. .. ..$ Title          : chr "gQTLstats: computationally efficient analysis for eQTL and allied studies" \cr
  .. .. .. ..$ Version        : chr "0.0.14" \cr
  .. .. .. ..$ Author         : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ Description    : chr "computationally efficient analysis of eQTL, mQTL, dsQTL, etc." \cr
  .. .. .. ..$ Suggests       : chr "geuvPack, geuvStore, Rsamtools, knitr, rmarkdown" \cr
  .. .. .. ..$ Depends        : chr "" \cr
  .. .. .. ..$ Imports        : chr "snpStats, GenomicRanges, VariantAnnotation, methods, Biobase,BiocParallel, BatchJobs, gQTLBase, GenomeInfoDb, S4Vectors, limm"| __truncated__ \cr
  .. .. .. ..$ Maintainer     : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad       : chr "yes" \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-26 13:08:50 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/gQTLstats/Meta/package.rds" \cr
  .. .. ..$ gQTLBase     :List of 14 \cr
  .. .. .. ..$ Package        : chr "gQTLBase" \cr
  .. .. .. ..$ Title          : chr "gQTLBase: infrastructure for eQTL, mQTL and similar studies" \cr
  .. .. .. ..$ Version        : chr "0.0.7" \cr
  .. .. .. ..$ Author         : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ Description    : chr "infrastructure for eQTL, mQTL and similar studies" \cr
  .. .. .. ..$ Suggests       : chr "geuvStore, knitr, rmarkdown, BiocStyle" \cr
  .. .. .. ..$ Imports        : chr "GenomicRanges, methods, BatchJobs, BiocParallel, BBmisc,S4Vectors, ff, ffbase, BiocGenerics" \cr
  .. .. .. ..$ Depends        : chr "" \cr
  .. .. .. ..$ Maintainer     : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad       : chr "yes" \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Packaged       : chr "2014-11-23 12:40:41 UTC; stvjc" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-25 15:31:33 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/gQTLBase/Meta/package.rds" \cr
  .. .. ..$ geuvStore    :List of 11 \cr
  .. .. .. ..$ Package    : chr "geuvStore" \cr
  .. .. .. ..$ Title      : chr "demonstrate storage discipline for eQTL enumerations" \cr
  .. .. .. ..$ Version    : chr "0.2.1" \cr
  .. .. .. ..$ Author     : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ Description: chr "demonstrate storage discipline for eQTL enumerations and analyses\tbased on a selection of GEUVADIS results" \cr
  .. .. .. ..$ Suggests   : chr "Homo.sapiens" \cr
  .. .. .. ..$ Depends    : chr "BatchJobs, GenomicRanges" \cr
  .. .. .. ..$ Maintainer : chr "VJ Carey <stvjc@channing.harvard.edu>" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; ; 2014-11-21 17:45:10 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/geuvStore/Meta/package.rds" \cr
  .. .. ..$ GenomicRanges:List of 15 \cr
  .. .. .. ..$ Package    : chr "GenomicRanges" \cr
  .. .. .. ..$ Title      : chr "Representation and manipulation of genomic intervals" \cr
  .. .. .. ..$ Description: chr "The ability to efficiently represent and manipulate genomic\tannotations and alignments is playing a central role when it com"| __truncated__ \cr
  .. .. .. ..$ Version    : chr "1.19.14" \cr
  .. .. .. ..$ Author     : chr "P. Aboyoun, H. Pages and M. Lawrence" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews  : chr "Genetics, Infrastructure, Sequencing, Annotation, Coverage,GenomeAnnotation" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.10), methods, BiocGenerics (>= 0.11.3), S4Vectors (>=0.5.10), IRanges (>= 2.1.19), GenomeInfoDb (>= 1.1.20)" \cr
  .. .. .. ..$ Imports    : chr "utils, stats, XVector" \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges" \cr
  .. .. .. ..$ Suggests   : chr "AnnotationDbi (>= 1.21.1), AnnotationHub, BSgenome,BSgenome.Hsapiens.UCSC.hg19, BSgenome.Scerevisiae.UCSC.sacCer2,Biostring"| __truncated__ \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "utils.R phicoef.R transcript-utils.R constraint.Rstrand-utils.R range-squeezers.R GenomicRanges-class.RGRanges-class.R Dele"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-22 03:21:36 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:10:22 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/GenomicRanges/Meta/package.rds" \cr
  .. .. ..$ GenomeInfoDb :List of 16 \cr
  .. .. .. ..$ Package        : chr "GenomeInfoDb" \cr
  .. .. .. ..$ Title          : chr "Utilities for manipulating chromosome and other 'seqname'identifiers" \cr
  .. .. .. ..$ Description    : chr "Contains data and functions that\tdefine and allow translation between different chromosome\tsequence naming conventions (e"| __truncated__ \cr
  .. .. .. ..$ Version        : chr "1.3.7" \cr
  .. .. .. ..$ Author         : chr "Sonali Arora, Martin Morgan, Marc Carlson, H. Pages" \cr
  .. .. .. ..$ Maintainer     : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews      : chr "Genetics, DataRepresentation, Annotation, GenomeAnnotation" \cr
  .. .. .. ..$ Depends        : chr "R (>= 3.1), methods, stats4, BiocGenerics, S4Vectors (>=0.2.0), IRanges (>= 1.99.26)" \cr
  .. .. .. ..$ Imports        : chr "methods, BiocGenerics, S4Vectors" \cr
  .. .. .. ..$ Suggests       : chr "GenomicRanges, Rsamtools, GenomicAlignments, BSgenome,GenomicFeatures, BSgenome.Scerevisiae.UCSC.sacCer2,BSgenome.Celegans."| __truncated__ \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate        : chr "utils.R assembly-utils.R fetchExtendedChromInfoFromUCSC.RrankSeqlevels.R seqinfo.R seqlevelsStyle.R seqlevels-wrappers.RSeq"| __truncated__ \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Video          : chr "http://youtu.be/wdEjCYSXa7w" \cr
  .. .. .. ..$ Packaged       : chr "2014-11-15 05:34:47 UTC; biocbuild" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-16 03:48:46 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/GenomeInfoDb/Meta/package.rds" \cr
  .. .. ..$ IRanges      :List of 16 \cr
  .. .. .. ..$ Package      : chr "IRanges" \cr
  .. .. .. ..$ Title        : chr "Infrastructure for manipulating intervals on sequences" \cr
  .. .. .. ..$ Description  : chr "The package provides efficient low-level and highly reusable S4\tclasses for storing ranges of integers, RLE vectors (Run-Len"| __truncated__ \cr
  .. .. .. ..$ Version      : chr "2.1.19" \cr
  .. .. .. ..$ Author       : chr "H. Pages, P. Aboyoun and M. Lawrence" \cr
  .. .. .. ..$ Maintainer   : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews    : chr "Infrastructure, DataRepresentation" \cr
  .. .. .. ..$ Depends      : chr "R (>= 3.1.0), methods, utils, stats, BiocGenerics (>= 0.11.3),S4Vectors (>= 0.5.10)" \cr
  .. .. .. ..$ Imports      : chr "stats4" \cr
  .. .. .. ..$ LinkingTo    : chr "S4Vectors" \cr
  .. .. .. ..$ Suggests     : chr "XVector, GenomicRanges, BSgenome.Celegans.UCSC.ce2, RUnit" \cr
  .. .. .. ..$ License      : chr "Artistic-2.0" \cr
  .. .. .. ..$ ExtraLicenses: chr "The following files in the 'src' directory are licensedfor all use by Jim Kent, in a manner compatible with theArtistic 2.0"| __truncated__ \cr
  .. .. .. ..$ Collate      : chr "subsetting-internals.R Vector-class-leftovers.RHits-class-leftovers.R List-class-leftovers.R List-comparison.RAtomicList-cl"| __truncated__ \cr
  .. .. .. ..$ Packaged     : chr "2014-11-22 02:23:42 UTC; biocbuild" \cr
  .. .. .. ..$ Built        : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:09:59 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/IRanges/Meta/package.rds" \cr
  .. .. ..$ S4Vectors    :List of 14 \cr
  .. .. .. ..$ Package    : chr "S4Vectors" \cr
  .. .. .. ..$ Title      : chr "S4 implementation of vectors and lists" \cr
  .. .. .. ..$ Description: chr "The S4Vectors package defines the Vector and List virtual classes\tand a set of generic functions that extend the semantic of"| __truncated__ \cr
  .. .. .. ..$ Version    : chr "0.5.10" \cr
  .. .. .. ..$ Author     : chr "H. Pages, M. Lawrence and P. Aboyoun" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure, DataRepresentation" \cr
  .. .. .. ..$ Depends    : chr "R (>= 3.1.0), methods, utils, stats, stats4, BiocGenerics (>=0.11.3)" \cr
  .. .. .. ..$ Imports    : chr "methods, utils, stats, stats4, BiocGenerics" \cr
  .. .. .. ..$ Suggests   : chr "IRanges, RUnit" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "S4-utils.R utils.R normarg-utils.R vector-utils.R isSorted.Rlogical-utils.R int-utils.R str-utils.R eval-utils.RDataTable-c"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-21 05:53:32 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:06:19 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/S4Vectors/Meta/package.rds" \cr
  .. .. ..$ BiocGenerics :List of 14 \cr
  .. .. .. ..$ Package    : chr "BiocGenerics" \cr
  .. .. .. ..$ Title      : chr "S4 generic functions for Bioconductor" \cr
  .. .. .. ..$ Description: chr "S4 generic functions needed by many Bioconductor packages." \cr
  .. .. .. ..$ Version    : chr "0.13.2" \cr
  .. .. .. ..$ Author     : chr "The Bioconductor Dev Team" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure" \cr
  .. .. .. ..$ Depends    : chr "methods, utils, graphics, stats, parallel" \cr
  .. .. .. ..$ Imports    : chr "methods, utils, graphics, stats, parallel" \cr
  .. .. .. ..$ Suggests   : chr "Biobase, S4Vectors, IRanges, GenomicRanges, AnnotationDbi,oligoClasses, oligo, affyPLM, flowClust, affy, RUnit, DESeq2" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "S3-classes-as-S4-classes.R append.R as.data.frame.R as.list.Ras.vector.R cbind.R do.call.R duplicated.R eval.R Extremes.Rfu"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-18 04:14:48 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; ; 2014-11-19 00:52:20 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BiocGenerics/Meta/package.rds" \cr
  .. .. ..$ BatchJobs    :List of 20 \cr
  .. .. .. ..$ Package         : chr "BatchJobs" \cr
  .. .. .. ..$ Title           : chr "Batch computing with R." \cr
  .. .. .. ..$ Description     : chr "Provides Map, Reduce and Filter variants to generate jobs on batch    computing systems like PBS/Torque, LSF, SLURM and Sun G"| __truncated__ \cr
  .. .. .. ..$ Author          : chr "Bernd Bischl <bernd_bischl@gmx.net>,    Michel Lang <michellang@gmail.com>" \cr
  .. .. .. ..$ Maintainer      : chr "Bernd Bischl <bernd_bischl@gmx.net>" \cr
  .. .. .. ..$ URL             : chr "https://github.com/tudo-r/BatchJobs" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/tudo-r/BatchJobs/issues" \cr
  .. .. .. ..$ MailingList     : chr "batchjobs@googlegroups.com" \cr
  .. .. .. ..$ License         : chr "BSD_2_clause + file LICENSE" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.15.0), BBmisc (>= 1.8), methods" \cr
  .. .. .. ..$ Imports         : chr "brew, checkmate (>= 1.5), DBI, digest, fail (>= 1.2), RSQLite(>= 1.0.0), sendmailR, stats, stringr, utils" \cr
  .. .. .. ..$ Suggests        : chr "MASS, testthat" \cr
  .. .. .. ..$ LazyData        : chr "yes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Version         : chr "1.5" \cr
  .. .. .. ..$ Packaged        : chr "2014-10-30 18:49:07 UTC; bischl" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-10-30 20:13:25" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-11-03 09:48:06 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BatchJobs/Meta/package.rds" \cr
  .. .. ..$ BBmisc       :List of 20 \cr
  .. .. .. ..$ Package         : chr "BBmisc" \cr
  .. .. .. ..$ Title           : chr "Miscellaneous helper functions for B. Bischl." \cr
  .. .. .. ..$ Description     : chr "Miscellaneous helper functions for and from B. Bischl and    some other guys at TU Dortmund, mainly for package development." \cr
  .. .. .. ..$ Authors@R       : chr "c(person(\"Bernd\", \"Bischl\", email = \"bernd_bischl@gmx.net\",    role = c(\"aut\", \"cre\")), person(\"Michel\", \"Lang\""| __truncated__ \cr
  .. .. .. ..$ URL             : chr "https://github.com/berndbischl/BBmisc" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/berndbischl/BBmisc/issues" \cr
  .. .. .. ..$ License         : chr "BSD_3_clause + file LICENSE" \cr
  .. .. .. ..$ Encoding        : chr "UTF-8" \cr
  .. .. .. ..$ Imports         : chr "stats, checkmate" \cr
  .. .. .. ..$ Suggests        : chr "testthat, microbenchmark" \cr
  .. .. .. ..$ LazyData        : chr "yes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Version         : chr "1.8" \cr
  .. .. .. ..$ Packaged        : chr "2014-10-30 17:28:43 UTC; bischl" \cr
  .. .. .. ..$ Author          : chr "Bernd Bischl [aut, cre],  Michel Lang [aut],  Jakob Bossek [aut],  Daniel Horn [aut],  Jakob Richter [aut]" \cr
  .. .. .. ..$ Maintainer      : chr "Bernd Bischl <bernd_bischl@gmx.net>" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-10-30 19:34:48" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 09:47:48 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BBmisc/Meta/package.rds" \cr
  .. .. ..$ BiocStyle    :List of 11 \cr
  .. .. .. ..$ Package        : chr "BiocStyle" \cr
  .. .. .. ..$ Title          : chr "Standard styles for vignettes and other Bioconductor documents" \cr
  .. .. .. ..$ Description    : chr "Provides standard formatting styles for Bioconductor PDF    and HTML documents. Package vignettes illustrate use and    fun"| __truncated__ \cr
  .. .. .. ..$ Version        : chr "1.5.3" \cr
  .. .. .. ..$ Author         : chr "Martin Morgan, Andrzej Oles, Wolfgang Huber" \cr
  .. .. .. ..$ Maintainer     : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Suggests       : chr "knitr (>= 1.7), rmarkdown, BiocGenerics, RUnit" \cr
  .. .. .. ..$ biocViews      : chr "Software" \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-25 15:28:08 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BiocStyle/Meta/package.rds" \cr
  .. .. ..$ rmarkdown    :List of 19 \cr
  .. .. .. ..$ Package           : chr "rmarkdown" \cr
  .. .. .. ..$ Type              : chr "Package" \cr
  .. .. .. ..$ Title             : chr "Dynamic Documents for R" \cr
  .. .. .. ..$ Version           : chr "0.3.3" \cr
  .. .. .. ..$ Date              : chr "2014-09-08" \cr
  .. .. .. ..$ Author            : chr "JJ Allaire, Jonathan McPherson, Yihui Xie, Hadley Wickham, Joe Cheng,    Jeff Allen" \cr
  .. .. .. ..$ Maintainer        : chr "JJ Allaire <jj@rstudio.com>" \cr
  .. .. .. ..$ Description       : chr "Convert R Markdown documents into a variety of formats    including HTML, MS Word, PDF, and Beamer." \cr
  .. .. .. ..$ Depends           : chr "R (>= 3.0)" \cr
  .. .. .. ..$ Imports           : chr "tools, utils, knitr (>= 1.6), yaml (>= 2.1.5), htmltools (>=0.2.4), caTools" \cr
  .. .. .. ..$ Suggests          : chr "shiny (>= 0.10.1), testthat, digest" \cr
  .. .. .. ..$ SystemRequirements: chr "pandoc (>= 1.12.3) -http://johnmacfarlane.net/pandoc" \cr
  .. .. .. ..$ URL               : chr "http://rmarkdown.rstudio.com" \cr
  .. .. .. ..$ License           : chr "GPL-3" \cr
  .. .. .. ..$ Packaged          : chr "2014-09-17 11:35:12 UTC; jjallaire" \cr
  .. .. .. ..$ NeedsCompilation  : chr "no" \cr
  .. .. .. ..$ Repository        : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication  : chr "2014-09-17 18:12:29" \cr
  .. .. .. ..$ Built             : chr "R 3.2.0; ; 2014-10-14 13:44:44 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/rmarkdown/Meta/package.rds" \cr
  .. .. ..$ knitr        :List of 22 \cr
  .. .. .. ..$ Package         : chr "knitr" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "A General-Purpose Package for Dynamic Report Generation in R" \cr
  .. .. .. ..$ Version         : chr "1.8" \cr
  .. .. .. ..$ Date            : chr "2014-11-11" \cr
  .. .. .. ..$ Authors@R       : chr "as.person(c(    \"Adam Vogt [ctb]\",    \"Alastair Andrew [ctb]\",    \"Alex Zvoleff [ctb]\",    \"Ashley Manton [ctb]\"| __truncated__ \cr
  .. .. .. ..$ Maintainer      : chr "Yihui Xie <xie@yihui.name>" \cr
  .. .. .. ..$ Description     : chr "This package provides a general-purpose tool for dynamic report    generation in R, which can be used to deal with any type o"| __truncated__ \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.14.1)" \cr
  .. .. .. ..$ Imports         : chr "evaluate (>= 0.5.5), digest, formatR (>= 1.0), highr (>= 0.4),markdown (>= 0.7.4), stringr (>= 0.6), tools" \cr
  .. .. .. ..$ Suggests        : chr "testit, rgl, codetools, rmarkdown, XML, RCurl" \cr
  .. .. .. ..$ License         : chr "GPL" \cr
  .. .. .. ..$ URL             : chr "http://yihui.name/knitr/" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/yihui/knitr/issues" \cr
  .. .. .. ..$ VignetteBuilder : chr "knitr" \cr
  .. .. .. ..$ Collate         : chr "'block.R' 'cache.R' 'citation.R' 'utils.R' 'plot.R''defaults.R' 'concordance.R' 'engine.R' 'themes.R''highlight.R' 'header."| __truncated__ \cr
  .. .. .. ..$ Packaged        : chr "2014-11-11 04:22:18 UTC; yihui" \cr
  .. .. .. ..$ Author          : chr "Adam Vogt [ctb],  Alastair Andrew [ctb],  Alex Zvoleff [ctb],  Ashley Manton [ctb],  Brian Diggs [ctb],  Cassio Perei"| __truncated__ \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-11-11 09:04:50" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-11-16 03:47:54 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/knitr/Meta/package.rds" \cr
  .. .. ..$ weaver       :List of 13 \cr
  .. .. .. ..$ Package    : chr "weaver" \cr
  .. .. .. ..$ Title      : chr "Tools and extensions for processing Sweave documents" \cr
  .. .. .. ..$ Version    : chr "1.33.0" \cr
  .. .. .. ..$ Author     : chr "Seth Falcon" \cr
  .. .. .. ..$ Description: chr "This package provides enhancements on the Sweave() function  in the base package.  In particular a facility for caching code "| __truncated__ \cr
  .. .. .. ..$ Maintainer : chr "Seth Falcon <seth@userprimary.net>" \cr
  .. .. .. ..$ License    : chr "GPL-2" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.5.0), digest, tools, utils, codetools" \cr
  .. .. .. ..$ Suggests   : chr "codetools" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ Packaged   : chr "2014-10-14 02:47:30 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; ; 2014-10-16 00:54:03 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/weaver/Meta/package.rds" \cr
  .. .. ..$ codetools    :List of 14 \cr
  .. .. .. ..$ Package         : chr "codetools" \cr
  .. .. .. ..$ Version         : chr "0.2-9" \cr
  .. .. .. ..$ Priority        : chr "recommended" \cr
  .. .. .. ..$ Author          : chr "Luke Tierney <luke-tierney@uiowa.edu>" \cr
  .. .. .. ..$ Description     : chr "Code analysis tools for R" \cr
  .. .. .. ..$ Title           : chr "Code Analysis Tools for R" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.1)" \cr
  .. .. .. ..$ Maintainer      : chr "Luke Tierney <luke-tierney@uiowa.edu>" \cr
  .. .. .. ..$ License         : chr "GPL" \cr
  .. .. .. ..$ Packaged        : chr "2014-08-18 19:30:24 UTC; luke" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-08-21 10:48:50" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-10-15 20:02:22 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/codetools/Meta/package.rds" \cr
  .. .. ..$ digest       :List of 15 \cr
  .. .. .. ..$ Package         : chr "digest" \cr
  .. .. .. ..$ Version         : chr "0.6.4" \cr
  .. .. .. ..$ Date            : chr "$Date: 2013-12-02 21:56:30 -0600 (Mon, 02 Dec 2013) $" \cr
  .. .. .. ..$ Author          : chr "Dirk Eddelbuettel <edd@debian.org> with contributions  by Antoine Lucas, Jarek Tuszynski, Henrik Bengtsson, Simon Urbanek, "| __truncated__ \cr
  .. .. .. ..$ Maintainer      : chr "Dirk Eddelbuettel <edd@debian.org>" \cr
  .. .. .. ..$ Title           : chr "Create cryptographic hash digests of R objects" \cr
  .. .. .. ..$ Description     : chr "The digest package provides a function 'digest()' for the  creation of hash digests of arbitrary R objects (using the md5, sh"| __truncated__ \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.4.1)" \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ URL             : chr "http://dirk.eddelbuettel.com/code/digest.html" \cr
  .. .. .. ..$ Packaged        : chr "2013-12-03 04:02:00.357418 UTC; edd" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2013-12-03 07:39:12" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:42:11 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/digest/Meta/package.rds" \cr
  .. .. ..$ BiocInstaller:List of 12 \cr
  .. .. .. ..$ Package    : chr "BiocInstaller" \cr
  .. .. .. ..$ Title      : chr "Install/Update Bioconductor and CRAN Packages" \cr
  .. .. .. ..$ Description: chr "Installs/updates Bioconductor and CRAN packages" \cr
  .. .. .. ..$ Version    : chr "1.17.1" \cr
  .. .. .. ..$ Author     : chr "Dan Tenenbaum and Biocore Team" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews  : chr "Software" \cr
  .. .. .. ..$ Depends    : chr "R (>= 3.2.0)" \cr
  .. .. .. ..$ Suggests   : chr "RUnit, BiocGenerics" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Packaged   : chr "2014-10-15 20:33:04 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; ; 2014-10-15 21:00:06 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BiocInstaller/Meta/package.rds" \cr
  .. ..$ loadedOnly:List of 40 \cr
  .. .. ..$ AnnotationDbi    :List of 16 \cr
  .. .. .. ..$ Package        : chr "AnnotationDbi" \cr
  .. .. .. ..$ Title          : chr "Annotation Database Interface" \cr
  .. .. .. ..$ Description    : chr "Provides user interface and database connection code    for annotation data packages using SQLite data storage." \cr
  .. .. .. ..$ Version        : chr "1.29.10" \cr
  .. .. .. ..$ Author         : chr "Herve Pages, Marc Carlson, Seth Falcon, Nianhua Li" \cr
  .. .. .. ..$ Maintainer     : chr "Bioconductor Package Maintainer<maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Depends        : chr "R (>= 2.7.0), methods, utils, stats4, BiocGenerics (>=0.11.2), Biobase (>= 1.17.0), GenomeInfoDb(>= 0.99.17)" \cr
  .. .. .. ..$ Imports        : chr "methods, utils, DBI, RSQLite, stats4, BiocGenerics, Biobase,S4Vectors" \cr
  .. .. .. ..$ Suggests       : chr "DBI (>= 0.2-4), RSQLite (>= 0.6-4), hgu95av2.db, GO.db,org.Sc.sgd.db, org.At.tair.db, KEGG.db, RUnit,TxDb.Hsapiens.UCSC.hg1"| __truncated__ \cr
  .. .. .. ..$ Collate        : chr "00RTobjs.R unlist2.R AllGenerics.R AllClasses.R FlatBimap.RSQL.R AnnDbObj-lowAPI.R Bimap.R GOTerms.R BimapFormatting.RBimap"| __truncated__ \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ biocViews      : chr "Annotation, Microarray, Sequencing, GenomeAnnotation" \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Packaged       : chr "2014-11-22 02:40:14 UTC; biocbuild" \cr
  .. .. .. ..$ Video          : chr "https://www.youtube.com/watch?v=8qvGNTVz3Ik" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-23 12:09:48 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/AnnotationDbi/Meta/package.rds" \cr
  .. .. ..$ base64enc        :List of 15 \cr
  .. .. .. ..$ Package         : chr "base64enc" \cr
  .. .. .. ..$ Version         : chr "0.1-2" \cr
  .. .. .. ..$ Title           : chr "Tools for base64 encoding" \cr
  .. .. .. ..$ Author          : chr "Simon Urbanek <Simon.Urbanek@r-project.org>" \cr
  .. .. .. ..$ Maintainer      : chr "Simon Urbanek <Simon.Urbanek@r-project.org>" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.9.0)" \cr
  .. .. .. ..$ Enhances        : chr "png" \cr
  .. .. .. ..$ Description     : chr "This package provides tools for handling base64 encoding. It is more flexible than the orphaned base64 package." \cr
  .. .. .. ..$ License         : chr "GPL-2 | GPL-3" \cr
  .. .. .. ..$ URL             : chr "http://www.rforge.net/base64enc" \cr
  .. .. .. ..$ Packaged        : chr "2014-06-26 14:47:44 UTC; svnuser" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-06-26 17:15:25" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:45:21 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/base64enc/Meta/package.rds" \cr
  .. .. ..$ Biobase          :List of 15 \cr
  .. .. .. ..$ Package    : chr "Biobase" \cr
  .. .. .. ..$ Title      : chr "Biobase: Base functions for Bioconductor" \cr
  .. .. .. ..$ Version    : chr "2.27.0" \cr
  .. .. .. ..$ Author     : chr "R. Gentleman, V. Carey, M. Morgan, S. Falcon" \cr
  .. .. .. ..$ Description: chr "Functions that are needed by many other packages or which        replace R functions." \cr
  .. .. .. ..$ Suggests   : chr "tools, tkWidgets, ALL, RUnit, golubEsets" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.10), BiocGenerics (>= 0.3.2), utils" \cr
  .. .. .. ..$ Imports    : chr "methods" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer<maintainer@bioconductor.org>" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "tools.R strings.R environment.R vignettes.R packages.RAllGenerics.R VersionsClass.R VersionedClasses.Rmethods-VersionsNull."| __truncated__ \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure" \cr
  .. .. .. ..$ Packaged   : chr "2014-10-14 01:23:01 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-10-16 00:54:06 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/Biobase/Meta/package.rds" \cr
  .. .. ..$ BiocParallel     :List of 18 \cr
  .. .. .. ..$ Package        : chr "BiocParallel" \cr
  .. .. .. ..$ Type           : chr "Package" \cr
  .. .. .. ..$ Title          : chr "Bioconductor facilities for parallel evaluation" \cr
  .. .. .. ..$ Version        : chr "1.1.6" \cr
  .. .. .. ..$ Authors@R      : chr "c(    person(\"Bioconductor Package Maintainer\",        email=\"maintainer@bioconductor.org\", role=\"cre\"),    person("| __truncated__ \cr
  .. .. .. ..$ Description    : chr "This package provides modified versions and novel    implementation of functions for parallel evaluation, tailored to    us"| __truncated__ \cr
  .. .. .. ..$ biocViews      : chr "Infrastructure" \cr
  .. .. .. ..$ License        : chr "GPL-2 | GPL-3" \cr
  .. .. .. ..$ Depends        : chr "methods" \cr
  .. .. .. ..$ Imports        : chr "parallel, foreach, tools, BatchJobs, BBmisc, BiocGenerics" \cr
  .. .. .. ..$ Suggests       : chr "doParallel, snow, Rmpi, GenomicRanges,RNAseqData.HNRNPC.bam.chr14, Rsamtools, GenomicAlignments,ShortRead, RUnit, BiocStyle"| __truncated__ \cr
  .. .. .. ..$ Collate.unix   : chr "AllGenerics.R BiocParallelParam-class.R ErrorHandling.Rbpbackend-methods.R bpisup-methods.R bplapply-methods.Rbpmapply-meth"| __truncated__ \cr
  .. .. .. ..$ Collate.windows: chr "AllGenerics.R BiocParallelParam-class.RErrorHandling.R bpbackend-methods.R bpisup-methods.Rbplapply-methods.R bpmapply-meth"| __truncated__ \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ Packaged       : chr "2014-11-20 04:56:05 UTC; biocbuild" \cr
  .. .. .. ..$ Author         : chr "Bioconductor Package Maintainer [cre],  Martin Morgan [aut],  Michel Lang [aut],  Ryan Thompson [aut]" \cr
  .. .. .. ..$ Maintainer     : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-23 12:05:30 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BiocParallel/Meta/package.rds" \cr
  .. .. ..$ biomaRt          :List of 15 \cr
  .. .. .. ..$ Package     : chr "biomaRt" \cr
  .. .. .. ..$ Version     : chr "2.23.5" \cr
  .. .. .. ..$ Title       : chr "Interface to BioMart databases (e.g. Ensembl, COSMIC ,Wormbaseand Gramene)" \cr
  .. .. .. ..$ Author      : chr "Steffen Durinck <durincks@gene.com>, Wolfgang Huber" \cr
  .. .. .. ..$ Contributors: chr "Sean Davis <sdavis2@mail.nih.gov>, Francois Pepin, VinceS. Buffalo" \cr
  .. .. .. ..$ Maintainer  : chr "Steffen Durinck <durincks@gene.com>" \cr
  .. .. .. ..$ Depends     : chr "methods" \cr
  .. .. .. ..$ Imports     : chr "utils, XML, RCurl, AnnotationDbi" \cr
  .. .. .. ..$ Suggests    : chr "annotate" \cr
  .. .. .. ..$ biocViews   : chr "Annotation" \cr
  .. .. .. ..$ Description : chr "In recent years a wealth of biological data has become available in public data repositories. Easy access to these valuable dat"| __truncated__ \cr
  .. .. .. ..$ License     : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad    : chr "yes" \cr
  .. .. .. ..$ Packaged    : chr "2014-11-22 02:27:50 UTC; biocbuild" \cr
  .. .. .. ..$ Built       : chr "R 3.2.0; ; 2014-11-23 12:10:13 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/biomaRt/Meta/package.rds" \cr
  .. .. ..$ Biostrings       :List of 17 \cr
  .. .. .. ..$ Package    : chr "Biostrings" \cr
  .. .. .. ..$ Title      : chr "String objects representing biological sequences, and matchingalgorithms" \cr
  .. .. .. ..$ Description: chr "Memory efficient string containers, string matching\talgorithms, and other utilities, for fast manipulation of large\tbiolo"| __truncated__ \cr
  .. .. .. ..$ Version    : chr "2.35.4" \cr
  .. .. .. ..$ Author     : chr "H. Pages, P. Aboyoun, R. Gentleman, and S. DebRoy" \cr
  .. .. .. ..$ Maintainer : chr "H. Pages <hpages@fhcrc.org>" \cr
  .. .. .. ..$ biocViews  : chr "SequenceMatching, Genetics, Sequencing, Infrastructure,DataImport, DataRepresentation" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.8.0), methods, BiocGenerics (>= 0.11.3), S4Vectors (>=0.2.2), IRanges (>= 2.1.2), XVector (>= 0.7.1)" \cr
  .. .. .. ..$ Imports    : chr "graphics, methods, stats, utils, BiocGenerics, IRanges,XVector, zlibbioc" \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges, XVector" \cr
  .. .. .. ..$ Enhances   : chr "Rmpi" \cr
  .. .. .. ..$ Suggests   : chr "BSgenome (>= 1.13.14), BSgenome.Celegans.UCSC.ce2 (>=1.3.11), BSgenome.Dmelanogaster.UCSC.dm3 (>= 1.3.11),BSgenome.Hsapiens"| __truncated__ \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ Collate    : chr "00datacache.R utils.R IUPAC_CODE_MAP.R AMINO_ACID_CODE.RGENETIC_CODE.R XStringCodec-class.R seqtype.R XString-class.RXStrin"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-12 02:28:16 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-16 03:48:31 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/Biostrings/Meta/package.rds" \cr
  .. .. ..$ bit              :List of 19 \cr
  .. .. .. ..$ Package         : chr "bit" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "A class for vectors of 1-bit booleans" \cr
  .. .. .. ..$ Version         : chr "1.1-12" \cr
  .. .. .. ..$ Date            : chr "2012-01-17" \cr
  .. .. .. ..$ Author          : chr "Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>" \cr
  .. .. .. ..$ Maintainer      : chr "Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.9.2)" \cr
  .. .. .. ..$ Description     : chr "bitmapped vectors of booleans (no NAs),   coercion from and to logicals, integers and integer subscripts;   fast boolean op"| __truncated__ \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ LazyLoad        : chr "yes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Encoding        : chr "latin1" \cr
  .. .. .. ..$ URL             : chr "http://ff.r-forge.r-project.org/" \cr
  .. .. .. ..$ Packaged        : chr "2014-04-08 21:51:43 UTC; jens" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-04-09 09:54:10" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:46:11 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/bit/Meta/package.rds" \cr
  .. .. ..$ bitops           :List of 13 \cr
  .. .. .. ..$ Package         : chr "bitops" \cr
  .. .. .. ..$ Version         : chr "1.0-6" \cr
  .. .. .. ..$ Date            : chr "2013-08-17" \cr
  .. .. .. ..$ Author          : chr "S original by Steve Dutky <sdutky@terpalum.umd.edu> initial R        port and extensions by Martin Maechler; revised and modi"| __truncated__ \cr
  .. .. .. ..$ Maintainer      : chr "Martin Maechler <maechler@stat.math.ethz.ch>" \cr
  .. .. .. ..$ Title           : chr "Bitwise Operations" \cr
  .. .. .. ..$ Description     : chr "Functions for bitwise operations on integer vectors." \cr
  .. .. .. ..$ License         : chr "GPL (>= 2)" \cr
  .. .. .. ..$ Packaged        : chr "2013-08-17 15:58:57 UTC; maechler" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2013-08-17 21:10:34" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:45:57 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/bitops/Meta/package.rds" \cr
  .. .. ..$ brew             :List of 13 \cr
  .. .. .. ..$ Package         : chr "brew" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "Templating Framework for Report Generation" \cr
  .. .. .. ..$ Version         : chr "1.0-6" \cr
  .. .. .. ..$ Date            : chr "2010-09-30" \cr
  .. .. .. ..$ Author          : chr "Jeffrey Horner" \cr
  .. .. .. ..$ Maintainer      : chr "Jeffrey Horner <jeffrey.horner@gmail.com>" \cr
  .. .. .. ..$ Description     : chr "brew implements a templating framework for mixing text and        R code for report generation. brew template syntax is simil"| __truncated__ \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ Packaged        : chr "2011-04-13 14:36:39 UTC; hornerj" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2011-04-13 15:16:08" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-08-13 10:45:25 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/brew/Meta/package.rds" \cr
  .. .. ..$ BSgenome         :List of 15 \cr
  .. .. .. ..$ Package    : chr "BSgenome" \cr
  .. .. .. ..$ Title      : chr "Infrastructure for Biostrings-based genome data packages" \cr
  .. .. .. ..$ Description: chr "Infrastructure shared by all the Biostrings-based genome data             packages " \cr
  .. .. .. ..$ Version    : chr "1.35.8" \cr
  .. .. .. ..$ Author     : chr "Herve Pages" \cr
  .. .. .. ..$ Maintainer : chr "H. Pages <hpages@fhcrc.org>" \cr
  .. .. .. ..$ biocViews  : chr "Genetics, Infrastructure, DataRepresentation,SequenceMatching, Annotation, SNP" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.8.0), methods, BiocGenerics (>= 0.1.2), S4Vectors (>=0.5.10), IRanges (>= 2.1.3), GenomeInfoDb (>= 1.1.4),GenomicRa"| __truncated__ \cr
  .. .. .. ..$ Imports    : chr "methods, stats, BiocGenerics, S4Vectors, IRanges, XVector,GenomeInfoDb, GenomicRanges, Biostrings, Rsamtools, rtracklayer" \cr
  .. .. .. ..$ Suggests   : chr "BiocInstaller, BSgenome.Celegans.UCSC.ce2 (>= 1.3.11),BSgenome.Hsapiens.UCSC.hg19 (>= 1.3.11),BSgenome.Hsapiens.UCSC.hg19.m"| __truncated__ \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ Collate    : chr "utils.R OnDiskLongTable-class.R OnDiskNamedSequences-class.RSNPlocs-class.R InjectSNPsHandler-class.R XtraSNPlocs-class.RBS"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-21 02:41:48 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; ; 2014-11-23 12:12:06 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/BSgenome/Meta/package.rds" \cr
  .. .. ..$ checkmate        :List of 19 \cr
  .. .. .. ..$ Package         : chr "checkmate" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "Fast and versatile argument checks" \cr
  .. .. .. ..$ Description     : chr "Tests and assertions to perform frequent argument checks. A    substantial part of the package was written in C to minimize a"| __truncated__ \cr
  .. .. .. ..$ Version         : chr "1.5.0" \cr
  .. .. .. ..$ Author          : chr "Michel Lang <michellang@gmail.com>,    Bernd Bischl <bernd_bischl@gmx.net>" \cr
  .. .. .. ..$ Maintainer      : chr "Michel Lang <michellang@gmail.com>" \cr
  .. .. .. ..$ URL             : chr "https://github.com/mllg/checkmate" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/mllg/checkmate/issues" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Encoding        : chr "UTF-8" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.15.0)" \cr
  .. .. .. ..$ Suggests        : chr "testthat" \cr
  .. .. .. ..$ License         : chr "BSD_3_clause + file LICENSE" \cr
  .. .. .. ..$ Packaged        : chr "2014-10-18 20:58:35 UTC; lang" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-10-19 01:41:06" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 09:46:42 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/checkmate/Meta/package.rds" \cr
  .. .. ..$ DBI              :List of 17 \cr
  .. .. .. ..$ Package         : chr "DBI" \cr
  .. .. .. ..$ Version         : chr "0.3.1" \cr
  .. .. .. ..$ Title           : chr "R Database Interface" \cr
  .. .. .. ..$ Author          : chr "R Special Interest Group on Databases (R-SIG-DB)" \cr
  .. .. .. ..$ Maintainer      : chr "Hadley Wickham <hadley@rstudio.com>" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.15.0), methods" \cr
  .. .. .. ..$ Suggests        : chr "testthat, RSQLite" \cr
  .. .. .. ..$ Description     : chr "A database interface (DBI) definition for communication    between R and relational database management systems.  All    cl"| __truncated__ \cr
  .. .. .. ..$ License         : chr "LGPL (>= 2)" \cr
  .. .. .. ..$ URL             : chr "https://github.com/rstats-db/DBI" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/rstats-db/DBI/issues" \cr
  .. .. .. ..$ Collate         : chr "'DBObject.R' 'DBConnection.R' 'DBDriver.R' 'DBResult.R''compliance.R' 'keywords.R' 'quote.R' 'util.R'" \cr
  .. .. .. ..$ Packaged        : chr "2014-09-23 21:42:02 UTC; hadley" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-09-24 07:27:12" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-10-14 13:46:57 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/DBI/Meta/package.rds" \cr
  .. .. ..$ evaluate         :List of 18 \cr
  .. .. .. ..$ Package         : chr "evaluate" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "Parsing and evaluation tools that provide more details than thedefault." \cr
  .. .. .. ..$ Version         : chr "0.5.5" \cr
  .. .. .. ..$ Date            : chr "2014-04-30" \cr
  .. .. .. ..$ Authors@R       : chr "c(person(\"Hadley\", \"Wickham\", role = \"aut\"), person(\"Yihui\", \"Xie\",    role = c(\"cre\", \"ctb\"), email = \"xie@yi"| __truncated__ \cr
  .. .. .. ..$ Description     : chr "Parsing and evaluation tools that make it easy to recreate the    command line behaviour of R." \cr
  .. .. .. ..$ License         : chr "GPL" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.14.0)" \cr
  .. .. .. ..$ Imports         : chr "stringr (>= 0.6.2)" \cr
  .. .. .. ..$ Suggests        : chr "testthat, lattice, ggplot2" \cr
  .. .. .. ..$ Packaged        : chr "2014-04-29 19:01:11 UTC; yihui" \cr
  .. .. .. ..$ Author          : chr "Hadley Wickham [aut],  Yihui Xie [cre, ctb],  Barret Schloerke [ctb]" \cr
  .. .. .. ..$ Maintainer      : chr "Yihui Xie <xie@yihui.name>" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-04-29 22:08:31" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-08-13 14:37:13 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/evaluate/Meta/package.rds" \cr
  .. .. ..$ fail             :List of 18 \cr
  .. .. .. ..$ Package         : chr "fail" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "File Abstraction Interface Layer (FAIL) mimicking a key-valuestore" \cr
  .. .. .. ..$ Description     : chr "More comfortable interface to work with a directory of R data or    source files" \cr
  .. .. .. ..$ Version         : chr "1.2" \cr
  .. .. .. ..$ Author          : chr "Michel Lang <michellang@gmail.com>" \cr
  .. .. .. ..$ Maintainer      : chr "Michel Lang <michellang@gmail.com>" \cr
  .. .. .. ..$ URL             : chr "https://github.com/mllg/fail" \cr
  .. .. .. ..$ License         : chr "BSD_3_clause + file LICENSE" \cr
  .. .. .. ..$ Imports         : chr "BBmisc" \cr
  .. .. .. ..$ Suggests        : chr "testthat" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Collate         : chr "'helper.R' 'methods.R' 'zzz.R' 'cache.R' 'fail.R' 'worker.R''sail.R' 'backends.R' 'constructor.R'" \cr
  .. .. .. ..$ Packaged        : chr "2013-09-19 11:29:30 UTC; lang" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2013-09-19 14:12:29" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-08-13 10:47:54 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/fail/Meta/package.rds" \cr
  .. .. ..$ fastmatch        :List of 12 \cr
  .. .. .. ..$ Package         : chr "fastmatch" \cr
  .. .. .. ..$ Version         : chr "1.0-4" \cr
  .. .. .. ..$ Title           : chr "Fast match() function" \cr
  .. .. .. ..$ Author          : chr "Simon Urbanek <simon.urbanek@r-project.org>" \cr
  .. .. .. ..$ Maintainer      : chr "Simon Urbanek <simon.urbanek@r-project.org>" \cr
  .. .. .. ..$ Description     : chr "Package providing a fast match() replacement for cases        that require repeated look-ups. It is slightly faster that R's\"| __truncated__ \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ URL             : chr "http://www.rforge.net/fastmatch" \cr
  .. .. .. ..$ Packaged        : chr "2012-01-21 10:09:18 UTC; svnuser" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2012-01-21 10:22:24" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 19:06:30 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/fastmatch/Meta/package.rds" \cr
  .. .. ..$ ff               :List of 19 \cr
  .. .. .. ..$ Package         : chr "ff" \cr
  .. .. .. ..$ Version         : chr "2.2-13" \cr
  .. .. .. ..$ Date            : chr "2012-03-29" \cr
  .. .. .. ..$ Title           : chr "memory-efficient storage of large data on disk and fast accessfunctions" \cr
  .. .. .. ..$ Author          : chr "Daniel Adler <dadler@uni-goettingen.de>, \tChristian Gläser <christian_glaeser@gmx.de>,\tOleg Nenadic <onenadi@uni-goetting"| __truncated__ \cr
  .. .. .. ..$ Maintainer      : chr "Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.10.1), bit (>= 1.1-12), utils" \cr
  .. .. .. ..$ Suggests        : chr "biglm" \cr
  .. .. .. ..$ Description     : chr "The ff package provides data structures that are stored on\tdisk but behave (almost) as if they were in RAM by transparently "| __truncated__ \cr
  .. .. .. ..$ License         : chr "GPL-2 | file LICENSE" \cr
  .. .. .. ..$ LazyLoad        : chr "yes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ Encoding        : chr "latin1" \cr
  .. .. .. ..$ URL             : chr "http://ff.r-forge.r-project.org/" \cr
  .. .. .. ..$ Packaged        : chr "2014-04-07 21:16:45 UTC; root" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-04-09 09:54:20" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:47:26 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/ff/Meta/package.rds" \cr
  .. .. ..$ ffbase           :List of 27 \cr
  .. .. .. ..$ Package       : chr "ffbase" \cr
  .. .. .. ..$ Maintainer    : chr "Edwin de Jonge <edwindjonge@gmail.com>" \cr
  .. .. .. ..$ License       : chr "GPL-3" \cr
  .. .. .. ..$ Title         : chr "Basic statistical functions for package ff" \cr
  .. .. .. ..$ Type          : chr "Package" \cr
  .. .. .. ..$ LazyLoad      : chr "yes" \cr
  .. .. .. ..$ Author        : chr "Edwin de Jonge, Jan Wijffels, Jan van der Laan" \cr
  .. .. .. ..$ Description   : chr "Basic (statistical) functionality for package ff" \cr
  .. .. .. ..$ Version       : chr "0.11.3" \cr
  .. .. .. ..$ URL           : chr "http://github.com/edwindj/ffbase" \cr
  .. .. .. ..$ Date          : chr "2013-11-1" \cr
  .. .. .. ..$ Depends       : chr "ff (>= 2.2-11),bit,R (>= 2.12.0)" \cr
  .. .. .. ..$ Imports       : chr "fastmatch" \cr
  .. .. .. ..$ Suggests      : chr "testthat,parallel,LaF,biglm" \cr
  .. .. .. ..$ RemoteType    : chr "github" \cr
  .. .. .. ..$ RemoteHost    : chr "api.github.com" \cr
  .. .. .. ..$ RemoteRepo    : chr "ffbase" \cr
  .. .. .. ..$ RemoteUsername: chr "edwindj" \cr
  .. .. .. ..$ RemoteRef     : chr "master" \cr
  .. .. .. ..$ RemoteSha     : chr "e55a420f1d14fa2674ebcdaae9571861737362b3" \cr
  .. .. .. ..$ RemoteSubdir  : chr "pkg" \cr
  .. .. .. ..$ GithubRepo    : chr "ffbase" \cr
  .. .. .. ..$ GithubUsername: chr "edwindj" \cr
  .. .. .. ..$ GithubRef     : chr "master" \cr
  .. .. .. ..$ GithubSHA1    : chr "e55a420f1d14fa2674ebcdaae9571861737362b3" \cr
  .. .. .. ..$ GithubSubdir  : chr "pkg" \cr
  .. .. .. ..$ Built         : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 19:06:31 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/ffbase/Meta/package.rds" \cr
  .. .. ..$ foreach          :List of 20 \cr
  .. .. .. ..$ Package                         : chr "foreach" \cr
  .. .. .. ..$ Type                            : chr "Package" \cr
  .. .. .. ..$ Title                           : chr "Foreach looping construct for R" \cr
  .. .. .. ..$ Version                         : chr "1.4.2" \cr
  .. .. .. ..$ Author                          : chr "Revolution Analytics, Steve Weston" \cr
  .. .. .. ..$ Maintainer                      : chr "Revolution Analytics <packages@revolutionanalytics.com>" \cr
  .. .. .. ..$ Description                     : chr "Support for the foreach looping construct.  Foreach is an        idiom that allows for iterating over elements in a collectio"| __truncated__ \cr
  .. .. .. ..$ Depends                         : chr "R (>= 2.5.0)" \cr
  .. .. .. ..$ Imports                         : chr "codetools, utils, iterators" \cr
  .. .. .. ..$ Suggests                        : chr "randomForest" \cr
  .. .. .. ..$ Enhances                        : chr "compiler, doMC, RUnit, doParallel" \cr
  .. .. .. ..$ License                         : chr "Apache License (== 2.0)" \cr
  .. .. .. ..$ Repository                      : chr "CRAN" \cr
  .. .. .. ..$ Repository/R-Forge/Project      : chr "foreach" \cr
  .. .. .. ..$ Repository/R-Forge/Revision     : chr "25" \cr
  .. .. .. ..$ Repository/R-Forge/DateTimeStamp: chr "2014-04-10 18:54:16" \cr
  .. .. .. ..$ Date/Publication                : chr "2014-04-11 07:52:42" \cr
  .. .. .. ..$ Packaged                        : chr "2014-04-10 20:16:29 UTC; rforge" \cr
  .. .. .. ..$ NeedsCompilation                : chr "no" \cr
  .. .. .. ..$ Built                           : chr "R 3.2.0; ; 2014-08-13 10:46:17 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/foreach/Meta/package.rds" \cr
  .. .. ..$ formatR          :List of 18 \cr
  .. .. .. ..$ Package         : chr "formatR" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "Format R Code Automatically" \cr
  .. .. .. ..$ Version         : chr "1.0" \cr
  .. .. .. ..$ Date            : chr "2014-08-25" \cr
  .. .. .. ..$ Author          : chr "Yihui Xie" \cr
  .. .. .. ..$ Maintainer      : chr "Yihui Xie <xie@yihui.name>" \cr
  .. .. .. ..$ Description     : chr "This package provides a function tidy_source() to format R source    code. Spaces and indent will be added to the code automa"| __truncated__ \cr
  .. .. .. ..$ Suggests        : chr "codetools, shiny, testit, knitr" \cr
  .. .. .. ..$ License         : chr "GPL" \cr
  .. .. .. ..$ URL             : chr "http://yihui.name/formatR" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/yihui/formatR/issues" \cr
  .. .. .. ..$ VignetteBuilder : chr "knitr" \cr
  .. .. .. ..$ Packaged        : chr "2014-08-24 21:04:31 UTC; yihui" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-08-25 00:55:00" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-10-15 20:02:23 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/formatR/Meta/package.rds" \cr
  .. .. ..$ GenomicAlignments:List of 16 \cr
  .. .. .. ..$ Package    : chr "GenomicAlignments" \cr
  .. .. .. ..$ Title      : chr "Representation and manipulation of short genomic alignments" \cr
  .. .. .. ..$ Description: chr "Provides efficient containers for storing and manipulating\tshort genomic alignments (typically obtained by aligning short re"| __truncated__ \cr
  .. .. .. ..$ Version    : chr "1.3.10" \cr
  .. .. .. ..$ Author     : chr "Herv\\'e Pag\\`es, Valerie Obenchain, Martin Morgan" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ biocViews  : chr "Genetics, Infrastructure, DataImport, Sequencing, RNASeq,SNP, Coverage, Alignment" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.10), methods, BiocGenerics (>= 0.11.3), S4Vectors (>=0.5.10), IRanges (>= 2.1.12), GenomeInfoDb (>= 1.1.20),Genomic"| __truncated__ \cr
  .. .. .. ..$ Imports    : chr "methods, stats, BiocGenerics, S4Vectors, IRanges,GenomicRanges, Biostrings, Rsamtools, BiocParallel" \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges" \cr
  .. .. .. ..$ Suggests   : chr "ShortRead, rtracklayer, BSgenome, GenomicFeatures,RNAseqData.HNRNPC.bam.chr14, pasillaBamSubset,TxDb.Hsapiens.UCSC.hg19.kno"| __truncated__ \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "utils.R cigar-utils.R GAlignments-class.RGAlignmentPairs-class.R GAlignmentsList-class.RGappedReads-class.R OverlapEncoding"| __truncated__ \cr
  .. .. .. ..$ Video      : chr "https://www.youtube.com/watch?v=2KqBSbkfhRo ,https://www.youtube.com/watch?v=3PK_jx44QTs" \cr
  .. .. .. ..$ Packaged   : chr "2014-11-21 05:37:33 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:12:23 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/GenomicAlignments/Meta/package.rds" \cr
  .. .. ..$ GenomicFeatures  :List of 15 \cr
  .. .. .. ..$ Package        : chr "GenomicFeatures" \cr
  .. .. .. ..$ Title          : chr "Tools for making and manipulating transcript centric annotations" \cr
  .. .. .. ..$ Version        : chr "1.19.6" \cr
  .. .. .. ..$ Author         : chr "M. Carlson, H. Pages, P. Aboyoun, S. Falcon, M. Morgan,        D. Sarkar, M. Lawrence" \cr
  .. .. .. ..$ License        : chr "Artistic-2.0" \cr
  .. .. .. ..$ Description    : chr "A set of tools and methods for making and manipulating        transcript centric annotations. With these tools the user can"| __truncated__ \cr
  .. .. .. ..$ Maintainer     : chr "Bioconductor Package Maintainer<maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Depends        : chr "BiocGenerics (>= 0.1.0), S4Vectors (>= 0.1.5), IRanges (>=1.99.1), GenomeInfoDb (>= 1.1.3), GenomicRanges (>= 1.17.12),Anno"| __truncated__ \cr
  .. .. .. ..$ Imports        : chr "methods, DBI (>= 0.2-5), RSQLite (>= 0.8-1), Biostrings (>=2.23.3), rtracklayer (>= 1.25.2), biomaRt (>= 2.17.1), RCurl,uti"| __truncated__ \cr
  .. .. .. ..$ Suggests       : chr "org.Mm.eg.db, BSgenome, BSgenome.Hsapiens.UCSC.hg19 (>=1.3.17), BSgenome.Celegans.UCSC.ce2,BSgenome.Dmelanogaster.UCSC.dm3 "| __truncated__ \cr
  .. .. .. ..$ Collate        : chr "utils.R Ensembl.utils.R findCompatibleMarts.R TxDb-class.RFeatureDb-class.R makeTxDb.R makeTxDbFromUCSC.RmakeTxDbFromBiomar"| __truncated__ \cr
  .. .. .. ..$ VignetteBuilder: chr "knitr" \cr
  .. .. .. ..$ biocViews      : chr "Genetics, Infrastructure, Annotation, Sequencing,GenomeAnnotation" \cr
  .. .. .. ..$ Packaged       : chr "2014-11-04 03:50:57 UTC; biocbuild" \cr
  .. .. .. ..$ Built          : chr "R 3.2.0; ; 2014-11-04 21:46:14 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/GenomicFeatures/Meta/package.rds" \cr
  .. .. ..$ grid             :List of 11 \cr
  .. .. .. ..$ Package    : chr "grid" \cr
  .. .. .. ..$ Version    : chr "3.2.0" \cr
  .. .. .. ..$ Priority   : chr "base" \cr
  .. .. .. ..$ Title      : chr "The Grid Graphics Package" \cr
  .. .. .. ..$ Author     : chr "Paul Murrell <paul@stat.auckland.ac.nz>" \cr
  .. .. .. ..$ Maintainer : chr "R Core Team <R-core@r-project.org>" \cr
  .. .. .. ..$ Description: chr "A rewrite of the graphics layout capabilities, plus some  support for interaction" \cr
  .. .. .. ..$ Imports    : chr "grDevices" \cr
  .. .. .. ..$ Suggests   : chr "lattice" \cr
  .. .. .. ..$ License    : chr "Part of R 3.2.0" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-10-14 13:30:58 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/grid/Meta/package.rds" \cr
  .. .. ..$ htmltools        :List of 18 \cr
  .. .. .. ..$ Package         : chr "htmltools" \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Title           : chr "Tools for HTML" \cr
  .. .. .. ..$ Version         : chr "0.2.6" \cr
  .. .. .. ..$ Date            : chr "2014-07-30" \cr
  .. .. .. ..$ Author          : chr "RStudio, Inc." \cr
  .. .. .. ..$ Maintainer      : chr "Joe Cheng <joe@rstudio.com>" \cr
  .. .. .. ..$ Description     : chr "Tools for HTML generation and output" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.14.1)" \cr
  .. .. .. ..$ Imports         : chr "utils, digest" \cr
  .. .. .. ..$ Suggests        : chr "markdown, testthat" \cr
  .. .. .. ..$ Enhances        : chr "knitr" \cr
  .. .. .. ..$ License         : chr "GPL (>= 2)" \cr
  .. .. .. ..$ Packaged        : chr "2014-09-08 16:45:00 UTC; jcheng" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-09-08 19:23:28" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-10-14 13:44:43 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/htmltools/Meta/package.rds" \cr
  .. .. ..$ iterators        :List of 18 \cr
  .. .. .. ..$ Package                         : chr "iterators" \cr
  .. .. .. ..$ Type                            : chr "Package" \cr
  .. .. .. ..$ Title                           : chr "Iterator construct for R" \cr
  .. .. .. ..$ Version                         : chr "1.0.7" \cr
  .. .. .. ..$ Author                          : chr "Revolution Analytics" \cr
  .. .. .. ..$ Maintainer                      : chr "Revolution Analytics <packages@revolutionanalytics.com>" \cr
  .. .. .. ..$ Description                     : chr "Support for iterators, which allow a programmer to traverse        through all the elements of a vector, list, or other colle"| __truncated__ \cr
  .. .. .. ..$ Depends                         : chr "R (>= 2.5.0), utils" \cr
  .. .. .. ..$ Suggests                        : chr "RUnit" \cr
  .. .. .. ..$ License                         : chr "Apache License (== 2.0)" \cr
  .. .. .. ..$ Repository                      : chr "CRAN" \cr
  .. .. .. ..$ Repository/R-Forge/Project      : chr "foreach" \cr
  .. .. .. ..$ Repository/R-Forge/Revision     : chr "25" \cr
  .. .. .. ..$ Repository/R-Forge/DateTimeStamp: chr "2014-04-10 18:54:16" \cr
  .. .. .. ..$ Date/Publication                : chr "2014-04-11 07:52:44" \cr
  .. .. .. ..$ Packaged                        : chr "2014-04-10 20:15:28 UTC; rforge" \cr
  .. .. .. ..$ NeedsCompilation                : chr "no" \cr
  .. .. .. ..$ Built                           : chr "R 3.2.0; ; 2014-08-13 10:46:16 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/iterators/Meta/package.rds" \cr
  .. .. ..$ lattice          :List of 22 \cr
  .. .. .. ..$ Package         : chr "lattice" \cr
  .. .. .. ..$ Version         : chr "0.20-29" \cr
  .. .. .. ..$ Date            : chr "2014/04/01" \cr
  .. .. .. ..$ Priority        : chr "recommended" \cr
  .. .. .. ..$ Title           : chr "Lattice Graphics" \cr
  .. .. .. ..$ Author          : chr "Deepayan Sarkar <deepayan.sarkar@r-project.org>" \cr
  .. .. .. ..$ Maintainer      : chr "Deepayan Sarkar <deepayan.sarkar@r-project.org>" \cr
  .. .. .. ..$ Description     : chr "Lattice is a powerful and elegant high-level data  visualization system, with an emphasis on multivariate data, that is  su"| __truncated__ \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.15.1)" \cr
  .. .. .. ..$ Suggests        : chr "KernSmooth, MASS" \cr
  .. .. .. ..$ Imports         : chr "grid, grDevices, graphics, stats, utils" \cr
  .. .. .. ..$ Enhances        : chr "chron" \cr
  .. .. .. ..$ LazyLoad        : chr "yes" \cr
  .. .. .. ..$ LazyData        : chr "yes" \cr
  .. .. .. ..$ License         : chr "GPL (>= 2)" \cr
  .. .. .. ..$ URL             : chr "http://lattice.r-forge.r-project.org/" \cr
  .. .. .. ..$ BugReports      : chr "http://r-forge.r-project.org/projects/lattice/" \cr
  .. .. .. ..$ Packaged        : chr "2014-04-03 11:25:19 UTC; deepayan" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-04-04 08:55:57" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:32:02 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/lattice/Meta/package.rds" \cr
  .. .. ..$ limma            :List of 15 \cr
  .. .. .. ..$ Package    : chr "limma" \cr
  .. .. .. ..$ Version    : chr "3.23.1" \cr
  .. .. .. ..$ Date       : chr "2014/10/20" \cr
  .. .. .. ..$ Title      : chr "Linear Models for Microarray Data" \cr
  .. .. .. ..$ Description: chr "Data analysis, linear models and differential expression for microarray data." \cr
  .. .. .. ..$ Author     : chr "Gordon Smyth [cre,aut], Matthew Ritchie [ctb], Jeremy Silver [ctb], James Wettenhall [ctb], Natalie Thorne [ctb], Davis McCarth"| __truncated__ \cr
  .. .. .. ..$ Maintainer : chr "Gordon Smyth <smyth@wehi.edu.au>" \cr
  .. .. .. ..$ License    : chr "GPL (>=2)" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.3.0), methods" \cr
  .. .. .. ..$ Suggests   : chr "statmod (>= 1.2.2), splines, locfit, MASS, ellipse, affy,vsn, AnnotationDbi, org.Hs.eg.db, GO.db, illuminaio, BiasedUrn" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ URL        : chr "http://bioinf.wehi.edu.au/limma" \cr
  .. .. .. ..$ biocViews  : chr "ExonArray, GeneExpression, Transcription,AlternativeSplicing, DifferentialExpression,DifferentialSplicing, GeneSetEnrichmen"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-10-21 01:24:05 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 09:47:06 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/limma/Meta/package.rds" \cr
  .. .. ..$ Matrix           :List of 26 \cr
  .. .. .. ..$ Package         : chr "Matrix" \cr
  .. .. .. ..$ Version         : chr "1.1-4" \cr
  .. .. .. ..$ Date            : chr "2014-06-14" \cr
  .. .. .. ..$ Priority        : chr "recommended" \cr
  .. .. .. ..$ Title           : chr "Sparse and Dense Matrix Classes and Methods" \cr
  .. .. .. ..$ Author          : chr "Douglas Bates <bates@stat.wisc.edu> and Martin Maechler" \cr
  .. .. .. ..$ Maintainer      : chr "Martin Maechler <mmaechler+Matrix@gmail.com>" \cr
  .. .. .. ..$ Contact         : chr "Doug and Martin <Matrix-authors@R-project.org>" \cr
  .. .. .. ..$ Description     : chr "Classes and methods for dense and sparse matrices and    operations on them using LAPACK and SuiteSparse." \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.15.2), methods" \cr
  .. .. .. ..$ Imports         : chr "graphics, grid, stats, utils, lattice" \cr
  .. .. .. ..$ Suggests        : chr "expm, MASS" \cr
  .. .. .. ..$ Enhances        : chr "MatrixModels, graph, SparseM, sfsmisc" \cr
  .. .. .. ..$ Encoding        : chr "UTF-8" \cr
  .. .. .. ..$ LazyData        : chr "no" \cr
  .. .. .. ..$ LazyDataNote    : chr "not possible, since we use data/*.R *and* our classes" \cr
  .. .. .. ..$ ByteCompile     : chr "yes" \cr
  .. .. .. ..$ BuildResaveData : chr "no" \cr
  .. .. .. ..$ License         : chr "GPL (>= 2)" \cr
  .. .. .. ..$ LicenseNote     : chr "The Matrix package includes libraries AMD, CHOLMOD,COLAMD, CSparse and SPQR from the SuiteSparse collection of TimDavis.  A"| __truncated__ \cr
  .. .. .. ..$ URL             : chr "http://Matrix.R-forge.R-project.org/" \cr
  .. .. .. ..$ Packaged        : chr "2014-06-14 16:22:31 UTC; maechler" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-06-15 20:26:23" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:32:13 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/Matrix/Meta/package.rds" \cr
  .. .. ..$ RCurl            :List of 21 \cr
  .. .. .. ..$ Package           : chr "RCurl" \cr
  .. .. .. ..$ Version           : chr "1.95-4.3" \cr
  .. .. .. ..$ Title             : chr "General network (HTTP/FTP/...) client interface for R" \cr
  .. .. .. ..$ Author            : chr "Duncan Temple Lang" \cr
  .. .. .. ..$ SystemRequirements: chr "libcurl (version 7.14.0 or higher)http://curl.haxx.se. On Linux systems, you will often have toexplicitly install libcurl-d"| __truncated__ \cr
  .. .. .. ..$ Maintainer        : chr "Duncan Temple Lang <duncan@r-project.org>" \cr
  .. .. .. ..$ Description       : chr "The package allows one to compose general HTTP requests        and provides convenient functions to fetch URIs, get & post "| __truncated__ \cr
  .. .. .. ..$ License           : chr "BSD" \cr
  .. .. .. ..$ Depends           : chr "R (>= 2.7.0), methods, bitops" \cr
  .. .. .. ..$ Imports           : chr "methods" \cr
  .. .. .. ..$ Suggests          : chr "Rcompression, XML" \cr
  .. .. .. ..$ URL               : chr "http://www.omegahat.org/RCurl" \cr
  .. .. .. ..$ Limitations       : chr "One doesn't yet have full control over the contents of aPOST form such as specifying files, content type. Errorhandling use"| __truncated__ \cr
  .. .. .. ..$ Note              : chr "This will be used generally to communicate with HTTP servers andspecifically in the SOAP package, HTML forms and for accessin"| __truncated__ \cr
  .. .. .. ..$ Check             : chr "To run R CMD check, one really needs network access. It is anHTTP client! This can be problematic as a Web site may be down"| __truncated__ \cr
  .. .. .. ..$ Collate           : chr "aclassesEnums.R bitClasses.R xbits.R base64.R binary.Sclasses.S curl.S curlAuthConstants.R curlEnums.R curlError.RcurlInfo."| __truncated__ \cr
  .. .. .. ..$ Packaged          : chr "2014-07-29 09:31:41 UTC; ligges" \cr
  .. .. .. ..$ NeedsCompilation  : chr "yes" \cr
  .. .. .. ..$ Repository        : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication  : chr "2014-07-29 11:53:00" \cr
  .. .. .. ..$ Built             : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:47:04 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/RCurl/Meta/package.rds" \cr
  .. .. ..$ Rsamtools        :List of 18 \cr
  .. .. .. ..$ Package    : chr "Rsamtools" \cr
  .. .. .. ..$ Type       : chr "Package" \cr
  .. .. .. ..$ Title      : chr "Binary alignment (BAM), FASTA, variant call (BCF), and tabixfile import" \cr
  .. .. .. ..$ Version    : chr "1.19.10" \cr
  .. .. .. ..$ Author     : chr "Martin Morgan, Herv\\'e Pag\\`es, Valerie Obenchain, Nathaniel Hayden" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer<maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Description: chr "This package provides an interface to the 'samtools',        'bcftools', and 'tabix' utilities (see 'LICENCE') for        m"| __truncated__ \cr
  .. .. .. ..$ URL        : chr "http://bioconductor.org/packages/release/bioc/html/Rsamtools.html" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0 | file LICENSE" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ Depends    : chr "methods, S4Vectors (>= 0.5.8), IRanges (>= 1.99.17),GenomeInfoDb (>= 1.1.3), GenomicRanges (>= 1.17.19), XVector(>= 0.5.3),"| __truncated__ \cr
  .. .. .. ..$ Imports    : chr "utils, BiocGenerics (>= 0.1.3), zlibbioc, bitops" \cr
  .. .. .. ..$ Suggests   : chr "GenomicAlignments, ShortRead (>= 1.19.10), GenomicFeatures,TxDb.Dmelanogaster.UCSC.dm3.ensGene, KEGG.db,TxDb.Hsapiens.UCSC."| __truncated__ \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges, XVector, Biostrings" \cr
  .. .. .. ..$ biocViews  : chr "DataImport, Sequencing, Coverage, Alignment, QualityControl" \cr
  .. .. .. ..$ Video      : chr "https://www.youtube.com/watch?v=Rfon-DQYbWA&list=UUqaMSQd_h-2EDGsU6WDiX0Q" \cr
  .. .. .. ..$ Packaged   : chr "2014-11-19 03:27:06 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:11:46 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/Rsamtools/Meta/package.rds" \cr
  .. .. ..$ RSQLite          :List of 18 \cr
  .. .. .. ..$ Package         : chr "RSQLite" \cr
  .. .. .. ..$ Version         : chr "1.0.0" \cr
  .. .. .. ..$ Title           : chr "SQLite Interface for R" \cr
  .. .. .. ..$ Authors@R       : chr "as.person(c(    \"Hadley Wickham <hadley@rstudio.com> [aut, cre]\",    \"David A. James [aut]\",    \"Seth Falcon [aut]\""| __truncated__ \cr
  .. .. .. ..$ Description     : chr "This package embeds the SQLite database engine in R and    provides an interface compliant with the DBI package.    The sou"| __truncated__ \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.10.0), DBI (>= 0.3.1), methods" \cr
  .. .. .. ..$ Suggests        : chr "testthat" \cr
  .. .. .. ..$ License         : chr "LGPL (>= 2)" \cr
  .. .. .. ..$ URL             : chr "https://github.com/rstats-db/RSQLite" \cr
  .. .. .. ..$ BugReports      : chr "https://github.com/rstats-db/RSQLite/issues" \cr
  .. .. .. ..$ Collate         : chr "'ConnectionExtensions.R' 'Connection.R' 'Driver.R' 'Connect.R''ConnectionRead.R' 'ConnectionTransactions.R''ConnectionWrite"| __truncated__ \cr
  .. .. .. ..$ Packaged        : chr "2014-10-24 18:38:41 UTC; hadley" \cr
  .. .. .. ..$ Author          : chr "Hadley Wickham [aut, cre],  David A. James [aut],  Seth Falcon [aut],  SQLite Authors [ctb] (for the included SQLite sour"| __truncated__ \cr
  .. .. .. ..$ Maintainer      : chr "Hadley Wickham <hadley@rstudio.com>" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-10-25 01:58:48" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-03 09:47:17 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/RSQLite/Meta/package.rds" \cr
  .. .. ..$ rtracklayer      :List of 16 \cr
  .. .. .. ..$ Package      : chr "rtracklayer" \cr
  .. .. .. ..$ Title        : chr "R interface to genome browsers and their annotation tracks" \cr
  .. .. .. ..$ Version      : chr "1.27.5" \cr
  .. .. .. ..$ Author       : chr "Michael Lawrence, Vince Carey, Robert Gentleman" \cr
  .. .. .. ..$ Depends      : chr "R (>= 2.10), methods, GenomicRanges (>= 1.19.3)" \cr
  .. .. .. ..$ Imports      : chr "XML (>= 1.98-0), BiocGenerics (>= 0.11.3), S4Vectors (>=0.5.4), IRanges (>= 2.1.6), XVector (>= 0.5.8), GenomeInfoDb(>= 1.1"| __truncated__ \cr
  .. .. .. ..$ Suggests     : chr "BSgenome (>= 1.33.4), humanStemCell, microRNA (>= 1.1.1),genefilter, limma, org.Hs.eg.db, hgu133plus2.db,BSgenome.Hsapiens."| __truncated__ \cr
  .. .. .. ..$ LinkingTo    : chr "S4Vectors, IRanges, XVector" \cr
  .. .. .. ..$ Description  : chr "Extensible framework for interacting with multiple genome   browsers (currently UCSC built-in) and manipulating   annotatio"| __truncated__ \cr
  .. .. .. ..$ Maintainer   : chr "Michael Lawrence <michafla@gene.com>" \cr
  .. .. .. ..$ License      : chr "Artistic-2.0" \cr
  .. .. .. ..$ ExtraLicenses: chr "The files in the 'src/ucsc' directory are licensed forall use by Jim Kent, in a manner compatible with the Artistic2.0 lice"| __truncated__ \cr
  .. .. .. ..$ Collate      : chr "io.R web.R ranges.R trackDb.R browser.R ucsc.R gff.R bed.Rwig.R utils.R bigWig.R chain.R quickload.R twobit.R fasta.Rtabix."| __truncated__ \cr
  .. .. .. ..$ biocViews    : chr "Annotation,Visualization,DataImport" \cr
  .. .. .. ..$ Packaged     : chr "2014-11-15 02:54:37 UTC; biocbuild" \cr
  .. .. .. ..$ Built        : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-16 03:49:19 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/rtracklayer/Meta/package.rds" \cr
  .. .. ..$ sendmailR        :List of 16 \cr
  .. .. .. ..$ Package         : chr "sendmailR" \cr
  .. .. .. ..$ Version         : chr "1.2-1" \cr
  .. .. .. ..$ Title           : chr "send email using R" \cr
  .. .. .. ..$ Description     : chr "Package contains a simple SMTP client which provides a        portable solution for sending email, including attachment, from"| __truncated__ \cr
  .. .. .. ..$ Authors@R       : chr "c(person(\"Olaf\", \"Mersmann\", role=c(\"aut\", \"cre\"),        email=\"olafm@p-value.net\"), person(\"Quinn\", \"Weber\",\"| __truncated__ \cr
  .. .. .. ..$ Depends         : chr "R (>= 3.0.0)" \cr
  .. .. .. ..$ Imports         : chr "base64enc" \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ LazyData        : chr "yes" \cr
  .. .. .. ..$ Packaged        : chr "2014-09-21 11:16:23 UTC; olafm" \cr
  .. .. .. ..$ Author          : chr "Olaf Mersmann [aut, cre],  Quinn Weber [ctb]" \cr
  .. .. .. ..$ Maintainer      : chr "Olaf Mersmann <olafm@p-value.net>" \cr
  .. .. .. ..$ NeedsCompilation: chr "no" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-09-21 13:32:55" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-10-14 13:47:56 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/sendmailR/Meta/package.rds" \cr
  .. .. ..$ snpStats         :List of 16 \cr
  .. .. .. ..$ Package    : chr "snpStats" \cr
  .. .. .. ..$ Title      : chr "SnpMatrix and XSnpMatrix classes and methods" \cr
  .. .. .. ..$ Version    : chr "1.17.0" \cr
  .. .. .. ..$ Date       : chr "2014-08-17" \cr
  .. .. .. ..$ Author     : chr "David Clayton <dc208@cam.ac.uk> " \cr
  .. .. .. ..$ Description: chr "Classes and statistical methods for large SNP association studies, extending the snpMatrix package" \cr
  .. .. .. ..$ Maintainer : chr "David Clayton <dc208@cam.ac.uk>" \cr
  .. .. .. ..$ Depends    : chr "R(>= 2.10.0), survival, Matrix, methods" \cr
  .. .. .. ..$ Imports    : chr "graphics, grDevices, stats, utils, BiocGenerics, zlibbioc" \cr
  .. .. .. ..$ Suggests   : chr "hexbin" \cr
  .. .. .. ..$ License    : chr "GPL-3" \cr
  .. .. .. ..$ Collate    : chr "ss.R contingency.table.R convert.R compare.R glm-test.Rimputation.R indata.R long.R misc.R ld.R mvtests.R pedfile.Routdata."| __truncated__ \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ biocViews  : chr "Microarray, SNP, GeneticVariability" \cr
  .. .. .. ..$ Packaged   : chr "2014-10-14 01:53:09 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-10-16 00:54:21 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/snpStats/Meta/package.rds" \cr
  .. .. ..$ splines          :List of 10 \cr
  .. .. .. ..$ Package    : chr "splines" \cr
  .. .. .. ..$ Version    : chr "3.2.0" \cr
  .. .. .. ..$ Priority   : chr "base" \cr
  .. .. .. ..$ Imports    : chr "graphics, stats" \cr
  .. .. .. ..$ Title      : chr "Regression Spline Functions and Classes" \cr
  .. .. .. ..$ Author     : chr "Douglas M. Bates <bates@stat.wisc.edu> and William N. Venables <Bill.Venables@csiro.au>" \cr
  .. .. .. ..$ Maintainer : chr "R Core Team <R-core@r-project.org>" \cr
  .. .. .. ..$ Description: chr "Regression spline functions and classes" \cr
  .. .. .. ..$ License    : chr "Part of R 3.2.0" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-10-14 13:30:58 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/splines/Meta/package.rds" \cr
  .. .. ..$ stringr          :List of 15 \cr
  .. .. .. ..$ Package         : chr "stringr" \cr
  .. .. .. ..$ Maintainer      : chr "Hadley Wickham <h.wickham@gmail.com>" \cr
  .. .. .. ..$ License         : chr "GPL-2" \cr
  .. .. .. ..$ Title           : chr "Make it easier to work with strings." \cr
  .. .. .. ..$ Type            : chr "Package" \cr
  .. .. .. ..$ Author          : chr "Hadley Wickham <h.wickham@gmail.com>" \cr
  .. .. .. ..$ Description     : chr "stringr is a set of simple wrappers that make R's string        functions more consistent, simpler and easier to use.  It doe"| __truncated__ \cr
  .. .. .. ..$ Version         : chr "0.6.2" \cr
  .. .. .. ..$ Depends         : chr "R (>= 2.14)" \cr
  .. .. .. ..$ Suggests        : chr "testthat (>= 0.3)" \cr
  .. .. .. ..$ Collate         : chr "'c.r' 'checks.r' 'count.r' 'detect.r' 'dup.r' 'extract.r''length.r' 'locate.r' 'match.r' 'modifiers.r' 'pad-trim.r''replace"| __truncated__ \cr
  .. .. .. ..$ Packaged        : chr "2012-12-05 21:47:03 UTC; hadley" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2012-12-06 08:39:59" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; ; 2014-08-13 10:45:29 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/stringr/Meta/package.rds" \cr
  .. .. ..$ survival         :List of 19 \cr
  .. .. .. ..$ Title           : chr "Survival Analysis" \cr
  .. .. .. ..$ Maintainer      : chr "Terry M Therneau <therneau.terry@mayo.edu>" \cr
  .. .. .. ..$ Priority        : chr "recommended" \cr
  .. .. .. ..$ Package         : chr "survival" \cr
  .. .. .. ..$ Version         : chr "2.37-7" \cr
  .. .. .. ..$ Depends         : chr "stats, utils, graphics, splines, R (>= 2.13.0)" \cr
  .. .. .. ..$ LazyData        : chr "Yes" \cr
  .. .. .. ..$ LazyLoad        : chr "Yes" \cr
  .. .. .. ..$ ByteCompile     : chr "Yes" \cr
  .. .. .. ..$ Authors@R       : chr "c(person(c(\"Terry\", \"M\"), \"Therneau\",                    email=\"therneau.terry@mayo.edu\",\t            role=c(\"aut"| __truncated__ \cr
  .. .. .. ..$ Description     : chr "survival analysis: descriptive statistics, two-sample  tests, parametric accelerated failure models, Cox model. Delayed entry"| __truncated__ \cr
  .. .. .. ..$ License         : chr "LGPL (>= 2)" \cr
  .. .. .. ..$ URL             : chr "http://r-forge.r-project.org" \cr
  .. .. .. ..$ Packaged        : chr "2014-01-22 13:49:28 UTC; therneau" \cr
  .. .. .. ..$ Author          : chr "Terry M Therneau [aut, cre],  Thomas Lumley [ctb, trl] (original S->R port and maintainer until 2009)" \cr
  .. .. .. ..$ NeedsCompilation: chr "yes" \cr
  .. .. .. ..$ Repository      : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication: chr "2014-01-22 17:48:51" \cr
  .. .. .. ..$ Built           : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:33:52 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/survival/Meta/package.rds" \cr
  .. .. ..$ VariantAnnotation:List of 18 \cr
  .. .. .. ..$ Package    : chr "VariantAnnotation" \cr
  .. .. .. ..$ Type       : chr "Package" \cr
  .. .. .. ..$ Title      : chr "Annotation of Genetic Variants" \cr
  .. .. .. ..$ Description: chr "Annotate variants, compute amino acid coding changes,        predict coding outcomes" \cr
  .. .. .. ..$ Version    : chr "1.13.12" \cr
  .. .. .. ..$ Authors@R  : chr "c(person(\"Valerie\", \"Obenchain\", role=c(\"aut\", \"cre\"),        email=\"vobencha@fhcrc.org\"),    person(\"Martin\", "| __truncated__ \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.8.0), methods, BiocGenerics (>= 0.7.7), GenomeInfoDb(>= 1.1.3), GenomicRanges (>= 1.19.7), Rsamtools (>= 1.19.4)" \cr
  .. .. .. ..$ Imports    : chr "utils, DBI, zlibbioc, Biobase, S4Vectors (>= 0.2.3), IRanges(>= 1.99.28), XVector (>= 0.5.6), Biostrings (>= 2.33.5),Annota"| __truncated__ \cr
  .. .. .. ..$ Suggests   : chr "RUnit, BSgenome.Hsapiens.UCSC.hg19,TxDb.Hsapiens.UCSC.hg19.knownGene,SNPlocs.Hsapiens.dbSNP.20110815,SNPlocs.Hsapiens.dbS"| __truncated__ \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges, XVector, Biostrings, Rsamtools" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ biocViews  : chr "DataImport, Sequencing, SNP, Annotation, Genetics" \cr
  .. .. .. ..$ Video      : chr "https://www.youtube.com/watch?v=Ro0lHQ_J--I&list=UUqaMSQd_h-2EDGsU6WDiX0Q" \cr
  .. .. .. ..$ Packaged   : chr "2014-11-20 04:00:42 UTC; biocbuild" \cr
  .. .. .. ..$ Author     : chr "Valerie Obenchain [aut, cre],  Martin Morgan [aut],  Michael Lawrence [aut],  Stephanie Gogarten [ctb]" \cr
  .. .. .. ..$ Maintainer : chr "Valerie Obenchain <vobencha@fhcrc.org>" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-23 12:12:37 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/VariantAnnotation/Meta/package.rds" \cr
  .. .. ..$ XML              :List of 19 \cr
  .. .. .. ..$ Package           : chr "XML" \cr
  .. .. .. ..$ Version           : chr "3.98-1.1" \cr
  .. .. .. ..$ Author            : chr "Duncan Temple Lang (duncan@r-project.org)" \cr
  .. .. .. ..$ Maintainer        : chr "Duncan Temple Lang <duncan@r-project.org>" \cr
  .. .. .. ..$ Title             : chr "Tools for parsing and generating XML within R and S-Plus." \cr
  .. .. .. ..$ Depends           : chr "R (>= 1.2.0), methods, utils" \cr
  .. .. .. ..$ Imports           : chr "methods" \cr
  .. .. .. ..$ Suggests          : chr "bitops, RCurl" \cr
  .. .. .. ..$ SystemRequirements: chr "libxml2 (>= 2.6.3)" \cr
  .. .. .. ..$ Description       : chr "This package provides many approaches for both reading and        creating XML (and HTML) documents (including DTDs), both lo"| __truncated__ \cr
  .. .. .. ..$ Note              : chr "In version 2.4.0 of this package, a new approach to garbagecollection has been implemented and it is experimental. You cand"| __truncated__ \cr
  .. .. .. ..$ URL               : chr "http://www.omegahat.org/RSXML" \cr
  .. .. .. ..$ License           : chr "BSD" \cr
  .. .. .. ..$ Collate           : chr "AAA.R DTD.R DTDClasses.R DTDRef.R SAXMethods.S XMLClasses.RapplyDOM.R assignChild.R catalog.R createNode.R dynSupports.Rerr"| __truncated__ \cr
  .. .. .. ..$ Packaged          : chr "2013-06-20 13:45:20 UTC; ripley" \cr
  .. .. .. ..$ NeedsCompilation  : chr "yes" \cr
  .. .. .. ..$ Repository        : chr "CRAN" \cr
  .. .. .. ..$ Date/Publication  : chr "2013-06-20 15:50:33" \cr
  .. .. .. ..$ Built             : chr "R 3.2.0; x86_64-apple-darwin13.3.0; 2014-08-13 10:46:02 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/XML/Meta/package.rds" \cr
  .. .. ..$ XVector          :List of 15 \cr
  .. .. .. ..$ Package    : chr "XVector" \cr
  .. .. .. ..$ Title      : chr "Representation and manpulation of external sequences" \cr
  .. .. .. ..$ Description: chr "Memory efficient S4 classes for storing sequences \"externally\"\t(behind an R external pointer, or on disk)." \cr
  .. .. .. ..$ Version    : chr "0.7.2" \cr
  .. .. .. ..$ Author     : chr "H. Pages and P. Aboyoun" \cr
  .. .. .. ..$ Maintainer : chr "H. Pages <hpages@fhcrc.org>" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure, DataRepresentation" \cr
  .. .. .. ..$ Depends    : chr "R (>= 2.8.0), methods, BiocGenerics (>= 0.11.3), S4Vectors (>=0.2.0), IRanges (>= 2.1.2)" \cr
  .. .. .. ..$ Imports    : chr "methods, BiocGenerics, S4Vectors, IRanges" \cr
  .. .. .. ..$ LinkingTo  : chr "S4Vectors, IRanges" \cr
  .. .. .. ..$ Suggests   : chr "Biostrings, drosophila2probe, RUnit" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0" \cr
  .. .. .. ..$ Collate    : chr "SharedVector-class.R SharedRaw-class.R SharedInteger-class.RSharedDouble-class.R XVector-class.R XRaw-class.RXInteger-class"| __truncated__ \cr
  .. .. .. ..$ Packaged   : chr "2014-11-06 05:08:02 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-11-07 15:41:56 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/XVector/Meta/package.rds" \cr
  .. .. ..$ zlibbioc         :List of 13 \cr
  .. .. .. ..$ Package    : chr "zlibbioc" \cr
  .. .. .. ..$ Type       : chr "Package" \cr
  .. .. .. ..$ Title      : chr "An R packaged zlib-1.2.5" \cr
  .. .. .. ..$ Version    : chr "1.13.0" \cr
  .. .. .. ..$ Author     : chr "Martin Morgan" \cr
  .. .. .. ..$ Maintainer : chr "Bioconductor Package Maintainer <maintainer@bioconductor.org>" \cr
  .. .. .. ..$ Description: chr "This package uses the source code of zlib-1.2.5 to create  libraries for systems that do not have these available via other"| __truncated__ \cr
  .. .. .. ..$ URL        : chr "http://bioconductor.org/packages/release/bioc/html/Zlibbioc.html" \cr
  .. .. .. ..$ License    : chr "Artistic-2.0 + file LICENSE" \cr
  .. .. .. ..$ LazyLoad   : chr "yes" \cr
  .. .. .. ..$ biocViews  : chr "Infrastructure" \cr
  .. .. .. ..$ Packaged   : chr "2014-10-14 02:39:15 UTC; biocbuild" \cr
  .. .. .. ..$ Built      : chr "R 3.2.0; x86_64-apple-darwin13.4.0; 2014-10-16 00:54:05 UTC; unix" \cr
  .. .. .. ..- attr(*, "class")= chr "packageDescription" \cr
  .. .. .. ..- attr(*, "file")= chr "/Users/stvjc/ExternalSoft/R-devel-dist/R.framework/Versions/3.2/Resources/library/zlibbioc/Meta/package.rds" \cr
  .. ..- attr(*, "class")= chr "sessionInfo" \cr
}
\details{
filtFDR was constructed on geuvStore contents, filtering to
MAF at least five percent and radius at most 500kbp.
rawFDR uses the entire geuvStore contents, with 1Mbp radius and
1 percent MAF lower bound
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(filtFDR)
filtFDR
}
\keyword{datasets}
