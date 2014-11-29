checkStats = function() {
# as of 29 XI 2014
# [1] "cisAssoc"          "clipPCs"           "directPlot"       
# [4] "enumerateByFDR"    "eqBox2"            "eqDesc2"          
# [7] "getFDRfunc"        "getTab"            "regressOut"       
#[10] "setFDRfunc"        "storeToFDR"        "storeToFDRByProbe"
#[13] "storeToHist"       "storeToQuantiles"  "txsPlot"        

#
# test cisAssoc
#
   require(GenomeInfoDb)
   require(geuvPack)
   require(Rsamtools)
   data(geuFPKM)
   lgeu = geuFPKM[ which(seqnames(geuFPKM)=="chr20"), ]
   seqlevelsStyle(lgeu) = "NCBI"
   tf20 = TabixFile(system.file("vcf/lit20.vcf.gz", package="gQTLstats"))
   lgeue = clipPCs(lgeu[,which(lgeu$popcode=="CEU")], 1:2)
   litc = cisAssoc(lgeue[c(162,201),], tf20, nperm=2, 
         lbmaf=.05, cisradius=50000)
   checkTrue(length(litc)==497)
   checkTrue(min(litc$MAF)>=.05)
   checkTrue(max(litc$mindist) <= 50000)
   checkTrue(max(litc$chisq) > 6.4 )
   checkTrue(median(litc$chisq) > 2.97 & median(litc$chisq) < 2.98 )
   checkTrue( all(c("sessInfo", "init.Random.seed") %in% names(metadata(litc)) ))

#
# test clipPCs
#
   gmns = apply(assay(lgeue),1,mean)
   checkTrue( max(abs(gmns)) < 1e-10 )

#
# test directPlot -- should have a plot=FALSE param
#

#
# test enumerateByFDR
#
     require(geuvStore)
     require(gQTLBase)
     re = partialRegistry()
     st = ciseStore(re, FALSE, FALSE)
     data(filtFDR)
     filtEnum = enumerateByFDR( st, filtFDR,
        filter=function(x)x[which(x$mindist <= 500000 & x$MAF >= 0.05)],
        ids=1:3 )
     checkTrue( all(c("enumCall", "enumSess", "fdrCall") %in% 
           names(metadata(filtEnum))))
     checkTrue( length(filtEnum) == 2387 )
     checkTrue( max(filtEnum$estFDR) < 0.05 )
     checkTrue( min(filtEnum$estFDR) >= 0.0 )
     checkTrue( length(unique(filtEnum$probeid)) == 29 )
     checkTrue( length(unique(filtEnum$snp)) == 2366 )

#
# test eqDesc2 (related eqBox2)
#

     require(GenomicRanges)
     mygr = GRanges("1", IRanges(54683925, width=1))
     gene = "ENSG00000231581.1"
     library(geuvPack)
     tf = TabixFile(system.file("vcf/small_1.vcf.gz", package="gQTLstats"))
     ed2 = eqDesc2(gene, se=geuFPKM, tf, mygr )
     checkTrue(all(names(ed2)==c("A/A", "A/B", "B/B")))
     checkTrue(all(as.numeric(ed2)==c(318, 94, 9)))

#
# test getFDRfunc
#
     checkTrue(is.function(getFDRfunc(filtFDR)))
     checkTrue(names(formals(getFDRfunc(filtFDR)))=="assoc")

}
