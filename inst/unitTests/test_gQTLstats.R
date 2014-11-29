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
#
# test getTab
#
     checkTrue(all(dim(getTab(filtFDR)) == c(27,4)))
#
#    test regressOut
#
     ro = regressOut(geuFPKM, ~popcode)
     checkTrue(max(abs(coef(lm(assay(ro)[1,]~geuFPKM$popcode))))<1e-16)
#
#    test setFDRfunc
#
     filtFDR@FDRfunc = NULL
     filtFDR = setFDRfunc( filtFDR )
     checkTrue(is.function(getFDRfunc(filtFDR)))
     checkTrue(max(abs(coef(filtFDR@FDRmodel)-
          c(-0.147156183359907, -0.0712895014838215))) < 1e-6 )
#
# test storeToFDR
#
      reg = partialRegistry()
      store = ciseStore(reg, addProbeMap=FALSE, addRangeMap=FALSE)
      stf = storeToFDR(store)
      checkTrue(all(dim(getTab(stf))==c(1004,4)))
      checkTrue(is.null(getFDRfunc(stf)))
      checkTrue(max(getTab(stf)$assoc)>252.89)
#
# skip storeToFDRByProbe as currently slow
#
# test storeToHist
#
      hh = storeToHist( store, breaks=c(0,1,2,4,8,350))
#      checkTrue(sum(hh$counts) == 3*5656290)
      checkTrue(all(hh$counts == c(11677429L, 2691328L, 1832875L, 667539L, 99699L)))
      dmfilt = function(x) x[ which(x$MAF >= .05 & x$mindist <= 5e4) ]
      hf = storeToHist( store, breaks = c(0,1,2,4,8,350), filter=dmfilt )
      checkTrue( sum(hf$counts) == 692397)
#
# test storeToQuantiles
#
      sq = storeToQuantiles(store, "chisq", seq(.1,.9,.1))
      qtarg = structure(c(0.0146886028850068, 0.0600514288350784, 0.13982186386647, 
0.261229560601213, 0.435622826467435, 0.682153612200419, 1.04084584125386, 
1.61960915615195, 2.76456101386141), .Names = c("10%", "20%", 
"30%", "40%", "50%", "60%", "70%", "80%", "90%"))
      checkTrue(max(abs(sq-qtarg))<1e-6)
      sqf = storeToQuantiles(store, "chisq", seq(.1,.9,.1), filter=dmfilt)

      targ2 = structure(c(0.0237787672262377, 0.0989828833729177, 0.233430426470715, 
0.434455200684701, 0.740987113540066, 1.1698714242342, 1.8634039013722, 
3.05695105975853, 6.1302713856499), .Names = c("10%", "20%", 
"30%", "40%", "50%", "60%", "70%", "80%", "90%"))
      checkTrue(max(abs(sqf-targ2))<1e-6)

}
