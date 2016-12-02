# symbols as of Dec 2 2016
#
# [1] "AllAssoc"             "cisAssoc"             "cisCount"            
# [4] "cisEsts"              "clipPCs"              "collapse_multiPerm"  
# [7] "describe"             "directPlot"           "distToGene"          
#[10] "enumerateByFDR"       "eqBox2"               "eqBox3"              
#[13] "eqDesc2"              "filteredDFwPerm"      "getFDRfunc"          
#[16] "getTab"               "gQTLs"                "gQTLswarm"           
#[19] "manhWngr"             "mixedVCFtoSnpMatrix"  "pifdr"               
#[22] "pifdr2"               "pifdr3"               "prep.cisAssocNB"     
#[25] "qqStore"              "queryVCF"             "regressOut"          
#[28] "senstab"              "setFDRfunc"           "storeToFDR"          
#[31] "storeToFDRByProbe"    "storeToHist"          "storeToMaxAssocBySNP"
#[34] "storeToQuantiles"     "table_sensobj_thresh" "tqbrowser"           
#[37] "transAssoc"           "transBrowse"          "transBrowse2"        
#[40] "TransChunk"           "TransStore"           "transTable"          
#[43] "tsByRank"             "tsByRankAccum"        "tsIndex.reg"         
#[46] "txsPlot"             

# declared tested at that date

# test cisAssoc
# test clipPCs
# test directPlot -- should have a plot=FALSE param
# test enumerateByFDR
# test eqDesc2 (related eqBox2)
# test getFDRfunc
# test getTab
#    test regressOut
#    test setFDRfunc
# test storeToFDR
# test storeToHist
# test storeToQuantiles


checkStats = function() {
#
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
   require(gQTLstats)
   data(geuFPKM)
   lgeu = geuFPKM[ which(seqnames(geuFPKM)=="chr20"), ]
   seqlevelsStyle(lgeu) = "NCBI"
   tf20 = TabixFile(system.file("vcf/lit20.vcf.gz", package="gQTLstats"))
   lgeue = clipPCs(lgeu[,which(lgeu$popcode=="CEU")], 1:2)
   litc = cisAssoc(lgeue[c(162,201),], tf20, nperm=2, 
         lbmaf=.05, cisradius=50000, lbgtf=-.01)
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
     require(geuvStore2)
     require(gQTLBase)
     st = makeGeuvStore2()
     data(filtFDR)
     filtEnum = enumerateByFDR( st, filtFDR,
        filter=function(x)x[which(x$mindist <= 500000 & x$MAF >= 0.05)],
        ids=1:3 )
     checkTrue( all(c("enumCall", "enumSess", "fdrCall") %in% 
           names(metadata(filtEnum))))
     checkTrue( length(filtEnum) == 989 )
     checkTrue( max(filtEnum$estFDR) < 0.05 )
     checkTrue( min(filtEnum$estFDR) >= 0.0 )
     checkTrue( length(unique(filtEnum$probeid)) == 18 )
     checkTrue( length(unique(filtEnum$snp)) == 985 )

#
# test eqDesc2 (related eqBox2)
#

     require(SummarizedExperiment)
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
newcoef = structure(c(0.57487958216754, 0.712226903412448, 1.21462315263198, 
0.158844281416502, 0.937680032995375, -0.0333469006067456, 0.52226548981198, 
0.277968914437912, 1.40712654986676, 0.53461376374711), .Names = c("(Intercept)", 
"s(assoc).1", "s(assoc).2", "s(assoc).3", "s(assoc).4", "s(assoc).5", 
"s(assoc).6", "s(assoc).7", "s(assoc).8", "s(assoc).9"))
     checkTrue(all.equal(coef(filtFDR@FDRmodel), newcoef))
#
# test storeToFDR
#
      reg = st@reg
      store = ciseStore(reg, 1:160, addProbeMap=FALSE, addRangeMap=FALSE)
      stf = storeToFDR(store)
      checkTrue(all(dim(getTab(stf))==c(1004,4)))
      checkTrue(is.null(getFDRfunc(stf)))
      checkTrue(max(getTab(stf)$assoc)>242.4)
#
# skip storeToFDRByProbe as currently slow
#
# test storeToHist
#
      hh = storeToHist( store, breaks=c(0,1,2,4,8,350))
      newstorecounts = c(25565239L, 5871108L, 4005041L, 1451170L, 206558L)
      checkTrue(sum(hh$counts) == sum(newstorecounts))
      checkTrue(all(hh$counts == newstorecounts))
      dmfilt = function(x) x[ which(x$MAF >= .05 & x$mindist <= 5e4) ]
      hf = storeToHist( store, breaks = c(0,1,2,4,8,350), filter=dmfilt )
      newhfcounts = c(1266283L, 297891L, 207421L, 74220L, 9781L)

      checkTrue( sum(hf$counts) == sum(newhfcounts) )
      checkTrue( all(hf$counts == newhfcounts) )
#
# test storeToQuantiles
#
      sq = storeToQuantiles(store, "chisq", seq(.1,.9,.1))
      newqtarg = structure(c(0.0146639630699681, 0.0601324043720933, 0.139374460408643, 
0.259530217530304, 0.432646613334439, 0.679787086216459, 1.04083586926541, 
1.61469326410928, 2.75963121745673), .Names = c("10%", "20%", 
"30%", "40%", "50%", "60%", "70%", "80%", "90%"))
      checkTrue(max(abs(sq-newqtarg))<1e-6)
      sqf = storeToQuantiles(store, "chisq", seq(.1,.9,.1), filter=dmfilt)
      targ2 = structure(c(0.024627608679227, 0.0977522914503348, 0.227309798172982, 
0.427344130590394, 0.72735521049526, 1.14331130249013, 1.80577209887343, 
2.92063094652981, 5.74538744443098), .Names = c("10%", "20%", 
"30%", "40%", "50%", "60%", "70%", "80%", "90%"))
      checkTrue(max(abs(sqf-targ2))<1e-6)

}
checkStats()

checkMixedVcfProc = function() {
 require("snpStats")
 require("VariantAnnotation")
 fn = system.file("vcf/polytypeSNV.vcf", package="gQTLstats")
 vv = readVcf(fn, genome="hg19")
 ans = mixedVCFtoSnpMatrix(vv, FALSE)$genotypes@.Data
 checkTrue(ans[1,4] == as.raw(0xfd))
}
checkMixedVcfProc()

checkTsByRank = function() {
 if (require(geuvStore2) && require(doParallel)) {
       registerDoSEQ()
       r17 = g17transRegistry()
       r18 = g18transRegistry()
       g1718 = TransStore(list(r17, r18))
       }
 tt1 = tsByRank(g1718, 1, mcol2keep=c("REF", "snp", "MAF")) # must limit
        # as geuvStore2 is legacy ...
 checkTrue(length(tt1) == 509728)
 keyfields = c("snp", "MAF", "feats", "scores", "permscores", 
        "obsdist", "permdist")  # eventually add z.HWE
 checkTrue(all(keyfields %in% names(mcols(tt1))))
 checkTrue(abs(mean(tt1$scores)-26.29626)<.001)
 checkTrue(abs(mean(tt1$MAF)-0.1718)<.001)
}
checkTsByRank()

