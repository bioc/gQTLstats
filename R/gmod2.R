
gmod2 = function (sym, genome = "hg19", orgDb,
   collector=exonsBy, verbose=FALSE) 
{
    if (missing(orgDb)) {
      message("assuming Homo.sapiens for gene models")
      orgDb = Homo.sapiens
     }
    txdb = TxDb(orgDb)
    num = AnnotationDbi::select(orgDb, keys=sym, keytype="SYMBOL", 
          columns="ENTREZID")$ENTREZID
    collector(txdb, by = "gene")[[num]]
}

