
gmod2 = function (sym, genome = "hg19", orgDb=Homo.sapiens,
   collector=exonsBy, verbose=FALSE) 
{
    rend = suppressPackageStartupMessages
    if (verbose) rend = force
    rend({
    require(txn <- gsub("%%G%%", genome, "TxDb.Hsapiens.UCSC.%%G%%.knownGene"),
      character.only=TRUE)
    require(deparse(substitute(orgDb)), character.only=TRUE)
    })
    txdb = get(txn)
    num = select(orgDb, keys=sym, keytype="SYMBOL", 
          columns="ENTREZID")$ENTREZID
    collector(txdb, by = "gene")[[num]]
}

