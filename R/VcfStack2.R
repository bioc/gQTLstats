VcfStack2 <- function(files=NULL, seqinfo=NULL, colData=NULL)
{
    if (is.null(files)) {
        files <- VcfFileList()
        header <- NULL
    } else {
        if (class(files) != "VcfFileList")
            files = VcfFileList(files)
#        files = indexVcf(files)
        header <- scanVcfHeader(files[[1]])
    
    }   

    if (is.null(seqinfo)) {
        seqinfo <- if (length(files)) {
            seqinfo(files)
        } else Seqinfo()
    }   

    if (is.null(colData) && length(files)) {
        colData <- DataFrame(row.names=samples(header))
    } else {
        colData <- as(colData, "DataFrame")
    }   

    if (is.null(rownames(colData)) && length(files))
         stop("specify rownames in 'colData'")

    new("VcfStack", files=files, colData=colData, seqinfo=seqinfo)
}
