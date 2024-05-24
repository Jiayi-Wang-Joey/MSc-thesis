suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("~/chromVAR2/code/utils.R")
    library(MotifDb)
})

m <- wcs$mot
se <- readRDS(args$dat)
source(args$fun)


Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
res <- fun(se, 
    genome = BSgenome.Hsapiens.UCSC.hg38, 
    motif = Hmotifs)

if (is.null(wcs$smt)) {
    df <- data.frame(res, motif=m, 
        dif=wcs$dif, mode="total",
        moderated="unmoderated",
        smooth="none", peakWeight="none",
        row.names=rownames(res))
} else {
    df <- data.frame(res, motif=m,
        dif=wcs$dif, mode="weight", 
        moderated=wcs$mdr,
        smooth=wcs$smt, peakWeight=wcs$pkw,
        row.names=rownames(res))
}

saveRDS(df, args$res)