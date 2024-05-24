suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    source("code/utils.R")
    library(MotifDb)
})


se <- readRDS(args$dat)
source(args$fun)

Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
res <- fun(se, genome = BSgenome.Hsapiens.UCSC.hg38, 
    motif = Hmotifs)

if (wcs$dif=="chromVAR>limma") res <- res$res

if (is.null(wcs$smt)) {
    df <- data.frame(res, test=wcs$mtf,  
      effect=wcs$eft, dif=wcs$dif, 
      mode="total", smooth="none", 
      peakWeight="none", row.names=NULL)
} else {
    df <- data.frame(res, test=wcs$mtf, 
      effect=wcs$eft, dif=wcs$dif, 
      mode="weight", smooth=wcs$smt, 
      peakWeight=wcs$pkw, row.names=NULL)
}

saveRDS(df, args$res)