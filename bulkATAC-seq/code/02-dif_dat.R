suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("code/utils.R")
    library(MotifDb)
})

mice <- c("BANP", "NR1H3", "NR1H4")
m <- wcs$mot
se <- readRDS(args$dat)
source(args$fun)

if (m %in% mice) {
    Mmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Mmusculus")
    if (m=="BANP") {
        banp <- readRDS("data/BANP.PFMatrix.rds")
        Mmotifs$BANP <- banp
    } else if (m=="NR1H3") {
        Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
        Mmotifs$NR1H3 <- Hmotifs$NR1H3
    }
    res <- fun(se, 
        genome = BSgenome.Mmusculus.UCSC.mm10, 
        motif = Mmotifs)
    if (!is.null(res$z)) {
        saveRDS(res$z, paste0("outs/dat/chromVAR-z-", wcs$mot, ".rds"))
        res <- res$res
    }
} else {
    Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    res <- fun(se, 
        genome = BSgenome.Hsapiens.UCSC.hg38, 
        motif = Hmotifs)
    if (!is.null(res$z)) {
        saveRDS(res$z, paste0("outs/dat/chromVAR-z-", wcs$mot, ".rds"))
        res <- res$res
    }
}

if (m=="MYC") { 
    target <- paste("MYC", "MAX", sep=",")
} else if (m=="NR1H4"){
    target <- paste("NR1H4","RXRA","RXRB",sep=",")
} else if (m=="ESR1"){
    target <- paste("ESR1", "ESR2", sep=",")
} else {
    target <- m
}
if (is.null(wcs$smt)) {
    df <- data.frame(res, truth=target, 
        dif=wcs$dif, mode="total", moderated="unmoderated",
        smooth="none", peakWeight="none", row.names=NULL)
} else {
    df <- data.frame(res, truth=target, 
        dif=wcs$dif, mode="weight", moderated=wcs$mdr,
        smooth=wcs$smt, peakWeight=wcs$pkw,row.names=NULL)
}

saveRDS(df, args$res)
