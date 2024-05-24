suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

se <- readRDS(args$dat)
source(args$fun)

if (wcs$sta!="asw") {
    sta <- fun(x=se, m=wcs$mot)
} else {
    sta <- fun(se)
}

if (!is.null(wcs$smt)) {
    res <- data.frame(motif=wcs$mot, mode="weight", peakWeight=wcs$pkw,
        status=wcs$mdr, smooth=wcs$smt, sta)
} else {
    res <- data.frame(motif=wcs$mot, mode="total", peakWeight="none",
        status="unmoderated", smooth="none", sta)
}

saveRDS(res, args$res)