suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

se <- readRDS(args$dat)
source(args$fun)
eva <- fun(se)

if (!is.null(wcs$smt)) {
    res <- data.frame(mode="weight", peakWeight=wcs$pkw,
        status=wcs$mdr, smooth=wcs$smt, eva)
} else {
    res <- data.frame(mode="total", peakWeight="none",
        status="unmoderated", smooth="none", eva)
}

saveRDS(res, args$res)