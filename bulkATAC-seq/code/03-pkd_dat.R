suppressPackageStartupMessages({
    source("code/utils.R")
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues=TRUE)
m <- wcs$mot
se <- readRDS(args$dat)
source(args$fun)

res <- fun(se)

if (is.null(wcs$smt)) {
    df <- data.frame(res, test=m, 
        dif=wcs$pkd, mode="total",  moderated="unmoderated",
        smooth="none", peakWeight="none", row.names=rownames(res))
} else {
    df <- data.frame(res, test=m, 
        dif=wcs$pkd, mode="weight", moderated=wcs$mdr,
        smooth=wcs$smt, peakWeight=wcs$pkw, row.names=rownames(res))
}

saveRDS(df, args$res)
