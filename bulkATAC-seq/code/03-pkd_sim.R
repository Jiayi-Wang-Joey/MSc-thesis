suppressPackageStartupMessages({
    source("code/utils.R")
})

se <- readRDS(args$dat)
source(args$fun)
res <- fun(se)

if (is.null(wcs$smt)) {
    df <- data.frame(res, test=wcs$mtf,  
        effect=wcs$eft, dif=wcs$pkd,
        mode="total", smooth="none", 
        peakWeight="none",
        row.names=NULL)
} else {
    df <- data.frame(res, test=wcs$mtf, 
        effect=wcs$eft, dif=wcs$pkd,
        mode="weight", smooth=wcs$smt, 
        peakWeight=wcs$pkw,
        row.names=NULL)
}

saveRDS(df, args$res)



