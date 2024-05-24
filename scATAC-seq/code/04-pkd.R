se <- readRDS(args$dat)
source(args$fun)

res <- fun(se)

if (is.null(wcs$smt)) {
    df <- data.frame(res, motif=wcs$mot, 
        dif=wcs$pkd, mode="total",  
        status="unmoderated",
        smooth="none",
        peakWeight="none",
        row.names=rownames(res))
} else {
    df <- data.frame(res, motif=wcs$mot, 
        dif=wcs$pkd, mode="weight", 
        status=wcs$mdr,
        smooth=wcs$smt, 
        peakWeight=wcs$pkw,
        row.names=rownames(res))
}

saveRDS(df, args$res)
