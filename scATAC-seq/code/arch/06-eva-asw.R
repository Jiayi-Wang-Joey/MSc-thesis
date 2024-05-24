suppressPackageStartupMessages({
    library(cluster)
    library(PCAtools)
})

fun <- \(x) {
    #y <- colData(x)
    idx <- grep("^LSI_(?:[2-9]|[12][0-9]|30)$", 
        colnames(x), 
        value = TRUE)
    y <- x[,idx]
    res <- silhouette(as.integer(factor(x$cell_type)), dist(y))
    data.frame(eva="ASW", eva_val=mean(res[,3]))
}