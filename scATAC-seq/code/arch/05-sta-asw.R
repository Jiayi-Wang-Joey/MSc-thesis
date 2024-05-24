suppressPackageStartupMessages({
    library(cluster)
    library(PCAtools)
})

fun <- \(x) {
    y <- colData(x)
    idx <- grep("^LSI_(?:[2-9]|[12][0-9]|30)$", 
        colnames(colData(x)), 
        value = TRUE)
    y <- y[,idx]
    res <- silhouette(as.integer(factor(x$group_id)), dist(y))
    data.frame(sta="ASW", sta_val=mean(res[,3]))
}
