suppressPackageStartupMessages({
    library(MLVSBM)
})

fun <- \(x) {
    res <- ARI(x$cell_type, x$cluster_id)
    data.frame(eva_val=res, eva="ARI", row.names=NULL)
}


