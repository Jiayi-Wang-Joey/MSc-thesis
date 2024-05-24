suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
})

fun <- function(scores, norm) {
    activityScore <- scores$activityScore
    motif_score <- t(sapply(activityScore, \(motif) motif$motif_score))
    sample_id <- colnames(motif_score)
    # if (norm=="TMM") {
    #     y <- calcNormFactors(DGEList(score))
    #     x <- cpm(y)
    # } else {
    #     x <- score*1e9
    # }
    x <- motif_score
    group_id <- substr(colnames(x), 1, nchar(colnames(x)) - 5)
    design <- model.matrix(~ group_id)
    x <- calcNormFactors(DGEList(x))
    #x <- cpm(x)
    x <- voom(x,design)
    fit <- eBayes(lmFit(x, design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res$name <- rownames(res)
    data.frame(res)
}
