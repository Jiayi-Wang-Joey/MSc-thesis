suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
})

fun <- function(x, genome, motif) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome = genome)
    y <- t(assay(motif_ix))%*%assay(x, "counts")
    mm <- model.matrix(~ x$group_id)
    y <- calcNormFactors(DGEList(y))
    v <- voom(y,mm)
    corfit <- duplicateCorrelation(v, mm, block=x$sample_id)
    v <- voom(y, mm, block=x$sample_id, correlation=corfit$consensus.correlation)
    fit <- lmFit(v, mm, block=x$sample_id, correlation=corfit$consensus.correlation)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res 
    
}