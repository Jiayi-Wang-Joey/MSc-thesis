suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
    library(pryr)
})

fun <- function(x, genome, motif) {
    #x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome = genome)
    y <- t(assay(motif_ix))%*%assay(x, "counts")
    group_id <- substr(colnames(y), 1, nchar(colnames(y)) - 5)
    design <- model.matrix(~ group_id)
    y <- calcNormFactors(DGEList(y))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res$name <- rownames(res)
    res 
    
}