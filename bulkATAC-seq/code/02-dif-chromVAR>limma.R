suppressPackageStartupMessages({
    library(chromVAR)
    library(limma)
    library(SummarizedExperiment)
    library(BiocParallel)
})
register(MulticoreParam(8))
fun <- function(x, genome, motif) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x,
        genome = genome)
    bg <- getBackgroundPeaks(object = x, niterations = 2000)
    motif_ix <- matchMotifs(motif,
        x,
        genome = genome)
    dev <- chromVAR::computeDeviations(object = x,
        annotations = motif_ix,
        expectation = computeExpectations(x),
        background_peaks = bg)
    
    group_id <- substr(colnames(dev), 1, nchar(colnames(dev)) - 5)
    design <- model.matrix(~ group_id)
    fit <- eBayes(lmFit(assay(dev, "z"), design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    ids <- match(rownames(res), rownames(dev))
    res$name <- rowData(dev)$name[ids]
    res$rank <- seq_len(nrow(res))
    list(res=res, z=assay(dev, "z"))
}

