suppressPackageStartupMessages({
    library(chromVAR)
    library(limma)
    library(SummarizedExperiment)
})

fun <- function(x, genome, motif) {
    counts <- assay(x, "counts")
    x <- filterPeaks(x, non_overlapping = TRUE)
    x <- addGCBias(x,
        genome = genome)
    bg <- getBackgroundPeaks(object = x)
    motif_ix <- matchMotifs(motif,
        x,
        genome = genome)
    dev <- chromVAR::computeDeviations(object = x,
        annotations = motif_ix,
        expectation = computeExpectations(x),
        background_peaks = bg)
    
    #group_id <- substr(colnames(dev), 1, nchar(colnames(dev)) - 5)
    design <- model.matrix(~ x$group_id)
    fit <- eBayes(lmFit(assay(dev, "z"), design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res 
}