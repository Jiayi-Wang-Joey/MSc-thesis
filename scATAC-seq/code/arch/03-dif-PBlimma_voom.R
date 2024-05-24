suppressPackageStartupMessages({
    library(edgeR)
    library(poolr)
    library(scuttle)
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

fun <- \(x, motif, genome) {
    x <- filterPeaks(x, non_overlapping = TRUE)
    ids <- colData(x)[c("sample_id")]
    z <- aggregateAcrossCells(x, ids)
    motif_ix <- matchMotifs(motif, x, genome)
    y <- t(assay(motif_ix))%*%assay(z, "counts")
    design <- model.matrix(~ z$group_id)
    y <- calcNormFactors(DGEList(y))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res 
    
}
