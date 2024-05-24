suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
})

fun <- function(x) {
    counts <- assay(x,"counts")
    design <- model.matrix(~ x$group_id)
    rownames(counts) <- seq_len(nrow(counts))
    y <- calcNormFactors(DGEList(counts))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf, sort.by = "none")
    res 
}
