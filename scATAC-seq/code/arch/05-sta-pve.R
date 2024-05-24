suppressPackageStartupMessages({
    library(matrixStats)
    library(variancePartition)
    library(edgeR)
    library(limma)
    library(motifmatchr)
    source("~/chromVAR2/code/utils.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(chromVAR)
})

fun <- \(x, m) {
    # subset selected features
    motif <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome=BSgenome.Hsapiens.UCSC.hg38)
    y <- t(assay(motif_ix))%*%assay(x, "counts")
    #counts <- assay(x,"counts")
    design <- model.matrix(~ x$group_id)
    #rownames(y) <- seq_len(nrow(y))
    y <- calcNormFactors(DGEList(y))
    y <- voom(y,design)
    # fit LMM to estimate fraction of variance
    # attributable to cluster, sample, group
    cd <- data.frame(colData(x))
    mod <- ~(1|group_id)
    res <- fitExtractVarPartModel(y, mod, cd, quiet=TRUE)
    rownames(res) <- rownames(y)
    rbind(
        data.frame(sta="PVE_all", sta_val=mean(res$group_id)),
        data.frame(sta="PVE_m", sta_val=res[m, "group_id"])
    )
}