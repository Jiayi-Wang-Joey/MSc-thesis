suppressPackageStartupMessages({
    library(chromVAR)
    library(motifmatchr)
    source("~/chromVAR2/code/utils.R")
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(matrixStats)
    library(edgeR)
    library(reshape2)
})


fun <- \(x, m) {
    motif <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome=BSgenome.Hsapiens.UCSC.hg38)
    y <- t(assay(motif_ix))%*%assay(x, "counts")
    y <- calcNormFactors(DGEList(y))
    y <- cpm(y)
    vars <- sapply(unique(x$sample_id), \(i) {
        ids <- which(x$sample_id==i)
        var_ctrl <- rowSds(y[,ids])
    })
    ctrl <- which(colData(x)$group_id=="NT")
    kd <- which(colData(x)$group_id==m)
    #ratio <- abs(mean(y[m,kd])-mean(y[m,ctrl]))/mean(vars[m,])
    ratio <- abs(rowMeans(y[,kd])-rowMeans(y[,ctrl]))/rowMeans(vars)
    df <- melt(vars)
    
    cvs <- sapply(unique(x$sample_id), \(i) {
        ids <- which(x$sample_id==i)
        rowSds(y[,ids])/rowMeans(y[,ids])
    })
    cv <- melt(cvs)
    rbind(
        data.frame(sta="within-sample_variance", 
            sta_val=df$value, sample_id=df$Var2, 
            motif_id=df$Var1),
        data.frame(sta="ratio", sta_val=ratio, sample_id=NA, 
            motif_id=names(ratio)),
        data.frame(sta="CV", sta_val=cv$value, sample_id=cv$Var2,
            motif_id=cv$Var1))
}