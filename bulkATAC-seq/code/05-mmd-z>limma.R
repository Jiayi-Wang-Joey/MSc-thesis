suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
    library(preprocessCore)
})

fun <- function(scores, norm) {
    activityScore <- scores$activityScore
    motif_score <- t(sapply(activityScore, \(motif) motif$motif_score))
    sample_id <- colnames(motif_score)
    if (norm=="TMM") {
        # y <- calcNormFactors(DGEList(motif_score))
        # nf <- y$samples$lib.size*y$samples$norm.factors
        # nf <- median(nf)/nf
        nf <- 1
    } else {
        nf <- 1
    }
    score <- t(sapply(activityScore, \(m) {
        ms <- m$motif_score[sample_id]*nf
        bg <- m$bg_scores[sample_id,]*nf
        dev_motif <- (ms - mean(ms))/mean(ms)
        dev_bg <- (bg-colMeans(bg))/colMeans(bg)
        (dev_motif - rowMeans(dev_bg))/rowSds(dev_bg)
    }))
    x <- score
    group_id <- substr(colnames(x), 1, nchar(colnames(x)) - 5)
    design <- model.matrix(~ group_id)
    fit <- eBayes(lmFit(x, design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res$name <- rownames(res)
    data.frame(res)
}

# dev_motif <- (motif_score-mean(motif_score))/mean(motif_score)
# dev_bg <- (bg_scores-colMeans(bg_scores))/colMeans(bg_scores)
# dev_z <- (dev_motif - rowMeans(dev_bg))/rowSds(dev_bg)