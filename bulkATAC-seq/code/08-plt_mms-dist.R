#args <- list(list.files("outs/dat", "^mms-.*", full.names=TRUE), "plts/mmd-heatmap.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(edgeR)
    library(reshape2)
    library(dplyr)
    library(matrixStats)
    library(patchwork)
})

res <- lapply(args[[1]], \(i) {
    x <- readRDS(i)
    wcs <- x[[1]]
    scores <- x[[2]]
    activityScore <- scores$activityScore
    motif_score <- t(sapply(activityScore, \(motif) motif$motif_score))
    sample_id <- colnames(motif_score)
    if (wcs$nrm=="TMM") {
        y <- calcNormFactors(DGEList(motif_score))
        nf <- y$samples$lib.size*y$samples$norm.factors
        nf <- median(nf)/nf
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
    df <- melt(score)
    colnames(df) <- c("motif", "sample", "z")
    data.frame(df, target=wcs$mot, method=paste(wcs$nrm,wcs$sym,sep=","))
})
mm <- do.call(rbind, res) 

cv <- lapply(unique(mm$target), \(m) {
    z <- readRDS(paste0("outs/dat/chromVAR-z-", m, ".rds"))
    res <- melt(z)
    colnames(res) <- c("motif","sample","z")
    data.frame(res, target=m, method="chromVAR>limma")
})
cv <- do.call(rbind,cv)
df <- rbind(mm,cv)

ps <- lapply(split(df,df$target), \(fd) { 
    # lapply(split(fd,fd$method), \(f) {
    #     ggplot(f, aes(x=z,col=sample)) + 
    #         geom_density() +
    #         scale_color_brewer(palette = "Paired") +
    #         theme_minimal(6) + 
    #         ggtitle(f$method[1]) +
    #         geom_abline(v=0, col="black", lty=2, alpha=0.5)
    # }) |>  wrap_plots(ncol=3) + 
    #     plot_layout(guides = "collect") + plot_annotation(fd$target[1]) 
        ggplot(fd, aes(x=z,col=sample)) +
            geom_density() +
            facet_wrap(~method, ncol=3) +
            scale_color_brewer(palette = "Paired") +
            theme_bw() + 
            geom_vline(xintercept = 0, col="black", lty=2, alpha=0.5)
}) 

pdf(args[[2]], onefile=TRUE, width=15.5, height=5)
for (p in ps) print(p); dev.off()
