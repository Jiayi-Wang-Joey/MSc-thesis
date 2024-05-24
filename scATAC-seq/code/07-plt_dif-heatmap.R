#args <- list(list.files("outs", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- lapply(res, \(x) x %>% filter(row.names(x) == motif)) 
df <- lapply(res, select, rank, t, motif, 
        dif, mode, moderated, smooth, peakWeight) |>
    do.call(what=rbind) |>
    mutate(method=paste(mode,moderated, smooth, peakWeight, dif, sep=","))




pr <- ggplot(df, 
    aes(reorder(motif,-sqrt(rank)), 
        reorder(method,-sqrt(rank)), fill=rank)) +
    geom_tile(col="white") +
    geom_text(aes(label = rank), size = 1.5) +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        n.breaks=3, direction=-1) +
    labs(x="motifs", y="method") +
    coord_fixed(expand=FALSE) + 
    theme_minimal(6) + theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank")

pt <- ggplot(df,
    aes(reorder(motif,-sqrt(rank)),
        reorder(method,-sqrt(rank)), fill=t)) +
    geom_tile(col="white") +
    geom_text(aes(label = round(t, 2)), size = 1.5) +
    scale_fill_distiller(NULL,
        palette="RdYlBu", na.value="lightgrey",
        n.breaks=3, direction=-1) +
    labs(x="motifs", y="method") +
    coord_fixed(expand=FALSE) +
    theme_minimal(6) + theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("t-value")

gg <- pt + pr + plot_annotation(tag_levels="a")
ggsave(args[[2]], gg, width=18, height=12, units="cm")