#args <- list(list.files("outs/dat", "^mmd-.*", full.names=TRUE), "plts/mmd-heatmap.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(data.table)
    library(tidyr)
    library(ggh4x)
    library(tidytext)
    library(cowplot)
})

res <- lapply(args[[1]], readRDS) 
mm <- do.call(rbind, res) |>
    mutate(method=paste(symmetric,normalization,dif,sep=",")) |>
    select(target, method, rank, t, name) |>
    group_by(method, target) |>
    filter(name %in% unlist(strsplit(target, ","))) |>
    mutate(min_rank = min(rank)) |>
    filter(rank == min_rank) |>
    ungroup() |>
    select(-min_rank)


motifs <- sapply(unique(mm$target), \(x) unlist(strsplit(x,","))[[1]])
bl <- lapply(motifs, \(m) {
    x <- readRDS(paste0("outs/dat/dif-total-",m,",chromVAR>limma.rds")) |>
        select(truth, rank, t, name) |>
        rename("target"="truth") |>
        mutate(method="chromVAR>limma") |>
        group_by(method, target) |>
        filter(name %in% unlist(strsplit(target, ","))) |>
        mutate(min_rank = min(rank)) |>
        filter(rank == min_rank) |>
        ungroup() |>
        select(-min_rank)
}) |> do.call(what=rbind)
df <- rbind(mm,bl)

pr <- ggplot(df, 
    aes(reorder(target,-sqrt(rank)), 
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
    aes(reorder(target,-sqrt(rank)), 
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

