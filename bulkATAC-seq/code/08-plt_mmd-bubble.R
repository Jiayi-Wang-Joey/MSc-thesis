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
    library(preprocessCore)
    library("viridis") 
})

res <- lapply(args[[1]], readRDS) 
mm <- do.call(rbind, res) |>
    mutate(method=paste(symmetric,normalization,dif,sep=",")) |>
    select(target, method, rank, t, name, adj.P.Val) |>
    group_by(method, target) |>
    filter(name %in% unlist(strsplit(target, ","))) |>
    mutate(min_rank = min(rank)) |>
    filter(rank == min_rank) |>
    ungroup() |>
    select(-min_rank) |>
    rename("label_rank"="rank") |>
    mutate(significance=ifelse(adj.P.Val < 0.05, "TRUE", "FALSE"),
        rank=case_when(label_rank > 100 ~ 100, TRUE ~ label_rank))
    


motifs <- sapply(unique(mm$target), \(x) unlist(strsplit(x,","))[[1]])
bl <- lapply(motifs, \(m) {
    x <- readRDS(paste0("outs/dat/dif-total-",m,",chromVAR>limma.rds")) |>
        select(truth, rank, t, name, adj.P.Val) |>
        rename("target"="truth") |>
        mutate(method="chromVAR>limma") |>
        group_by(method, target) |>
        filter(name %in% unlist(strsplit(target, ","))) |>
        mutate(min_rank = min(rank)) |>
        filter(rank == min_rank) |>
        ungroup() |>
        select(-min_rank) |>
        rename("label_rank"="rank") |>
        mutate(significance=ifelse(adj.P.Val < 0.05, "TRUE", "FALSE"),
            rank=case_when(label_rank > 100 ~ 100, TRUE ~ label_rank))
}) |> do.call(what=rbind)
df <- rbind(mm,bl) |>
    mutate(method = case_when(
        method == "sym,TMM,z>limma" ~ "sym>z-score>limma",
        method == "asy,TMM,z>limma" ~ "asy>z-score>limma",
        method == "asy,TMM,motifScore>limma_voom" ~ "asy>maxMS>limma_voom",
        method == "sym,TMM,motifScore>limma_voom" ~ "sym>maxMS>limma_voom",
        method == "asy,TMM,sumMS>limma_voom" ~ "asy>sumMS>limma_voom",
        method == "sym,TMM,sumMS>limma_voom" ~ "sym>sumMS>limma_voom",
        TRUE ~ method  
    ))


df$col1 <- ifelse(df$rank>=100,"black","white")       
df$col2 <- ifelse(abs(df$t)<15,"black","white")
aes_r <- list(geom_point(aes( 
    color=rank, size=significance)),
    scale_size_manual(values = c("TRUE"=9, "FALSE"=6.5, 1)),
    geom_text(aes(label=label_rank),size=2.5, color=df$col1),
    scale_color_viridis(discrete=FALSE, option = "D",
        breaks=c(5, 50, 100),
        labels=c("5", "50", ">= 100")),
    theme_bw(),
    labs(x="perturbed TFs", y="method"),
    theme(legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
    guides())


aes_t <- list(geom_point(aes(color=abs(t), size=significance)),
    scale_size_manual(values = c("TRUE"=9, "FALSE"=6.5, 1)),
    geom_text(aes(label=round(t,2)), size=2.2, color=df$col2),
    scale_color_viridis(discrete=FALSE, option = "magma",direction=-1),
    theme_bw(),
    labs(x="perturbed TFs", y="method"),
    theme(legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.position="bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
    guides()) 

.p <- \(df) {
    ggplot(df, aes(reorder(target,-sqrt(rank)),
        reorder(method,-sqrt(rank)))) + aes_r + ggtitle("rank") +
        ggplot(df, aes(reorder(target,-sqrt(rank)),
            reorder(method,-sqrt(rank)))) + aes_t + ggtitle("t-statistics") +
        plot_layout(nrow=1, guides="collect", axes="collect_x") &
        theme(legend.position="bottom")
}

p1 <- .p(df)

# pdf(args[[2]], onefile=TRUE, width=8, height=5)
# for (p in list(p1,p2)) print(p); dev.off()
ggsave(args[[2]], p1, width=23, height=13, units="cm")

