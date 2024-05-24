#args <- list(list.files("outs/dat", "^dif-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(tibble)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, rank, t, mode, moderated, 
    truth, name, smooth, peakWeight, dif, adj.P.Val) |>
    do.call(what=rbind) |>
    mutate(method=paste(mode,moderated,smooth,peakWeight,dif,sep=",")) |>
    group_by(method, truth) |>
    filter(name %in% unlist(strsplit(truth, ","))) |>
    mutate(min_rank = min(rank)) |>
    filter(rank == min_rank) |>
    ungroup() |>
    select(-min_rank) |>
    rename("label_rank"="rank") |>
    mutate(significance=ifelse(adj.P.Val < 0.05, "TRUE", "FALSE"),
        mode=ifelse(mode=="weight","fragWeight","none"),
        peakWeight=ifelse(peakWeight=="loess","peakWeight","none"),
        rank=case_when(label_rank > 100 ~ 100, TRUE ~ label_rank),
        method=paste(mode,smooth,moderated,peakWeight,dif,sep=">")) 


aes_r <- list(
    geom_tile(col="white"),
    geom_text(aes(label=label_rank),size=2.5, color="white"),
    scale_fill_gradientn(
        #colors=c("navajowhite1", "gold", "red", "navy"),
        colors=c("tomato", "lightsalmon", "wheat1","paleturquoise2","paleturquoise3"),
        #colors=c("moccasin", "pink", "red", "firebrick", "firebrick4","black"),
        breaks=c(5, 50, 100),
        labels=c("5", "50", ">= 100")
    ),
    theme_bw(),
    labs(x="perturbed TFs", y="method"),
    theme(legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
    guides())

aes_t <- list(
    geom_tile(col="white"),
    geom_text(aes(label=round(t,2)), size=2.2, color="white"),
    scale_fill_gradientn( 
        #colors=c("moccasin", "pink", "red", "firebrick", "firebrick4","black")
        colors=c("tomato","lightsalmon", "wheat1","paleturquoise2","paleturquoise3"),
    ),
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
    ggplot(df, aes(reorder(truth,-sqrt(rank)),
        reorder(method,-sqrt(rank)), fill=rank)) + aes_r + ggtitle("rank") +
        ggplot(df, aes(reorder(truth,-sqrt(rank)),
            reorder(method,-sqrt(rank)), fill=abs(t))) + aes_t + ggtitle("t-statistics") +
        plot_layout(nrow=1, guides="collect", axes="collect_x") &
        theme(legend.position="bottom")
}

p1 <- .p(df)



ggsave(args[[2]], p1, width=25, height=15, units="cm")

