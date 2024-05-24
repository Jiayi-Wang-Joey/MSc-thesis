#args <- list(list.files("outs/dat", "^dif-.*", full.names=TRUE), "plts/dif-buuble_pairwsie.pdf")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
    library("viridis")
})


res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, rank, t, mode, moderated, 
    truth, name, smooth, peakWeight, dif, adj.P.Val) |>
    do.call(what=rbind) |>
    mutate(method=paste(mode,moderated,smooth,peakWeight,dif,sep=",")) |>
    group_by(method, truth) |>
    filter(name %in% unlist(strsplit(truth, ","))) |>
    #filter(name==unlist(strsplit(truth, ","))[[1]]) |>
    mutate(min_rank = min(rank)) |>
    filter(rank == min_rank) |>
    ungroup() |>
    select(-min_rank) |>
    rename("label_rank"="rank") |>
    mutate(significance=ifelse(adj.P.Val < 0.05, "TRUE", "FALSE"),
        mode=ifelse(mode=="weight","fragWeight","none"),
        peakWeight=ifelse(peakWeight=="loess","peakWeight","none"),
        rank=case_when(label_rank > 100 ~ 100, TRUE ~ label_rank))
df$col1 <- ifelse(df$rank>=100,"black","white")       
df$col2 <- ifelse(abs(df$t)<15,"black","white")
# test fragWeight
df1 <- df %>% filter(peakWeight=="none", smooth=="none",
    moderated=="unmoderated", !grepl("nf", dif)) %>%
  select(t, adj.P.Val, label_rank, truth, dif, mode, significance, rank, col1, col2) %>%
  mutate(method = paste(mode, dif, sep = ">"))

# fragWeight + smooth
df2 <- df %>% filter(moderated=="unmoderated",
  !grepl("nf", dif), peakWeight=="none") %>%
  select(t, adj.P.Val, label_rank, truth, dif, mode, significance, rank, smooth, col1, col2) %>%
  mutate(method = paste(mode, smooth, dif, sep = ">"))
# 
df3 <-  df %>% filter(
    !grepl("nf", dif), peakWeight=="none") %>%
    select(t, adj.P.Val, label_rank, truth, dif, 
        mode, significance, rank, smooth, moderated, col1,col2) %>%
    mutate(method = paste(mode, smooth, moderated, dif, sep = ">"))
# 
df4 <- df %>% filter(moderated=="unmoderated",
    !grepl("nf", dif)) %>%
    select(t, adj.P.Val, label_rank, truth, dif, 
        mode, significance, rank, smooth, peakWeight, col1,col2) %>%
    mutate(method = paste(mode, smooth, peakWeight, dif, sep = ">"))

aes_r <- list(geom_point(aes( 
    color=rank, size=significance)),
    scale_size_manual(values = c("TRUE"=9, "FALSE"=6.5, 1)),
    #geom_text(aes(label=label_rank),size=2.5, color="white"),
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

# dfs <- list(df1, df2, df3, df4)
# pr <- lapply(dfs, \(df) {
#   ggplot(df, aes(reorder(truth,-sqrt(label_rank)), 
#       reorder(method,-sqrt(label_rank)))) + aes_r 
# }) |> wrap_plots(ncol=2) & plot_annotation(title="rank", tag_levels="a") & 
#   theme(plot.tag=element_text(size=9, face="bold")) 
# 
# pt <- lapply(dfs, \(df) {
#   ggplot(df, aes(reorder(truth,-sqrt(label_rank)), 
#       reorder(method,-sqrt(label_rank)))) + aes_t 
# }) |> wrap_plots(ncol=2) & plot_annotation(title="t-value", tag_levels="a") & 
#   theme(plot.tag=element_text(size=9, face="bold")) 

.p <- \(df) {
    ggplot(df, aes(reorder(truth,-sqrt(rank)),
        reorder(method,-sqrt(rank)))) + aes_r + 
        geom_text(aes(label=label_rank),size=2.5, color=df$col1)+
        ggtitle("rank") +
        ggplot(df, aes(reorder(truth,-sqrt(rank)),
            reorder(method,-sqrt(rank)))) + aes_t +
        geom_text(aes(label=round(t,2)), size=2.2, color=df$col2) +
        ggtitle("t-statistics") +
        plot_layout(nrow=1, guides="collect", axes="collect_x") &
        theme(legend.position="bottom")
}

p1 <- .p(df1)

p2 <- .p(df2)


p3 <- .p(df3)

p4 <- .p(df4)

# pdf(args[[2]], onefile=TRUE, width=8, height=5)
# for (p in list(p1,p2)) print(p); dev.off()
ggsave(args[[2]], p4, width=25, height=12, units="cm")
