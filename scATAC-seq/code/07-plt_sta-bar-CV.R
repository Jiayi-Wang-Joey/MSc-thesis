#args <- list(list.files("outs", "^sta-.*", full.names=TRUE), "plts/sta-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(dplyr)
    library(tidytext)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res,select, motif, mode, smooth, sample_id,
    status, sta, sta_val, motif_id, peakWeight) |>
    do.call(what=rbind) |>
    mutate(mode=ifelse(mode=="weight","fragWeight","none")) |>
    mutate(method=paste(mode, status, smooth, sep = ">")) |>
    filter(sta=="CV") |>
    mutate(group_id=gsub("^sg|_[0-9]$", "", sample_id)) |>
    mutate(group_id=ifelse(group_id=="NT", "CTRL", "TREAT")) |>
    group_by(group_id, motif, method) |>
    summarise_at("sta_val", mean) 



gg <- ggplot(df, aes(x=reorder_within(method, sta_val, list(motif, group_id)), 
    y=sta_val, fill=method)) +
    geom_bar(stat="identity") +
    facet_grid2(group_id~motif, scales="free", independent="all") +
    scale_fill_brewer(palette="Paired") +
    labs(x="method",y="Coefficient of variation") +
    theme_minimal() +
    theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x=element_blank()) 

ggsave(args[[2]], gg, width=35, height=12, units="cm")