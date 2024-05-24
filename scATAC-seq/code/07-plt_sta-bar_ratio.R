#args <- list(list.files("outs", "^sta-.*", full.names=TRUE), "plts/sta-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(dplyr)
    library(tidytext)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res,select, motif, mode, smooth, peakWeight,
    status, sta, sta_val, motif_id) |>
    do.call(what=rbind) |>
    filter(sta=="ratio" & motif==motif_id) |>
    mutate(mode=ifelse(mode=="weight","fragWeight","none")) |>
    mutate(method=paste(mode, status, smooth, sep = ">"))


gg <- ggplot(df, aes(x=reorder_within(method, sta_val, motif), y=sta_val, fill=method)) +
    geom_bar(stat="identity") +
    facet_wrap2(~motif, scales="free", ncol=6) +
    theme_bw() +
    scale_fill_brewer(palette="Paired") +
    labs(x="method",y="ratio") +
    theme_minimal() +
    theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x=element_blank()) 
    
    
ggsave(args[[2]], gg, width=35, height=8, units="cm")