# args <- list(list.files("outs", "^eva-.*", full.names=TRUE), "plts/eva-bar.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggh4x)
    library(tidytext)
})
res <- lapply(args[[1]], readRDS)
df <- do.call(rbind, res) |> 
    mutate(method=paste(mode,status,smooth,peakWeight, sep=","))
    
gg <- ggplot(df, aes(x=reorder_within(method, eva_val, eva), 
    y=eva_val, fill=method)) +
    geom_bar(stat="identity") +
    facet_grid2(~eva, scales="free", independent="all") +
    theme_bw() +
    scale_fill_brewer(palette="Paired") +
    labs(x="method",y="Value") +
    theme_minimal(6) +
    theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x=element_blank()) 

ggsave(args[[2]], gg, width=20, height=8, units="cm")