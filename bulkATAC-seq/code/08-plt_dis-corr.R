#args <- list(list.files("/mnt/plger/jwang/data/dat/00-dist", full.names=TRUE), "plts/dif-rankHM.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(data.table)
    library(tidyr)
    library(ggrastr)
})

res <- lapply(args[[1]], readRDS)
res <- lapply(res,\(dt) {
    dt[,freq_bin:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
    dt[,global_freq_bin:=(count_bin+1L)/sum(count_bin+1L)]
    dt[,sample_id:=as.integer(as.factor(sample))]
    dt
})
df <- lapply(res, select, motif, sample_id, weight, moderated_weight, bin) |>
    do.call(what=rbind)


gg <- ggplot(df, aes(weight,moderated_weight,col=factor(sample_id))) +
    facet_wrap(~motif, ncol=4, scales = "free") +
    geom_point_rast(alpha=0.3) +
    scale_color_brewer(palette = "Paired") +
    theme_bw()


ggsave(args[[2]], gg, width=30, height=15, units="cm")
