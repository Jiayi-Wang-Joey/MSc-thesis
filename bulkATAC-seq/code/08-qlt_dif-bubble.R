#args <- list(list.files("outs/sim", "^dif-.*", full.names=TRUE), "plts/sim/dif-bubble.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, t, rank, mode, smooth, 
    dif, effect, test, name, peakWeight, adj.P.Val) |>
    do.call(what=rbind) |>
    mutate(test=ifelse(test!="ZNF143", test, "ZN143"),
        significance=ifelse(adj.P.Val < 0.05, TRUE, FALSE),
        method=paste(mode,smooth,peakWeight,dif,sep=",")) |>
    filter(name==test)



ps <- lapply(split(df, df$test), \(fd) {
    rng <- range(fd$rank, na.rm=TRUE)
    rng <- c(
        floor(rng[1]*10)/10, 
        ceiling(rng[2]*10)/10)
    ggplot(fd, aes(x = effect, y = reorder(method, significance), 
        size = significance, color = rank)) +
        geom_point() +
        scale_color_gradientn(
            limits=rng, breaks=rng, na.value="lightgrey", 
            colors=c("black", "firebrick", "red", "pink", "moccasin")) +
        theme_minimal() + 
        ggtitle(fd$test[1]) +
        ylab("methods")
}) |> wrap_plots(ncol = 2)

ggsave(args[[2]], ps, width=30, height=20, units="cm")
# pdf(args[[2]], onefile=TRUE, width=6, height=4)
# for (p in ps) print(p); dev.off()