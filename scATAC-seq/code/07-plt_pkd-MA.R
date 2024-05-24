#args <- list(list.files("outs", "^pkd-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    library(patchwork)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, logFC, AveExpr,
    motif, mode, smooth, status, dif, peakWeight) |> 
    do.call(what=rbind) |>
    filter(dif=="limma_voom")

gg <- list(
    geom_point_rast(alpha=0.3),
    geom_density2d(alpha=0.5),
    geom_smooth(alpha=0.1, colour="darkred", linewidth=0.5),
    theme(plot.title = element_text(size = 7)),
    theme_bw()
)

ps <- lapply(split(df, df$motif), \(fd) {
    ggplot(fd, aes(x=AveExpr, y=logFC)) +
        gg +
        facet_grid(status~smooth+mode+peakWeight, scales="free", 
            labeller = \(.) label_both(.)) +
        ggtitle(fd$motif[1])
})


pdf(args[[2]], onefile=TRUE, width=12, height=7)
for (p in ps) print(p); dev.off()