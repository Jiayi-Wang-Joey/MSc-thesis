#args <- list(list.files("outs/sim", "^pkd_", full.names=TRUE), "plts/sim/pkd-MA.R")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    library(patchwork) 
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind, res) |>
    filter(effect=="GC") |>
    rename("fragWeight"="mode") |>
    mutate(fragWeight=ifelse(fragWeight=="fragWeight", "Yes", "No")) |>
    mutate(peakWeight=ifelse(peakWeight=="loess", "Yes", "No"))

aes <- list(
  geom_point_rast(alpha=0.3),
  geom_density2d(alpha=0.8),
  geom_smooth(alpha=0.1, colour="darkred", linewidth=1),
  theme(plot.title = element_text(size = 7)),
  theme_bw(15),
  xlab("mean accessibility")
)

gg <-   ggplot(df, aes(x=AveExpr, y=logFC)) +
    aes +
    facet_grid(motif~peakWeight+fragWeight, scales="free",
        labeller=\(.) label_both(.))


ggsave(args[[2]], gg, width=25, height=18, units="cm")