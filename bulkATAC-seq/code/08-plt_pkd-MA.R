#args <- list(list.files("outs/dat", "^pkd-.*", full.names=TRUE), "plts/dif-rankHM.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(patchwork)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res, select, logFC, AveExpr,
    test, mode, smooth, peakWeight) |> 
    do.call(what=rbind) |>
    rename("fragWeight"="mode",
        "motif"="test") |>
    mutate(
    fragWeight=ifelse(fragWeight=="weight","Yes","No"),
    peakWeight=ifelse(peakWeight=="loess","Yes","No")) 

gg <- list(
  geom_abline(intercept=0, slope=0, col="black", lty=2, alpha=0.5),
  geom_point_rast(alpha=0.3, size=0.3),
  geom_density2d(alpha=0.5),
  geom_smooth(alpha=0.1, colour="darkred", linewidth=0.7),
  theme(plot.title = element_text(size = 7)),
  theme_bw()
)


ps <- lapply(split(df, df$motif), \(fd) {
    ggplot(fd, aes(x=AveExpr, y=logFC)) +
      gg +
      facet_grid(.~fragWeight+peakWeight, scales="free",
        labeller = \(.) label_both(.)) +
      ggtitle(fd$test[1])
})

gg <- ggplot(df, aes(x=AveExpr, y=logFC)) +
      gg +
      facet_grid(motif~fragWeight+peakWeight, scales="free",
        labeller = \(.) label_both(.))

#ggsave(args[[2]], gg, width=26, height=15, units="cm")
pdf(args[[2]], onefile=TRUE, width=14, height=5)
for (p in ps) print(p); dev.off()