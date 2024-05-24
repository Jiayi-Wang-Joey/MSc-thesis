#args <- list(list.files("outs/sim", "^pkd_", full.names=TRUE), "plts/sim/pkd-MA.R")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    library(patchwork)
    library(data.table)
    library(RColorBrewer)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- do.call(rbind, res) |>
    filter(effect=="FLD") |>
    rename("fragWeight"="mode") |>
    mutate(fragWeight=ifelse(fragWeight=="fragWeight", "Yes", "No")) |>
    mutate(peakWeight=ifelse(peakWeight=="loess", "Yes", "No"))
df <- data.table(df)

FLIntervals <-  unique(quantile(df$median_width,
     probs = seq(0,1,by=1/10)))
#FLIntervals <- c(min(df$median_width), 90, 120, 300, 500, max(df$median_width))
df[,FLBin:=cut(median_width, 
    breaks=FLIntervals, 
    include.lowest=TRUE)]

aes <- list(
    geom_violin_rast(), 
    geom_boxplot(width=0.1),
    theme_bw(15),
    scale_color_manual(values=colorRampPalette(brewer.pal(9, "RdYlBu"))(10)),
    geom_abline(intercept=0, slope=0, col="black", lty=2),
    geom_smooth(method="gam",se=FALSE, color="darkred", aes(group=1), lwd=0.8, alpha=0.7),
    guides(col = FALSE),
    theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)

gg <- ggplot(df, aes(x=FLBin, y=logFC, col=FLBin)) +
    aes +
    facet_grid(motif~peakWeight+fragWeight, scales="free",
        labeller=\(.) label_both(.)) +
    xlab("FL bin")

ggsave(args[[2]], gg, width=25, height=18, units="cm")