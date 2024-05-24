#args <- list(list.files("outs/sim", "^pkd-.*", full.names=TRUE), "plts/pkd-MA.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrastr)
    library(ggpointdensity)
    library(patchwork)
    library(dplyr)
})

res <- lapply(args[[1]], \(x) {
  df <- readRDS(x)
  if (!is.null(df) & df$dif[1]=="limma_voom") df
})
res <- res[sapply(res, \(x) !is.null(x))]
df <- lapply(res, select, logFC, AveExpr, effect,
  test, mode, smooth, peakWeight) |> 
  do.call(what=rbind)



gg <- list(
    geom_point_rast(alpha=0.5),
    geom_smooth(method="lm"),
    theme(plot.title = element_text(size = 7))
)

# df1 <- df[df$mode=="weight", ]
# df2 <- df[df$mode=="total",]
ps <- lapply(split(df, df$test), \(fd) {
    #dd <- df2[df2$test==fd$test[1],]
    p1 <- ggplot(fd, aes(x=AveExpr, y=logFC)) +
        gg +
        facet_grid(effect~peakWeight+smooth+mode, scales="free") 
    
    # p2 <- ggplot(dd, aes(x=baseMeanLog2, y=log2FoldChange)) + 
    #     gg + ggtitle("original counts")
    # list(p1,p2) |> wrap_plots(ncol=1) + plot_annotation(fd$test[1]) +
    #     plot_layout(heights=c(2,1))
})


pdf(args[[2]], onefile=TRUE, width=13, height=8)
for (p in ps) print(p); dev.off()