#args <- list(list.files("outs/sim", "^dif-.*", full.names=TRUE), "plts/sim/dif-heatmap.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(dplyr)
})



# res <- lapply(args[[1]], function(x) {
#     df <- readRDS(x)
#     test <- ifelse(df$test[1]!="ZNF143", df$test[1], "ZN143")
#     df <- df[row.names(df)==test,]
#     df <- df[,c("t", "rank", "mode", "dif", "effect", 
#         "test", "smooth", "peakWeight")]
# })

#res <- res[vapply(res, function(x) nrow(x)!=0, logical(1))]

res <- lapply(args[[1]], readRDS)
#res <- res[vapply(res, \(x) nrow(x)!=0, logical(1))]

df <- lapply(res, select, t, rank, test, effect, 
    name, dif, mode, smooth, peakWeight) |>
    do.call(what=rbind) |>
    #filter(name==test) |>
    mutate(test=ifelse(test!="ZNF143", test, "ZN143"),
        method=paste(mode,smooth,peakWeight,dif,sep=",")) |>
    filter(name==test)


# aesthetics
gg <- list(
    scale_fill_distiller(NULL,
            palette="RdYlBu", na.value="lightgrey",
            n.breaks=3, direction=-1),
    labs(x="effect", y="method"),
    coord_fixed(expand=FALSE),
    theme_minimal(6), 
    theme(legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines")))

df$effect <- factor(df$effect, 
    levels=c("0", "0.25", "0.5", "1", "3"))

ps <- lapply(split(df, df$test), \(fd) {
    p1 <- ggplot(fd, aes(effect, reorder(method,-sqrt(rank)), fill=rank)) +
        geom_tile(col="white") +
        geom_text(aes(label = rank), size = 1.5) +
        gg +
        ggtitle("rank")
    
    p2 <- ggplot(fd, aes(effect, reorder(method,abs(t)), fill=t)) +
        geom_tile(col="white") +
        geom_text(aes(label = round(t, 2)), size = 1.5) +
        gg +
        ggtitle("t")
    
    p2 + p1 + plot_annotation(title = fd$test[1], tag_levels="a") 
    
})

pdf(args[[2]], onefile=TRUE, width=5, height=5)
for (p in ps) print(p); dev.off()

