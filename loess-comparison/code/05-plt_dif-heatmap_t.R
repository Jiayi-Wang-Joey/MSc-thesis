#args <- list(list.files("outs", "^dif-", full.names=TRUE), "plts/dif-stability.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
})

res <- lapply(args[[1]], readRDS)
res <- res[!vapply(res, is.null, logical(1))]

df <- lapply(res, select,
    logFC, adj.P.Val, t, span, nSample, 
    nBin, family, seed, motif, rank, name) |>
    do.call(what=rbind) |>
    filter(motif==name) |>
    group_by(motif, span, nSample, nBin, family) |>
    slice(which.max(abs(t))) |>
    ungroup()


df$parameter <- paste0("span_", df$span, ",family_", df$family)
df$sampleSize <- paste0("nSample_", df$nSample, ",nBin_", df$nBin)

# aesthetics
aes <- list(
    # scale_fill_distiller(NULL,
    #     palette="RdYlBu", na.value="lightgrey",
    #     n.breaks=5, direction=-1),
    scale_fill_gradientn( 
        #colors=c("moccasin", "pink", "red", "firebrick", "firebrick4","black")
        colors=c("lightsalmon", 
            "wheat1","moccasin","paleturquoise3"),
    ),
    labs(x="motif", y="parameters"),
    theme_minimal(), 
    theme(legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))


p1 <- ggplot(df, aes(reorder(motif, abs(t)), reorder(parameter, abs(t)),
        fill=abs(t))) +
        geom_tile(col="white") +
        geom_text(aes(label = round(t, 2)), size=2.5) +
        aes +
        facet_grid(nSample~nBin, scales="free") +
        ggtitle("nSample~nBin")

p2 <- ggplot(df, aes(reorder(motif, abs(t)), reorder(sampleSize, abs(t)),
    fill=abs(t))) +
    geom_tile(col="white") +
    geom_text(aes(label = round(t, 2)), size=2.5) +
    aes +
    facet_grid(span~family, scales="free") +
    ggtitle("span~family")


pdf(args[[2]], onefile=TRUE, width=8, height=10)
for (p in list(p1,p2)) print(p); dev.off()







