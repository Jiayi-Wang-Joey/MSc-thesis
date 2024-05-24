#args <- list(list.files("outs", "^sta-.*", full.names=TRUE), "plts/sta-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(dplyr)
})

res <- lapply(args[[1]], readRDS)
df <- lapply(res,select, motif, mode, smooth, sample_id,
    status, sta, sta_val, motif_id, peakWeight) |>
    do.call(what=rbind) |>
    mutate(method=paste(mode, status, smooth, peakWeight, sep = ",")) |>
    filter(sta=="within-sample_variance") |>
    mutate(group_id=gsub("^sg|_[0-9]$", "", sample_id)) |>
    mutate(group_id=ifelse(group_id=="NT", "CTRL", "TREAT"))


aes <- list(xlab("within-sample variances (log10)"),
    geom_density(),
    scale_color_brewer(palette = "Set1"),
    theme_minimal(),
    theme(
        legend.position="bottom",
        panel.border=element_rect(fill=NA),
        legend.key.width=unit(1, "lines"),
        legend.key.height=unit(0.5, "lines")))

# p1 <- ggplot(df, aes(x=log(sta_val), col=group_id)) +
#     facet_grid(method~motif, scales="free") +
#     aes

p2 <- ggplot(df, aes(x=log10(sta_val), col=method)) +
    facet_grid(group_id~motif, scales="free") +
    aes

ggsave(args[[2]], p2, width=35, height=13, units="cm")

# pdf(args[[2]], onefile=TRUE, width=20, height=12)
# for (p in list(p1,p2)) print(p); dev.off()