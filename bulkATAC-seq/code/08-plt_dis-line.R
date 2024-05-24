#args <- list(list.files("/mnt/plger/jwang/data/dat/00-dist", full.names=TRUE), "plts/dif-rankHM.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(data.table)
    library(tidyr)
})

res <- lapply(args[[1]], readRDS)
res <- lapply(res,\(dt) {
  dt[,freq_bin:=(count_bin+1L)/(sum(count_bin)+1L),by=sample]
  dt[,global_freq_bin:=(count_bin+1L)/sum(count_bin+1L)]
  dt[,sample_id:=as.integer(as.factor(sample))]
  dt
  })
df <- lapply(res, select, motif, sample_id, mean_freq, moderated_mean_freq,
    uncorrected_Freq, uncorrected_moderatedFreq, 
    weight, moderated_weight, bin, count_bin) |>
  do.call(what=rbind)

df[,corrected_Freq:=(weight*count_bin+1L)/(sum(weight*count_bin)+1L),
    by=c("motif", "sample_id")]
df[,corrected_moderatedFreq:=
        (moderated_weight*count_bin+1L)/(sum(moderated_weight*count_bin)+1L),
    by=c("motif", "sample_id")]

df <- df %>% 
    select(motif, sample_id, 
        uncorrected_Freq, uncorrected_moderatedFreq,
        corrected_Freq, corrected_moderatedFreq, bin) %>%
    group_by(motif, sample_id, bin) %>%
    pivot_longer(
        cols = starts_with(c("corrected", "uncorrected")), 
        names_to = c(".value", "type"), 
        names_sep = "_", 
        names_ptypes = list(value = numeric(), type = character())) %>%
    pivot_longer(
        cols = c(corrected, uncorrected),
        names_to = "status",
        values_to = "value")
df$status <- factor(df$status, levels=c("uncorrected", "corrected"))

aes <- list(geom_line(alpha=0.6),
  theme_minimal(6),
  scale_color_brewer(palette = "Paired"))

ps <- lapply(split(df, df$motif), \(fd) {
    ggplot(fd, aes(x=bin, y=value, col=factor(sample_id))) + aes +
        facet_grid(status~type, scales="free", labeller=\(.) label_both(.)) +
        ggtitle(fd$motif[1]) 
})



pdf(args[[2]], onefile=TRUE, width=12, height=7)
for (p in ps) print(p); dev.off()


