suppressPackageStartupMessages({
    library(SummarizedExperiment)
    source("code/utils.R")
})

set.seed(seed <- as.numeric(wcs$seed))
print(seed)
motif <- wcs$mot
se <- readRDS(args$dat)
counts <- assay(se, "counts")

newCounts <- .weightPeaks(counts, 
    nBins = as.numeric(wcs$nb), family=wcs$family,
    nSample = as.numeric(wcs$ns), span=as.numeric(wcs$span), 
    degree=1)

assay(se,"counts") <- newCounts

saveRDS(se, args[[2]])
