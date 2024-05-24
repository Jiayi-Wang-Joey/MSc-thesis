suppressPackageStartupMessages({
    library(edgeR)
    library(SummarizedExperiment)
})
source(args$fun)
x <- readRDS(args$dat)
scores <- x[[2]]
res <- fun(scores, wcs$nrm)
m <- wcs$mot
if (m=="MYC") { 
    target <- paste("MYC", "MAX", sep=",")
} else if (m=="NR1H4"){
    target <- paste("NR1H4","RXRA","RXRB",sep=",")
    #target <- "NR1H4"
} else if (m=="ESR1"){
    target <- paste("ESR1", "ESR2", sep=",")
} else {
    target <- m
}

df <- data.frame(res, target=target, 
    symmetric=wcs$sym, normalization=wcs$nrm,
    dif=wcs$mmd)

saveRDS(df, args$res)