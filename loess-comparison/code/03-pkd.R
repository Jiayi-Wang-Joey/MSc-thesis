suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(BSgenome.Hsapiens.UCSC.hg38)
    source("~/chromVAR/R/getCounts.R")
    library(edgeR)
    library(data.table)
})

.limma_voom <- function(x) {
  counts <- assay(x, "counts")
  group_id <- rep(LETTERS[1:2],each=ncol(counts)/2)
  design <- model.matrix(~ group_id)
  rownames(counts) <- seq_len(nrow(counts))
  y <- calcNormFactors(DGEList(counts))
  y <- voom(y,design)
  fit <- eBayes(lmFit(y, design))
  res <- topTable(fit, n = Inf, sort.by = "none")
  res 
}

se <- readRDS(args$dat)
peaks <- .getGCContent(rowRanges(se), genome=BSgenome.Hsapiens.UCSC.hg38)
peaks <- as.data.table(peaks)
mean_width <- rowMeans(assay(se,"mean_width"))
median_width <- rowMeans(assay(se,"median_width"))

res <- .limma_voom(se)
if (!is.null(wcs$pkw)) {
  df <- data.frame(res, mode="fragWeight",
        motif=wcs$mtf, peakWeight=wcs$pkw, 
        effect=wcs$eft, dsg=wcs$dsg, 
        gc=peaks$gc,mean_width=mean_width, 
      median_width=median_width)
} else {
  df <- data.frame(res, mode="origin",
        motif=wcs$mtf, peakWeight="none",
        effect=wcs$eft, dsg=wcs$dsg,
        gc=peaks$gc, mean_width=mean_width, 
        median_width=median_width)
}



saveRDS(df, args$res)