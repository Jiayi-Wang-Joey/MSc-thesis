suppressPackageStartupMessages({
    library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("code/utils.R")
    library(MotifDb)
    library(SummarizedExperiment)
    library(chromVAR)
    library(limma)
    library(edgeR)
})

m <- wcs$mot
mice <- c("BANP", "NR1H3", "NR1H4")
se <- readRDS(args$dat)
fun <- function(x, genome, motif) {
    x <- x[which(!is.infinite(rowSums(counts(x)))),]
    x <- filterPeaks(x, non_overlapping = TRUE)
    motif_ix <- matchMotifs(motif, x, genome = genome)
    y <- t(assay(motif_ix))%*%assay(x, "counts")
    group_id <- substr(colnames(y), 1, nchar(colnames(y)) - 5)
    design <- model.matrix(~ group_id)
    y <- calcNormFactors(DGEList(y))
    y <- voom(y,design)
    fit <- eBayes(lmFit(y, design))
    res <- topTable(fit, n = Inf)
    res <- res[order(res$P.Value, -abs(res$logFC)),]
    res$rank <- seq_len(nrow(res))
    res$name <- rownames(res)
    res 
  
}

if (m %in% mice) {
    Mmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Mmusculus")
    if (m=="BANP") {
        banp <- readRDS("../chromVAR2/data/BANP.PFMatrix.rds")
        Mmotifs$BANP <- banp
    } else if (m=="NR1H3") {
        Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
        Mmotifs$NR1H3 <- Hmotifs$NR1H3
    }
    res <- fun(se, 
        genome = BSgenome.Mmusculus.UCSC.mm10, 
        motif = Mmotifs)
} else {
    Hmotifs <- getNonRedundantMotifs(format="PFMatrix", species="Hsapiens")
    res <- fun(se, 
        genome = BSgenome.Hsapiens.UCSC.hg38, 
        motif = Hmotifs)
}

df <- data.frame(res, motif=m, 
  seed=wcs$seed, span=wcs$span, 
  nSample=wcs$ns,family=wcs$family,
  nBin=wcs$nb, row.names=NULL)

saveRDS(df, args$res)