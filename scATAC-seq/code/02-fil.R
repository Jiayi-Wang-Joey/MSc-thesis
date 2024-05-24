suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(chromVAR)
    library(Seurat)
    library(Signac)
})

se <- readRDS(args$dat)
colData(se)$depth <- colSums(counts(se))
se <- filterSamples(se, min_depth = 2000,
   min_in_peaks = 0.15, shiny = FALSE)
srt <- CreateSeuratObject(
    counts = counts(se),
    assay = "peaks",
    meta.data = as.data.frame(colData(se))
)

srt <- RunTFIDF(srt)
srt <- FindTopFeatures(srt, min.cutoff = 'q0')
srt <- RunSVD(srt)
srt <- RunUMAP(object = srt, reduction = 'lsi', dims = 2:30)

umap <- srt@reductions$umap@cell.embeddings
lsi <- srt@reductions$lsi@cell.embeddings

if (is.null(wcs$smt)) {
    cd <- data.frame(colData(se), umap,
        lsi, motif=wcs$mot,
        mode="total",status="none",
        smooth="none", peakWeight="none")
    colData(se) <- DataFrame(cd)
} else {
    cd <- data.frame(colData(se), umap,
        lsi, motif=wcs$mot,
        mode="weight",status=wcs$mdr,
        smooth=wcs$smt, peakWeight=wcs$pkw)
    colData(se) <- DataFrame(cd)
}

saveRDS(se, args$res)
saveRDS(cd, args$cd)


#' TODO: do I need these?
# compute nucleosome signal score per cell
#srt <- NucleosomeSignal(object = srt)

# compute TSS enrichment score per cell
#srt <- TSsenrichment(object = srt, fast = FALSE)


 
