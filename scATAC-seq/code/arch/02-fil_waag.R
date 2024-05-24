suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(chromVAR)
    library(Seurat)
    library(Signac)
})

se <- readRDS(args$dat)
sce <- readRDS("/mnt/plger/jwang/scATAC-seq/cd.rds")
srt <- CreateSeuratObject(
    counts = assay(se, "counts"),
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

idx <- match(cd$barcode, rownames(sce))
cd$cell_type <- sce$cell_type_label_joint[idx]

saveRDS(se, args$res)
saveRDS(cd, args$cd)