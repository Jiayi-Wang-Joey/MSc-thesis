suppressPackageStartupMessages({
    library(BSgenome.Mmusculus.UCSC.mm10)
    source("~/chromVAR/R/getCounts.R")
    library(data.table)
    library(SummarizedExperiment)
    library(fields)
    library(SingleCellExperiment)
    library(chromVAR)
    library(Seurat)
    library(Signac)
})

smt <- strsplit(wcs$smt, "_")[[1]]
mdr <- ifelse(wcs$mdr=="unmoderated", FALSE, TRUE)
if (smt[1]=="none") aRange <- 0 else aRange <- as.numeric(smt[2])

atacFrag <- readRDS("/mnt/plger/jwang/scATAC-seq/fragDt_waag.rds")
peaks <- import.bed(con="/mnt/germain/datasets/2022_waag_inVivo_scMultiome/o294832_sample_1/outs/atac_peaks.bed")

se <- getCounts(atacFrag = atacFrag,
    ranges = peaks,
    mode = "weight",
    genome = BSgenome.Mmusculus.UCSC.mm10,
    species = "Mus_musculus",
    width = 300,
    smooth = smt[1],
    nGCBins = 10,
    nWidthBins = 30,
    aRange = aRange,
    moderating = mdr, 
    peakWeight = wcs$pkw,
    singleCell = TRUE)

colData(se)$depth <- colSums(counts(se))
se <- filterPeaks(se, non_overlapping = TRUE)
se <- filterSamples(se, min_depth = 2000,
    min_in_peaks = 0.15, shiny = FALSE)
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
srt <- FindNeighbors(object = srt, reduction = 'lsi', dims = 2:30)
srt <- FindClusters(object = srt, verbose = FALSE, algorithm = 3, resolution=0.4)
umap <- srt@reductions$umap@cell.embeddings
lsi <- srt@reductions$lsi@cell.embeddings
cd <- data.frame(colData(se), umap, lsi, 
    cluster_id = srt$peaks_snn_res.0.4,
    mode="weight",status=wcs$mdr,
    smooth=wcs$smt, peakWeight=wcs$pkw)
idx <- match(cd$barcode, rownames(sce))
cd$cell_type <- sce$cell_type_label_joint[idx]

colData(se) <- DataFrame(cd)


saveRDS(cd, args$cd)
saveRDS(se, args$res)

