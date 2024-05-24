suppressPackageStartupMessages({
    library(data.table)
    source("~/chromVAR/R/getCounts.R")
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(SummarizedExperiment)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(fields)
})

bams <- list.files("/mnt/plger/datasets/CusanovichATACprotocols/aligned", 
    pattern = "\\.bam$",
    full.names = TRUE)
sample <- vapply(bams, \(x) {
    name <- basename(x)
    unlist(strsplit(name, ".", fixed = TRUE))[2]
}, character(1))
names(bams) <- sample

meta <- read.table("/mnt/plger/datasets/CusanovichATACprotocols/SraRunTable.txt",
    sep = ",", header=TRUE)

#human <- meta[meta$Organism==""]
# origin <- lapply(seq_along(bams), \(x){
#   lst <- as.list(bams[[x]])
#   names(lst) <- names(bams)[[x]]
#   res <- tryCatch ({
#     .importFragments(lst)[[1]]
#   }, error = function(e){
#     NULL
#   })
# })

# Homo sapiens
atacFrags <- readRDS("/mnt/plger/jwang/Cusanovich/origin.rds")
names(atacFrags) <- names(bams)
atacFrags <- atacFrags[!vapply(atacFrags, is.null, logical(1))]
atacFrags <- atacFrags[names(atacFrags) %in% meta[meta$Organism=="Homo sapiens", "Sample.Name"]]

peaks <- import.bed(con="/mnt/plger/datasets/CusanovichATACprotocols/peaks/SRX10861468.GSM5290888_summits.bed")
counts <- getCounts(atacFrag = atacFrags, ranges = peaks, 
                    genome = BSgenome.Hsapiens.UCSC.hg38, 
                    width = 300,
                    species = "Homo sapiens")
saveRDS(counts, "/mnt/plger/jwang/Cusanovich/counts_HS.rds")

# weights
weights <- getCounts(atacFrag = atacFrags,
  ranges = peaks,
  mode = "weight",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "Homo sapiens",
  width = 300,
  smooth = "smooth.2d",
  aRange = 2)

saveRDS(weights, "/mnt/plger/jwang/Cusanovich/weights_HS_s2.rds")

# weights with bias correction
weights <- getCounts(atacFrag = atacFrags,
  ranges = peaks,
  mode = "weight",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "Homo sapiens",
  width = 300,
  peakWeight = "loess2"
  )

saveRDS(weights, "/mnt/plger/jwang/Cusanovich/weights_HS_loess2.rds")

# Mus_musculus
atacFrags <- readRDS("/mnt/plger/jwang/Cusanovich/origin.rds")
names(atacFrags) <- names(bams)
atacFrags <- atacFrags[!vapply(atacFrags, is.null, logical(1))]
atacFrags <- atacFrags[names(atacFrags) %in% meta[meta$Organism=="Mus musculus", "Sample.Name"]]

peaks <- import.bed(con="/mnt/plger/datasets/CusanovichATACprotocols/peaks/SRX14660758.GSM5983291_summits.bed")
counts <- getCounts(atacFrag = atacFrags, ranges = peaks, 
                    genome = BSgenome.Mmusculus.UCSC.mm10, 
                    width = 300,
                    species = "Mus_musculus")
saveRDS(counts, "/mnt/plger/jwang/Cusanovich/counts_MM.rds")

weights <- getCounts(atacFrag = atacFrags,
                     ranges = peaks,
                     mode = "weight",
                     genome = BSgenome.Mmusculus.UCSC.mm39,
                     species = "Mus_musculus",
                     width = 300)
saveRDS(weights, "/mnt/plger/jwang/Cusanovich/weights_MM.rds")





