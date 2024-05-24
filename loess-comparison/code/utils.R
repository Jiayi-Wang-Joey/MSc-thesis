.weightPeaks <- function(counts, 
    nBins,
    nSample,
    span=0.3,
    family="symmetric",
    ...) {
    #nBins <- nBins
    #nSample <- nSample
    allAvg <- log2(rowMeans(counts+1L))
    if (nrow(counts)<=1e4) {
        idx <- 1:nrow(counts) 
    } else {
    bins <- cut(allAvg,
        breaks = seq(min(allAvg), max(allAvg), (max(allAvg)-min(allAvg))/nBins),
        include.lowest = TRUE)
    idx <- unlist(sapply(seq_len(length(levels(bins))), \(i) {
        ids <- which(bins==levels(bins)[i])
        if (length(ids) > 1)
            sample(ids,
                size = min(length(ids), round(nSample/nBins)))
        else if (length(ids)==1)
            ids
        else NULL
    }))}
    res <- affy::normalize.loess(counts+1L, span=span, 
        subset=idx, family.loess=family)
    
    return(res)
}


getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
  species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
    TF=gsub(pat,"",names(motifs)),
    grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
    universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
    PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
      universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
      modf$TF)),
    PWMatrix=do.call(TFBSTools::PWMatrixList, 
      setNames(universalmotif::convert_motifs(motifs, 
        class="TFBSTools-PWMatrix"), modf$TF))
  )
}