#All the functions made for analysis

#Function to read and generate SCE objects
xread <- function(sce, x) {
  sce <- read10xCounts(sce, col.names = TRUE)
  rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$ID, rowData(sce)$Symbol)
  sce$T2D <- x
  return (sce)
}
# QC Function
dataqc <- function(sce) {
  location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID,
                     column="SEQNAME", keytype="GENEID")
  
  #remove empty droplets
  set.seed(100)
  e.out <- emptyDrops(counts(sce))
  sce <- sce[,which(e.out$FDR <= 0.001)]
  unfiltered <- sce
  #remove mt-contaminated cells
  stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  sce <- sce[,!high.mito]
  colData(unfiltered) <- cbind(colData(unfiltered), stats)
  unfiltered$discard <- high.mito
  return (sce) # returns QC'd SCE object
}

#Clustering and Normalization Function
datanorm <- function(sce) {
  set.seed(1000)
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  #sce <- computeSpikeFactors(sce, character(0))
  sce <- logNormCounts(sce)
  return (sce)
}
#Graph cluster (not used)
clusterplot <- function(sce,dec){
  top.hvgs <- getTopHVGs(dec, prop=0.1)
  set.seed(10000)
  sce <- denoisePCA(sce, subset.row=top.hvgs, technical=dec)
  
  set.seed(100000)
  sce <- runTSNE(sce, dimred="PCA")
  g <- buildSNNGraph(sce, k=10, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(sce) <- factor(clust)
  #plotTSNE(sce, colour_by="label")
  return (sce)
}
#All of the above (Not used)
biocon <- function(sce){
  sce <- xread(sce)
  sce <- dataqc(sce)
  sce <- datanorm(sce)
  dec <- modelGeneVar(sce)
  sce <- clusterplot(sce,dec)
  return (sce)
}
#MNN + Marker Detection (Not used)
mnnmark <- function(sce){
  sce <- pairwiseWilcox(sce, direction="up")
  markers <- getTopMarkers(sce[[1]], sce[[2]], n=10)
  return (markers)
}
