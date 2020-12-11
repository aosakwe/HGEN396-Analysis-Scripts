## Labelling cells by batch
sce59$batch <- "ND59"
sce24$batch <- "ND24"
sce26$batch <- "ND26"
sce27$batch <- "ND27"
sce29$batch <- "ND29"
sce35$batch <- "ND35"
sce36$batch <- "ND36"
sce37$batch <- "ND37"
sce40$batch <- "ND40"
sce45$batch <- "ND45"
sce49$batch <- "ND49"
sce50$batch <- "ND50"
sce52$batch <- "ND52"
sce56$batch <- "ND56"
sce51$batch <- "T2D51"
sce57$batch <- "T2D57"
sce58$batch <- "T2D58"
uncorrected <- cbind(sce24,sce26,sce27,sce29,sce35,sce36,sce37,sce40,sce45,sce49,sce50,sce52,sce56,sce59,sce51,sce57,sce58)
###################
#######################
uncorrected <- logNormCounts(uncorrected)
uncorrected
dec <- modelGeneVar(uncorrected, block=uncorrected$batch)
dec
chosen.hvgs <- rownames(dec)[dec$bio > 0]
chosen.hvgs

############################### Uncorrected UMAP
set.seed(0010101010)
uncorrected <- runPCA(uncorrected, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())
snn.gr <- buildSNNGraph(uncorrected, use.dimred="PCA")
groups <- igraph::cluster_walktrap(snn.gr)$membership
colLabels(uncorrected) <- factor(groups)
set.seed(1111001)
uncorrected <- runUMAP(uncorrected, dimred="PCA")
plotUMAP(uncorrected, colour_by="batch")
################################## Corrected UMAP
set.seed(01001001)
merged <- correctExperiments(uncorrected, batch=uncorrected$batch, subset.row = chosen.hvgs, PARAM = FastMnnParam())
merged
g <- buildSNNGraph(merged, use.dimred ="corrected")
clusters <- igraph::cluster_louvain(g)
clusters <- factor(clusters$membership)
colLabels(merged) <- clusters
merged <- runUMAP(merged, dimred="corrected")
merged <- runTSNE(merged, dimred="corrected")
table(colLabels(merged), merged$T2D)
merged$label

table(Cluster=clusters, Batch=merged$batch)
plotTSNE(merged, colour_by="batch", text_by="label")#Ad
plotUMAP(merged, colour_by="batch", text_by="label")
plotUMAP(merged, colour_by="label", text_by="label")
plotUMAP(merged, colour_by="T2D", text_by="label")
##############################################Cell Annotation
### Note: See Annotation Conversion script for converting Gene IDs to match dataset

sce.muraro <- MuraroPancreasData()
sce.muraro <- logNormCounts(sce.muraro)

sce.muraro <- sce.muraro[,!is.na(sce.muraro$label) & 
                           sce.muraro$label!="unclear"]
table(sce.muraro$label)
library(AnnotationHub)
?AnnotationHub
test <- AnnotationHub()
test[test$species == "Homo sapiens"]
hs.db <- AnnotationHub()[["AH73881"]]
hs.db
hs.exons <- exonsBy(hs.db, by="gene")
hs.exons <- reduce(hs.exons)
hs.len <- sum(width(hs.exons))

library(scuttle)
available <- intersect(rownames(merged), names(hs.len))

fpkm <- calculateFPKM(merged[available,], hs.len[available])
pred <- SingleR(test=fpkm, ref=sce.muraro, 
                      labels=sce.muraro$label, de.method="wilcox")
table(pred$labels)
merged
sce.muraro 
merged$label
tab <- table(pred$labels, merged$label)
table(pred$labels, merged$batch)
tab
library(pheatmap)
pheatmap(log2(tab+10), color = colorRampPalette(c("white", "blue"))(101))
##########################

##Heatmap for all clusters

markers <- findMarkers(merged)
markers
chosen <- 15
interesting <- markers[[chosen]]
colnames(interesting)
interesting[1:10,1:4]
best.set <- interesting[interesting$Top <= 6,]
logFCs <- getMarkerEffects(best.set)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))
#Clusters of Interest


#######
emergency <- merged ## stored object in case manipulation generated errors
merged <- emergency

#merged <- mapIds(edb, keys=rownames(merged), keytype = "GENEID", column = "SEQNAME")
t2d <- merged[,merged$T2D == TRUE & merged$label == 15]
t2d$label <- "T2D Beta"
t2d
nd <- merged[,merged$T2D == FALSE & merged$label == 15]
nd$label <- "ND Beta"
nd
t2da <-merged[,merged$T2D == TRUE & merged$label == 11]
t2da$label <- "T2D Alpha"
t2da
nda <- merged[,merged$T2D == FALSE & merged$label == 11]
nda$label <- "ND Alpha"
nda
acnd <- merged[,merged$T2D == FALSE & merged$label == 10]
acnd$label <- "ND Acinar"
acnd
act2d <- merged[,merged$T2D == TRUE & merged$label == 10]
act2d$label <- "T2D Acinar"
act2d
cc <- cbind(t2d,nd,t2da,nda,acnd,act2d)
cc <- cbind(t2d,nd,t2da,nda)
cc$label
table( Cluster=cc$label, Batch=cc$batch)
ccmark <- findMarkers(cc)
chosen <- "ND Beta"
inter <-ccmark[[chosen]]
colnames(inter)
bs <- inter[inter$Top <= 6,]
logs <- getMarkerEffects(bs)
logs
pheatmap(logs, breaks = seq(-5, 5,length.out = 101))
#########################
beta <- cbind(t2d,nd)
findMarkers(beta)
bmark <- findMarkers(beta)
btest <- bmark[[chosen]]
colnames(btest)
BS <- btest[btest$Top <= 30,]
BS
beta
rownames(BS)
rownames(beta)
write.table(rownames(BS), file = "Beta_markers.txt", sep = " ")
###############
alpha <- cbind(t2da,nda)
amark <- findMarkers(alpha)
atest <- amark[["ND Alpha"]]
AS <- atest[atest$Top <= 30,]
write.table(rownames(AS), file = "Alpha_markers.txt", sep = " ")
write.table(rownames(alpha), file = "Alpha.txt", sep = " ")
?findMarkers
### both 
table(merged$label, merged$batch)
R
##############################
acinar <- cbind(acnd, act2d)
table(acinar$label)
acmark <- findMarkers(acinar)
actest <- acmark[["ND Acinar"]]
actest
ACS <-actest[actest$Top <= 30,]
rownames(ACS)
write.table(rownames(ACS), file = "Acinar_markers.txt", sep = " ")
write.table(rownames(acinar), file = "Acinar.txt", sep = " ")
table(pred$labels, merged$T2D)

#######Doublets?
dbl <- findDoubletClusters(merged)
chosen.doublet <- rownames(dbl)[isOutlier(dbl$num.de, type="lower", log=TRUE)]
chosen.doublet
dbl.dens <- computeDoubletDensity(merged, subset.row=chosen.hvgs, d=ncol(reducedDim(merged)))
merged$DoubletScore <- dbl.dens
plotUMAP(merged, colour_by = "DoubletScore", text_by="label")
####test to plot acinar ######
rownames(merged)
goodclusts <- merged[,merged$label != 1 & merged$label != 2 & merged$label != 5]
merged$Cluster <- merged$label
plotExpression(merged, features ="PRSS1", x="Cluster", colour_by = "Cluster")
"PRSS1" %in% rownames(merged)
