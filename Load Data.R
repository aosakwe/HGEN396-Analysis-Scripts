##Loading Bioconductor Packages
library(BiocManager)
library(DropletUtils)
library(scater)
library(scuttle)
library(scran)
library(batchelor)
library(EnsDb.Hsapiens.v86)
library(limma)
library(celldex)
library(SingleR)
library(scRNAseq)
library(AnnotationHub)
library(pheatmap)
library(scDblFinder)
# loading individual datasets as "sceXXX" where XXX is the donor number
sce24 <- "~/Desktop/datasets/ND/cellranger24"
sce26 <- "~/Desktop/datasets/ND/cellranger26"
sce27 <- "~/Desktop/datasets/ND/cellranger27"
sce29 <- "~/Desktop/datasets/ND/cellranger29"
sce35 <- "~/Desktop/datasets/ND/cellranger35"
sce36 <- "~/Desktop/datasets/ND/cellranger36"
sce37 <- "~/Desktop/datasets/ND/cellranger37"
sce40 <- "~/Desktop/datasets/ND/cellranger40"
sce45 <- "~/Desktop/datasets/ND/cellranger45"
sce49 <- "~/Desktop/datasets/ND/cellranger49"
sce50 <- "~/Desktop/datasets/ND/cellranger50"
sce52 <- "~/Desktop/datasets/ND/cellranger52"
sce56 <- "~/Desktop/datasets/ND/cellranger56"
sce59 <- "~/Desktop/datasets/ND/cellranger59"

## T2D   patients
sce51 <- "~/Desktop/datasets/T2D/cellranger51"
sce57 <- "~/Desktop/datasets/T2D/cellranger57"
sce58 <- "~/Desktop/datasets/T2D/cellranger58"
