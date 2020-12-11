## to generate the annotated muraro
edb <- AnnotationHub()[["AH73881"]]
edb
gene.symb <- sub("__chr.*$", "", rownames(sce.muraro))
gene.symb
gene.ids <- mapIds(edb, keys=gene.symb, keytype = "SYMBOL", column = "GENEID")
gene.ids
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
keep <- !is.na(gene.symb) & !duplicated(gene.symb)
sce.muraro <- sce.muraro[keep,]
rownames(sce.muraro) <- gene.ids[keep]
rownames(sce.muraro)
test <- sce.muraro[keep,]
rownames(sce.muraro) <- gene.symb[keep]
rownames(test)
rownames(sce.muraro)
sce.muraro 


merged
gene.ids <- mapIds(edb, keys=rownames(merged), keytype = "SYMBOL", column = "GENEID")
gene.ids
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
merged <- merged[keep,]
rownames(merged) <- gene.ids[keep]
rownames(merged)

chosen.hvgs
annocell
gid <- mapIds(edb, keys=chosen.hvgs, keytype = "SYMBOL", column = "GENEID")
keep <-  !is.na(gid) & !duplicated(gid)
annocell <- gid[keep]
annocell

anno
