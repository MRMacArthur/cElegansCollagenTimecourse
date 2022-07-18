library(Seurat)
library(dittoSeq)

load("ageS.RData")
wormColl <- read.csv("wormMatreosome3.csv")
annData <- read.csv("annotations.csv")

matGenes <- wormColl$symbol

set_collagen <- wormColl[wormColl$Matrisome.Category == "Collagens" &
                           wormColl$symbol %in% rownames(ageS), ]$symbol
set_cuticlins <- wormColl[wormColl$Matrisome.Category == "Cuticlins" &
                            wormColl$symbol %in% rownames(ageS), ]$symbol
set_cuticularCollagens <- wormColl[wormColl$Matrisome.Category == "Cuticular Collagens" &
                                     wormColl$symbol %in% rownames(ageS), ]$symbol
set_ecmAffiliated <- wormColl[wormColl$Matrisome.Category == "ECM-affiliated" &
                                wormColl$symbol %in% rownames(ageS), ]$symbol
set_ecmGlycoproteins <- wormColl[wormColl$Matrisome.Category == "ECM Glycoproteins" &
                                   wormColl$symbol %in% rownames(ageS), ]$symbol
set_ecmRegulators <- wormColl[wormColl$Matrisome.Category == "ECM Regulators" &
                                wormColl$symbol %in% rownames(ageS), ]$symbol
set_proteoglycans <- wormColl[wormColl$Matrisome.Category == "Proteoglycans" &
                                wormColl$symbol %in% rownames(ageS), ]$symbol
set_secretedFactors <- wormColl[wormColl$Matrisome.Category == "Secreted Factors" &
                                  wormColl$symbol %in% rownames(ageS), ]$symbol
set_adhesionECM <- wormColl[wormColl$Matrisome.Category == "AdhesionReceptorsToECM" &
                              wormColl$symbol %in% rownames(ageS), ]$symbol
set_adhesionSig <- wormColl[wormColl$Matrisome.Category == "AdhesionSignalling" &
                              wormColl$symbol %in% rownames(ageS), ]$symbol


ageS$cellType <- plyr::mapvalues(
  x = ageS$annotate_name,
  from = annData$annotation,
  to = annData$cellType
)

ageS2 <- ScaleData(ageS)

ageMat <- subset(ageS, features = matGenes)
matGenes <- unique(matGenes[matGenes %in% rownames(ageS)])

dittoHeatmap(ageS,
             genes = set_cuticularCollagens,
             annot.by = c("timepoint", "cellType"))

DoHeatmap(ageS2,
          features = set_cuticularCollagens,
          group.by = "cellType")

DimPlot(ageS, group = "timepoint")

levels(ageS@meta.data$annotate_name)

VlnPlot(ageS,
        features = c("col-106"),
        group = "cellType",
        split.by = "timepoint")

head(wormColl, 25)


VlnPlot(m1, features = c("annotate_name-10"),
        group.by = "annotate_name", split.by = "time")

FeaturePlot(ageS, features = c(""))

