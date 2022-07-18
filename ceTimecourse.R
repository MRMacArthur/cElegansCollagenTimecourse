library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse)
library(org.Ce.eg.db)
library(DEGreport)

myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 legend.key=element_blank())

wormColl <- read.csv("wormMatreosome3.csv")
table(wormColl$Matrisome.Category)

set_collagen <- wormColl[wormColl$Matrisome.Category == "Collagens", ]$WormBase.ID
set_cuticlins <- wormColl[wormColl$Matrisome.Category == "Cuticlins", ]$WormBase.ID
set_cuticularCollagens <- wormColl[wormColl$Matrisome.Category == "Cuticular Collagens", ]$WormBase.ID
set_ecmAffiliated <- wormColl[wormColl$Matrisome.Category == "ECM-affiliated", ]$WormBase.ID
set_ecmGlycoproteins <- wormColl[wormColl$Matrisome.Category == "ECM Glycoproteins", ]$WormBase.ID
set_ecmRegulators <- wormColl[wormColl$Matrisome.Category == "ECM Regulators", ]$WormBase.ID
set_proteoglycans <- wormColl[wormColl$Matrisome.Category == "Proteoglycans", ]$WormBase.ID
set_secretedFactors <- wormColl[wormColl$Matrisome.Category == "Secreted Factors", ]$WormBase.ID
set_adhesionECM <- wormColl[wormColl$Matrisome.Category == "AdhesionReceptorsToECM", ]$WormBase.ID
set_adhesionSig <- wormColl[wormColl$Matrisome.Category == "AdhesionSignalling", ]$WormBase.ID

allWBList <- setNames(list(set_collagen, set_cuticlins, set_cuticularCollagens, set_ecmAffiliated,
           set_ecmGlycoproteins, set_ecmRegulators, set_proteoglycans, set_secretedFactors,
           set_adhesionECM, set_adhesionSig),
           c('collagen', 'cuticlins', 'cuticularCollagens', 'ecmAffiliated',
             'ecmGlycoproteins', 'ecmRegulators', 'proteoglycans', 'secretedFactors',
             'adhesionECM', 'adhesionSignal'))

allWBVec <- c(set_collagen, set_cuticlins, set_cuticularCollagens, set_ecmAffiliated,
                  set_ecmGlycoproteins, set_ecmRegulators, set_proteoglycans, set_secretedFactors,
              set_adhesionECM, set_adhesionSig)

allWBmelt <- reshape2::melt(allWBList)
colnames(allWBmelt) <- c("variable", "type")

tc15h <- read.csv("elegansTC_15h.csv")
tc24h <- read.csv("elegansTC_24h.csv")
tc48h <- read.csv("elegansTC_48h.csv")

tc15h_g <- tc15h[,1]
rownames(tc15h) <- tc15h[,1]
tc15h <- tc15h[-1]
tc15_time <- factor(0:15)

tc24h_g <- tc24h[,1]
rownames(tc24h) <- tc24h[,1]
tc24h <- tc24h[,2:25]
tc24_time <- factor(c(1:24))

tc48h_g <- tc48h[,1]
rownames(tc48h) <- tc48h[,1]
tc48h <- tc48h[,2:45]
tc48_time <- factor(5:48)

tc15x <- DGEList(counts = tc15h, genes = tc15h_g)
tc15x$genes$symbol <- mapIds(org.Ce.eg.db, rownames(tc15x),
                         keytype = "ENSEMBL", column = "SYMBOL")

tc24x <- DGEList(counts = tc24h, genes = tc24h_g)
tc24x$genes$symbol <- mapIds(org.Ce.eg.db, rownames(tc24x),
                             keytype = "ENSEMBL", column = "SYMBOL")

tc48x <- DGEList(counts = tc48h, genes = tc48h_g)
tc48x$genes$symbol <- mapIds(org.Ce.eg.db, rownames(tc48x),
                             keytype = "ENSEMBL", column = "SYMBOL")

tc15_cpm <- data.frame(t(cpm(tc15x)))
tc15_cpm <- data.frame(scale(tc15_cpm))
tc15_cpm$time <- tc15_time
tc15_cpm <- reshape2::melt(tc15_cpm)
tc15_cpm <- tc15_cpm[tc15_cpm$variable %in% allWBVec, ]
tc15_cpm <- merge(tc15_cpm, allWBmelt)
tc15_cpm$time <- as.numeric(tc15_cpm$time)

ggplot(tc15_cpm, aes(x = time, y = value, group = variable)) +
  geom_line(alpha = 0.15) + facet_wrap(~type, ncol = 4)

tc24_cpm <- data.frame(t(cpm(tc24x)))
tc24_cpm <- data.frame(scale(tc24_cpm))
tc24_cpm$time <- tc24_time
tc24_cpm <- reshape2::melt(tc24_cpm)
tc24_cpm <- tc24_cpm[tc24_cpm$variable %in% allWBVec, ]
tc24_cpm <- merge(tc24_cpm, allWBmelt)
tc24_cpm$time <- as.numeric(tc24_cpm$time)

ggplot(tc24_cpm, aes(x = time, y = value, group = variable)) +
  geom_line(alpha = 0.15) + facet_wrap(~type, ncol = 4)

tc48_cpm <- data.frame(t(cpm(tc48x)))
tc48_cpm <- data.frame(scale(tc48_cpm))
tc48_cpm$time <- tc48_time
tc48_cpm <- reshape2::melt(tc48_cpm)
tc48_cpm <- tc48_cpm[tc48_cpm$variable %in% allWBVec, ]
tc48_cpm <- merge(tc48_cpm, allWBmelt)
tc48_cpm$time <- as.numeric(tc48_cpm$time)

ggplot(tc48_cpm, aes(x = time, y = value, group = variable)) +
  geom_line(alpha = 0.15) + facet_wrap(~type, ncol = 4)

tc15Frame <- data.frame(cpm(tc15x, log = T))
tc15Meta <- data.frame(timepoint = (tc15_time))
rownames(tc15Meta) <- colnames(tc15Frame)

tc24Frame <- data.frame(cpm(tc24x, log = T))
tc24Meta <- data.frame(timepoint = (tc24_time))
rownames(tc24Meta) <- colnames(tc24Frame)

tc48Frame <- data.frame(cpm(tc48x, log = T))
tc48Meta <- data.frame(timepoint = (tc48_time))
rownames(tc48Meta) <- colnames(tc48Frame)

# Cuticle collagens
deg15_cutColl <- degPatterns(ma = tc15Frame[rownames(tc15Frame) %in% set_cuticularCollagens,], 
                     metadata = tc15Meta,
                     minc = 5, time = "timepoint")

deg24_cutColl <- degPatterns(ma = tc24Frame[rownames(tc24Frame) %in% set_cuticularCollagens,], 
                     metadata = tc24Meta,
                     minc = 5, time = "timepoint")

deg48_cutColl <- degPatterns(ma = tc48Frame[rownames(tc48Frame) %in% set_cuticularCollagens,], 
                     metadata = tc48Meta,
                     minc = 8, time = "timepoint")

# Cuticle collagens 
# 15h
deg15_cutColl_clust <- deg15_cutColl$df
deg15_cutColl_frame <- data.frame(t(scale(t(tc15Frame[rownames(tc15Frame) %in% 
                                   deg15_cutColl_clust$genes,]))))
#deg15_cutColl_frame <- tc15Frame[rownames(tc15Frame) %in% 
#                                   deg15_cutColl_clust$genes,]

deg15_cutColl_frame <- merge(deg15_cutColl_frame, 
                             deg15_cutColl_clust,
                             by = 0)

deg15_cutColl_melt <- reshape2::melt(deg15_cutColl_frame[,2:19], 
                                     id.vars = c("genes", "cluster"))
deg15_cutColl_melt$timeNum <- rep(0:15, each = 58)
deg15_cutColl_melt$cluster <- factor(deg15_cutColl_melt$cluster)

ggplot(deg15_cutColl_melt, aes(x = timeNum, y = value, 
                               group = cluster, fill = cluster,
                               color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg15_genes <- merge(deg15_cutColl_clust, tc15x$genes, by = 0)

# 24h
deg24_cutColl_clust <- deg24_cutColl$df
deg24_cutColl_frame <- data.frame(t(scale(t(tc24Frame[rownames(tc24Frame) %in% 
                                   deg24_cutColl_clust$genes,]))))
#deg24_cutColl_frame <- tc24Frame[rownames(tc24Frame) %in% 
#                                   deg24_cutColl_clust$genes,]))))

deg24_cutColl_frame <- merge(deg24_cutColl_frame, 
                             deg24_cutColl_clust,
                             by = 0)

deg24_cutColl_melt <- reshape2::melt(deg24_cutColl_frame[,2:27], 
                                     id.vars = c("genes", "cluster"))
deg24_cutColl_melt$timeNum <- rep(1:24, each = 63)
deg24_cutColl_melt$cluster <- factor(deg24_cutColl_melt$cluster)

ggplot(deg24_cutColl_melt, aes(x = timeNum, y = value, 
                               group = cluster, fill = cluster,
                               color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg24_genes <- merge(deg24_cutColl_clust, tc24x$genes, by = 0)

# 48h
deg48_cutColl_clust <- deg48_cutColl$df
deg48_cutColl_clust[deg48_cutColl_clust$cluster == 9, ]$cluster <- 7
deg48_cutColl_frame <- data.frame(t(scale(t(tc48Frame[rownames(tc48Frame) %in% 
                                   deg48_cutColl_clust$genes,]))))
#deg48_cutColl_frame <- tc48Frame[rownames(tc48Frame) %in% 
#                                   deg48_cutColl_clust$genes,]))))

deg48_cutColl_frame <- merge(deg48_cutColl_frame, 
                             deg48_cutColl_clust,
                             by = 0)

deg48_cutColl_melt <- reshape2::melt(deg48_cutColl_frame[,2:47], 
                                     id.vars = c("genes", "cluster"))
deg48_cutColl_melt$timeNum <- rep(5:48, each = 65)
deg48_cutColl_melt$cluster <- factor(deg48_cutColl_melt$cluster)

ggplot(deg48_cutColl_melt, aes(x = timeNum, y = value, 
                               group = cluster, fill = cluster,
                               color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg48_genes <- merge(deg48_cutColl_clust, tc48x$genes, by = 0)

# Adhesion Receptors to ECM
deg15_Adhes <- degPatterns(ma = tc15Frame[rownames(tc15Frame) %in% c(set_adhesionSig, set_adhesionECM),], 
                             metadata = tc15Meta,
                             minc = 3, time = "timepoint")

deg24_Adhes <- degPatterns(ma = tc24Frame[rownames(tc24Frame) %in% c(set_adhesionSig, set_adhesionECM),], 
                             metadata = tc24Meta,
                             minc = 3, time = "timepoint")

deg48_Adhes <- degPatterns(ma = tc48Frame[rownames(tc48Frame) %in% c(set_adhesionSig, set_adhesionECM),], 
                             metadata = tc48Meta,
                             minc = 2, time = "timepoint")


# Adhesion 
# 15h
deg15_Adhes_clust <- deg15_Adhes$df
#deg15_Adhes_frame <- data.frame(t(scale(t(tc15Frame[rownames(tc15Frame) %in% 
#                                   deg15_Adhes_clust$genes,]))))
deg15_Adhes_frame <- tc15Frame[rownames(tc15Frame) %in%
                                deg15_Adhes_clust$genes,]

deg15_Adhes_frame <- merge(deg15_Adhes_frame, 
                             deg15_Adhes_clust,
                             by = 0)

deg15_Adhes_melt <- reshape2::melt(deg15_Adhes_frame[,2:19], 
                                     id.vars = c("genes", "cluster"))
deg15_Adhes_melt$timeNum <- rep(0:15, each = 18)
deg15_Adhes_melt$cluster <- factor(deg15_Adhes_melt$cluster)

ggplot(deg15_Adhes_melt, aes(x = timeNum, y = value, 
                               group = cluster, fill = cluster,
                               color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg15_genes <- merge(deg15_Adhes_clust, tc15x$genes, by = 0)

# 24h
deg24_Adhes_clust <- deg24_Adhes$df
#deg24_Adhes_frame <- data.frame(t(scale(t(tc24Frame[rownames(tc24Frame) %in% 
#                                                      deg24_Adhes_clust$genes,]))))
deg24_Adhes_frame <- tc24Frame[rownames(tc24Frame) %in% 
                                 deg24_Adhes_clust$genes,]

deg24_Adhes_frame <- merge(deg24_Adhes_frame, 
                           deg24_Adhes_clust,
                           by = 0)

deg24_Adhes_melt <- reshape2::melt(deg24_Adhes_frame[,2:27], 
                                   id.vars = c("genes", "cluster"))
deg24_Adhes_melt$timeNum <- rep(1:24, each = 20)
deg24_Adhes_melt$cluster <- factor(deg24_Adhes_melt$cluster)

ggplot(deg24_Adhes_melt, aes(x = timeNum, y = value, 
                             group = cluster, fill = cluster,
                             color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg24_genes <- merge(deg24_Adhes_clust, tc24x$genes, by = 0)

# 48h
deg48_Adhes_clust <- deg48_Adhes$df
#deg48_Adhes_frame <- data.frame(t(scale(t(tc48Frame[rownames(tc48Frame) %in% 
#                                                      deg48_Adhes_clust$genes,]))))
deg48_Adhes_frame <- tc48Frame[rownames(tc48Frame) %in%
                                 deg48_Adhes_clust$genes,]

deg48_Adhes_frame <- merge(deg48_Adhes_frame, 
                           deg48_Adhes_clust,
                           by = 0)

deg48_Adhes_melt <- reshape2::melt(deg48_Adhes_frame[,2:47], 
                                   id.vars = c("genes", "cluster"))
deg48_Adhes_melt$timeNum <- rep(5:48, each = 26)
deg48_Adhes_melt$cluster <- factor(deg48_Adhes_melt$cluster)

ggplot(deg48_Adhes_melt, aes(x = timeNum, y = value, 
                             group = cluster, fill = cluster,
                             color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme

deg48_genes <- merge(deg48_Adhes_clust, tc48x$genes, by = 0)


##### ECM Regulators
set_ecmRegulators
# Adhesion Receptors to ECM
deg15_ecmReg <- degPatterns(ma = tc15Frame[rownames(tc15Frame) %in% set_ecmRegulators,], 
                           metadata = tc15Meta,
                           minc = 8, time = "timepoint")

deg24_ecmReg <- degPatterns(ma = tc24Frame[rownames(tc24Frame) %in% set_ecmRegulators,], 
                           metadata = tc24Meta,
                           minc = 8, time = "timepoint")

deg48_ecmReg <- degPatterns(ma = tc48Frame[rownames(tc48Frame) %in% set_ecmRegulators,], 
                           metadata = tc48Meta,
                           minc = 7, time = "timepoint")


# ECM Reg
# 15h
deg15_ecmReg_clust <- deg15_ecmReg$df
deg15_ecmReg_frame <- data.frame(t(scale(t(tc15Frame[rownames(tc15Frame) %in% 
                                                      deg15_ecmReg_clust$genes,]))))
#deg15_ecmReg_frame <- tc15Frame[rownames(tc15Frame) %in%
#                                deg15_ecmReg_clust$genes,]

deg15_ecmReg_frame <- merge(deg15_ecmReg_frame, 
                           deg15_ecmReg_clust,
                           by = 0)

deg15_ecmReg_melt <- reshape2::melt(deg15_ecmReg_frame[,2:19], 
                                   id.vars = c("genes", "cluster"))
deg15_ecmReg_melt$timeNum <- rep(0:15, each = 86)
deg15_ecmReg_melt$cluster <- factor(deg15_ecmReg_melt$cluster)

ggplot(deg15_ecmReg_melt, aes(x = timeNum, y = value, 
                             group = cluster, fill = cluster,
                             color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme #+ facet_wrap(~cluster, ncol = 3)

deg15_genes <- merge(deg15_ecmReg_clust, tc15x$genes, by = 0)

# ECM Reg
# 24h
deg24_ecmReg_clust <- deg24_ecmReg$df
deg24_ecmReg_frame <- data.frame(t(scale(t(tc24Frame[rownames(tc24Frame) %in% 
                                                      deg24_ecmReg_clust$genes,]))))
#deg24_ecmReg_frame <- tc24Frame[rownames(tc24Frame) %in%
#                                  deg24_ecmReg_clust$genes,]

deg24_ecmReg_frame <- merge(deg24_ecmReg_frame, 
                            deg24_ecmReg_clust,
                            by = 0)

deg24_ecmReg_melt <- reshape2::melt(deg24_ecmReg_frame[,2:27], 
                                    id.vars = c("genes", "cluster"))
deg24_ecmReg_melt$timeNum <- rep(1:24, each = 90)
deg24_ecmReg_melt$cluster <- factor(deg24_ecmReg_melt$cluster)

ggplot(deg24_ecmReg_melt, aes(x = timeNum, y = value, 
                              group = cluster, fill = cluster,
                              color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme #+ facet_wrap(~cluster, ncol = 3)

deg24_genes <- merge(deg24_ecmReg_clust, tc24x$genes, by = 0)


# ECM Reg
# 48h
deg48_ecmReg_clust <- deg48_ecmReg$df
deg48_ecmReg_frame <- data.frame(t(scale(t(tc48Frame[rownames(tc48Frame) %in% 
                                                      deg48_ecmReg_clust$genes,]))))
#deg48_ecmReg_frame <- tc48Frame[rownames(tc48Frame) %in%
#                                  deg48_ecmReg_clust$genes,]

deg48_ecmReg_frame <- merge(deg48_ecmReg_frame, 
                            deg48_ecmReg_clust,
                            by = 0)

deg48_ecmReg_melt <- reshape2::melt(deg48_ecmReg_frame[,2:47], 
                                    id.vars = c("genes", "cluster"))
deg48_ecmReg_melt$timeNum <- rep(5:48, each = 83)
deg48_ecmReg_melt$cluster <- factor(deg48_ecmReg_melt$cluster)

ggplot(deg48_ecmReg_melt, aes(x = timeNum, y = value, 
                              group = cluster, fill = cluster,
                              color = cluster)) +
  geom_smooth(se = T, span = 0.1) +
  geom_point(alpha = 0.2) + labs(x = "Time (h)", y = "Expression z-score") +
  myTheme #+ facet_wrap(~cluster, ncol = 3)

deg48_genes <- merge(deg48_ecmReg_clust, tc48x$genes, by = 0)

table(deg15_ecmReg_clust$cluster)
clust_ecmReg_15h_1 <- deg15_ecmReg_clust[deg15_ecmReg_clust$cluster == 1, ]$genes
clust_ecmReg_15h_2 <- deg15_ecmReg_clust[deg15_ecmReg_clust$cluster == 2, ]$genes
clust_ecmReg_15h_3 <- deg15_ecmReg_clust[deg15_ecmReg_clust$cluster == 3, ]$genes
clust_ecmReg_15h_4 <- deg15_ecmReg_clust[deg15_ecmReg_clust$cluster == 5, ]$genes
clust_ecmReg_15h_5 <- deg15_ecmReg_clust[deg15_ecmReg_clust$cluster == 11, ]$genes




