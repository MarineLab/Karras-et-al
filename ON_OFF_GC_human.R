library(nichenetr)
GC_malignant_subset_BT <- readRDS("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Malignant_BT.rds")
#################################### Karras et all ON OFF
################ Signatures from Landscape paper ###########################################
##########################        Neural_like        ##########################
genes<-c("Tfap2b","Prnp","Mef2c","Gfra1","Cd200","Syt11","Thsd7a","Cxxc4","Sema5a","Tbx3","Kif26b","Efhd1","Neto2","Hmcn1","Igsf10","Olfml3","Rgs2","Nt5e","Morc4","Aqp1","Gfra2","Wnt4","Abcg2","Elovl5","Emilin1","Fibin")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Neural_like")
geneSets
cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Neural_like_mouse_NEW.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Neural_like<-getAUC(cells_AUC)
Neural_like<-t(Neural_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Neural_like)
FeaturePlot(GC_malignant_subset_BT, features = "Neural_like", label = T) +NoAxes()
VlnPlot(GC_malignant_subset_BT, features = "Neural_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Neural_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Neural_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Neural_like", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()

##########################        Melanocytic_OXPHOS        ##########################
genes<-c("Mlph","Mitf","Mlana","Pmel","Slc45a2","Apoe","Dct","Gpnmb","Tyr","Cited1","Bace2","Cox17","Cox7a2","Ndufb2","Ndufb4","Uqcr10","Uqcrb")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Melanocytic_OXPHOS")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Melanocytic_OXPHOS_mouse_NEW.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Melanocytic_OXPHOS<-getAUC(cells_AUC)
Melanocytic_OXPHOS<-t(Melanocytic_OXPHOS)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Melanocytic_OXPHOS)
FeaturePlot(GC_malignant_subset_BT, features = "Melanocytic_OXPHOS", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Melanocytic_OXPHOS", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Melanocytic_OXPHOS", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Melanocytic_OXPHOS", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Melanocytic_OXPHOS", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Trans_reg        ##########################
genes<-c("Pabpn1","Ccnl1","Pnisr","Sfpq","Ccnl2","Bclaf1","Hnrnpu","Rbm25","Luc7l2","Ddx17","Cpsf6","Snrnp70","Srek1","Prpf4b","Srrm2","Rbm39","Son","Ddx5","Prpf38b","Pan3","Tra2a","Hnrnph1","Fus")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Trans_reg")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Trans_reg_mouse_NEW.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Trans_reg<-getAUC(cells_AUC)
Trans_reg<-t(Trans_reg)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Trans_reg)
FeaturePlot(GC_malignant_subset_BT, features = "Trans_reg", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Trans_reg", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Trans_reg", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Trans_reg", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Trans_reg", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Glycolytic_UPR        ##########################
genes<-c("Bnip3","Tpi1","Slc2a1","Mif","Vldlr","Hk2","Vegfa","Ldha","Pfkl","P4ha1","Fam162a","Bhlhe40","Pgk1","Aldoa","Pfkp","Pdk1","Hspa9","Pgam1","Pkm","Atf4")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Glycolytic_UPR")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Glycolytic_UPR_mouse_NEW.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Glycolytic_UPR<-getAUC(cells_AUC)
Glycolytic_UPR<-t(Glycolytic_UPR)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Glycolytic_UPR)
FeaturePlot(GC_malignant_subset_BT, features = "Glycolytic_UPR", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Glycolytic_UPR", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Glycolytic_UPR", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Glycolytic_UPR", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Glycolytic_UPR", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Stem_like         ##########################
genes<-c("Gtse1","Cep170b","Tap1","Ercc5","Kank3","Rap2a","Ei24","Zmat3","Fas","Igdcc4","Notch3","Vcan","Mybl1","Dusp15","Dgka","Nrp1","Slc27a3","Slc19a2","Dgki","Siva1","Rbm38","Cad","Mcm7","Nme4","Thyn1","Gpx3","Arl4c","Ass1","Efnb2","Frrs1","Sdc1","Itga6","Icam1","Chst3","Nes","Sox4","Fat1","Angptl2")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Stem_like")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Stem_like_mouse_NEW.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stem_like<-getAUC(cells_AUC)
Stem_like<-t(Stem_like)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Stem_like)
FeaturePlot(GC_malignant_subset_BT, features = "Stem_like", label = T) +NoAxes()
VlnPlot(GC_malignant_subset_BT, features = "Stem_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Stem_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Stem_like", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()



##########################        Immune        ##########################
GC_malignant_subset_BT$Immune <- NULL
genes<-c("Gbp4","Gbp6","Irf1","Nlrc5","Stat2","Cd74","Il12rb1","B2m","Gbp2","Ifit3","Gbp3","Psmb8","Stat1","Ifit1","Isg15","Psmb9","Gbp7","Xaf1","Tap1","Gbp5","Usp18","Ciita","Irf7","Ube2l6","Tap2","Psme1","Ifitm3","Psme2","Tapbp","Ifi35","Eif2ak2","Irf9","Ddx58","Trim25","Ifit2","Irf8","H2-DMa","H2-Aa","H2-DMb1","H2-Eb1","H2-Ab1","H2-T23")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
genes
geneSets <- GeneSet(genes, setName="Immune")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Immune_mouse_NEW.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Immune)
FeaturePlot(GC_malignant_subset_BT, features = "Immune", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Immune", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Immune", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Immune", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


##########################        Mesenchymal        ##########################
genes<-c("Fap","Lama2","Prrx1","Slit3","Abi3bp","Loxl1","Cdh11","Col6a3","Col6a2","Loxl2","Mfap5","Fbn1","Col4a2","Pcolce","Lum","Col4a1","Col5a2","Thy1","Fbn2","Bgn","Tgfbi","Pdgfrb","Sulf1","Inhba","Col1a1","Col3a1","Col1a2","Ctsk","Dcn","Serpinf1","Sparc","Fstl1")
genes <- convert_mouse_to_human_symbols(genes)
genes <- unique(genes)
geneSets <- GeneSet(genes, setName="Mesenchymal")
geneSets
#cells_rankings <- AUCell_buildRankings(GC_malignant_subset_BT@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Mesenchymal_mouse_NEW.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mesenchymal<-getAUC(cells_AUC)
Mesenchymal<-t(Mesenchymal)
GC_malignant_subset_BT@meta.data<-cbind(GC_malignant_subset_BT@meta.data, Mesenchymal)
FeaturePlot(GC_malignant_subset_BT, features = "Mesenchymal", label = T)
VlnPlot(GC_malignant_subset_BT, features = "Mesenchymal", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Mesenchymal", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
VlnPlot(GC_malignant_subset_BT, features = "Mesenchymal", pt.size = 0, group.by = "GC number") + theme(legend.position = 'none')
dev.off()


############## heatmap 
####################################heatmap#########################################################
reglist <- c("Malignant_clusters","Immune","Melanocytic_OXPHOS","Mesenchymal",  "Neural_like","Stem_like", "Trans_reg", "Glycolytic_UPR")
matr <- FetchData(GC_malignant_subset_BT, vars = reglist, cells = NULL) # %>%
#apply(2, function(x) (x - min(x)) / (max(x) - min(x))) # %>%
#t
library(scales)
matr <- matr[ , !duplicated(colnames(matr))]
matr <- matr %>% group_by(Malignant_clusters) %>% summarise_all(mean)
#
Malignant_clusters <- levels(GC_malignant_subset_BT$Malignant_clusters)
Malignant_clusters
color_1 <- setNames(hue_pal()(11), Malignant_clusters)
color_list <- list(color_1 = Malignant_clusters) 
color_list
#names(color_list) <- list(Malignant_clusters_old)
Malignant_clusters <- matr$Malignant_clusters
ha = ComplexHeatmap::HeatmapAnnotation(Malignant_clusters = Malignant_clusters,
                                       annotation_name_side = "left",
                                       annotation_name_rot = 180,
                                       show_annotation_name = FALSE
                                       # col = color_list
)

#head(matr)
#matr <-apply(matr, 1, as.numeric)
head(matr)

matr$Malignant_clusters <- NULL
mat_num <- matrix(as.numeric(as.character((matr),    # Convert to numeric matrix
                                          ncol = ncol(matr))))
matr <- t(matr)
mat_scaled = t(apply(matr, 1, scale))
ht<-ComplexHeatmap::Heatmap(mat_scaled, 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            top_annotation = ha,
                            #col = color_list,
                            column_split = Malignant_clusters)
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/scores_heatmap_av_Karras_et_al_NEW.pdf", width = 10, height = 6)
draw(ht)
dev.off()



################## ON OFF
GC_malignant_subset_BT$NEURAL_ON_OFF <- ifelse(GC_malignant_subset_BT$Neural_like > 0.1,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "NEURAL_ON_OFF")
NEURAL_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(NEURAL_ON_OFF)
Neural_like_ON <- NEURAL_ON_OFF %>% subset(NEURAL_ON_OFF =="ON")
Neural_like_ON$ON_ALL <- paste("Neural_like", Neural_like_ON$NEURAL_ON_OFF, sep="_")

GC_malignant_subset_BT$MELANOCYTIC_ON_OFF <- ifelse(GC_malignant_subset_BT$Melanocytic_OXPHOS > 0.5,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "MELANOCYTIC_ON_OFF")
MELANOCYTIC_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(MELANOCYTIC_ON_OFF)
Melanocytic_ON <- MELANOCYTIC_ON_OFF %>% subset(MELANOCYTIC_ON_OFF =="ON")
Melanocytic_ON$ON_ALL <- paste("Melanocytic", Melanocytic_ON$MELANOCYTIC_ON_OFF, sep="_")

GC_malignant_subset_BT$TRANS_REG_ON_OFF <- ifelse(GC_malignant_subset_BT$Trans_reg > 0.5,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "TRANS_REG_ON_OFF")
TRANS_REG_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(TRANS_REG_ON_OFF)
Trans_reg_ON <- TRANS_REG_ON_OFF %>% subset(TRANS_REG_ON_OFF =="ON")
Trans_reg_ON$ON_ALL <- paste("Trans_reg", Trans_reg_ON$TRANS_REG_ON_OFF, sep="_")


GC_malignant_subset_BT$UPR_ON_OFF <- ifelse(GC_malignant_subset_BT$Glycolytic_UPR > 0.5,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "UPR_ON_OFF")
UPR_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(UPR_ON_OFF)
UPR_ON <- UPR_ON_OFF %>% subset(UPR_ON_OFF =="ON")
UPR_ON$ON_ALL <- paste("UPR", UPR_ON$UPR_ON_OFF, sep="_")


GC_malignant_subset_BT$STEM_ON_OFF <- ifelse(GC_malignant_subset_BT$Stem_like > 0.12,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "STEM_ON_OFF")
STEM_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(STEM_ON_OFF)
preEMT_ON <- STEM_ON_OFF %>% subset(STEM_ON_OFF =="ON")
preEMT_ON$ON_ALL <- paste("preEMT", preEMT_ON$STEM_ON_OFF, sep="_")


GC_malignant_subset_BT$MESENCHYMAL_ON_OFF <- ifelse(GC_malignant_subset_BT$Mesenchymal > 0.2,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "MESENCHYMAL_ON_OFF")
MESENCHYMAL_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(MESENCHYMAL_ON_OFF)
Mesenchymal_ON <- MESENCHYMAL_ON_OFF %>% subset(MESENCHYMAL_ON_OFF =="ON")
Mesenchymal_ON$ON_ALL <- paste("Mesenchymal", Mesenchymal_ON$MESENCHYMAL_ON_OFF, sep="_")

DimPlot(GC_malignant_subset_BT, group.by = c("STEM_ON_OFF"), cols = c("grey", "blue")) +NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/STEM_ON_OFF_markers.pdf",useDingbats=FALSE,family="sans")
DimPlot(GC_malignant_subset_BT, group.by = c("NEURAL_ON_OFF"), cols = c("grey", "blue"))+NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/NEURAL_ON_OFF_markers.pdf",useDingbats=FALSE,family="sans")
DimPlot(GC_malignant_subset_BT, group.by = c("MELANOCYTIC_ON_OFF"), cols = c("grey", "blue"))+NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/MELANOCYTIC_ON_OFF_markers.pdf",useDingbats=FALSE,family="sans")
DimPlot(GC_malignant_subset_BT, group.by = c("MESENCHYMAL_ON_OFF"), cols = c("grey", "blue"))+NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/MESENCHYMAL_ON_OFF_ON_OFF_markers.pdf",useDingbats=FALSE,family="sans")

pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/MALIGNANT/Stem_like_ON_OFF_markers.pdf", width = 10, height = 10)
VlnPlot(GC_malignant_subset_BT, features = c("GTSE1", "NOTCH3", "VCAN", "MYBL1", "CEP170B", "DUSP15", "DGKA", "SDC1", "NRP1", "TAP1", "ERCC5", "SLC27A3"), group.by = "STEM_ON_OFF")
VlnPlot(GC_malignant_subset_BT, features = c("KANK3", "SLC19A2", "ITGA6", "DGKI", "SIVA1", "RAP2A", "RBM38", "ICAM1", "CAD", "CHST3", "MCM7", "AKAP1"), group.by = "STEM_ON_OFF")
VlnPlot(GC_malignant_subset_BT, features = c("NME4", "THYN1", "GPX3", "ARL4C", "ASS1", "EFNB2", "FRRS1", "EI24", "ZMAT3", "NES"), group.by = "STEM_ON_OFF")
dev.off()



GC_malignant_subset_BT$Immune_ON_OFF <- ifelse(GC_malignant_subset_BT$Immune > 0.3,  "ON", "OFF")
DimPlot(GC_malignant_subset_BT, group.by = "Immune_ON_OFF")
Immune_ON_OFF <- GC_malignant_subset_BT@meta.data %>% dplyr::select(Immune_ON_OFF)
Immune_ON <- Immune_ON_OFF %>% subset(Immune_ON_OFF =="ON")
Immune_ON$ON_ALL <- paste("Immune", Immune_ON$Immune_ON_OFF, sep="_")


ALL_cells <- bind_rows(Neural_like_ON, Melanocytic_ON, Mesenchymal_ON, Trans_reg_ON, preEMT_ON, Immune_ON, UPR_ON) #renamed duplicate names - good then merging with metadata they will be omitted (the cells that had 2x or more ON)
dim(ALL_cells)
ALL_cells[1:5, 1:5]
ALL_cells <- ALL_cells %>% dplyr::select(ON_ALL)
ALL_cells[1:5, 1:1]
ALL_cells$x_merge <- rownames(ALL_cells)

gtest <- inner_join(ALL_cells, GC_malignant_subset_BT@meta.data, by= "x_merge")
#ON_OFF_matrix <- GC_malignant_subset_BT@meta.data[, 219:228]
#ON_OFF_matrix <- filter_all(ON_OFF_matrix, any_vars(. =="ON"))
gtest$"Mutation" <- plyr::revalue(as.character(gtest$`Mut type`),
                                  c("N/A" = "NA"
                                  ))
cell_num <- gtest %>%
  mutate(sample_id = as.factor(paste(orig.ident,Mutation,Malignant_clusters, sep="_"))) %>% # Flo here you can exchange Mutation by the origin primary or met from the Tirosh data
  mutate(ON_ALL = as.factor(ON_ALL)) %>%
  group_by(sample_id, ON_ALL, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("orig.ident", "Mutation", "Malignant_clusters"))
cell_num
total_cells<- gtest %>%
  group_by(orig.ident) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage$Mutation <- factor(cell_percentage$Mutation, levels = c("NRAS", "BRAF", "WT", "NA"))
pdf("/Users/u0128760/Documents/PROJECTS/Resource_paper_BT/Percentage_per_sample_mutation.pdf", height = 3.5, width = 10)
ggplot(cell_percentage, aes(orig.ident, y=percentage, fill=ON_ALL))+geom_bar(stat="identity", position = 'stack', width = 0.9)+RotatedAxis() +  labs(fill = "NRAS_predicted_cells") + theme(axis.title.x=element_blank()) +facet_wrap(~Mutation, ncol = 4, scales = "free") + theme_classic() + RotatedAxis()
dev.off()














