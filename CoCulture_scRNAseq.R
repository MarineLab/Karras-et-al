# By Joanna Pozniak
library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(DoubletFinder)
library(harmony)
library(nichenetr)
Ada_0 <- Read10X(data.dir = "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Ada_0h/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
Ada_0 <- CreateSeuratObject(counts = Ada_0, project = "Time_0h", min.cells = 10, min.features = 500)
Ada_0
Ada_0[["percent.mt"]] <- PercentageFeatureSet(Ada_0, pattern = "^mt-")
Ada_0 <- SCTransform(Ada_0, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(Ada_0, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results")
pdf("DF_Ada_0.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(Ada_0@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi
Ada_0 <- doubletFinder_v3(Ada_0, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
Ada_0$Doublets <-  Ada_0@meta.data[, grep("DF.", colnames(Ada_0@meta.data))]



Ada_48h <- Read10X(data.dir = "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Ada_48h/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
Ada_48h <- CreateSeuratObject(counts = Ada_48h, project = "Time_48h", min.cells = 10, min.features = 500)
Ada_48h
Ada_48h[["percent.mt"]] <- PercentageFeatureSet(Ada_48h, pattern = "^mt-")
Ada_48h <- SCTransform(Ada_48h, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(Ada_48h, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results")
pdf("DF_Ada_48h.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(Ada_48h@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi
Ada_48h <- doubletFinder_v3(Ada_48h, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
Ada_48h$Doublets <-  Ada_48h@meta.data[, grep("DF.", colnames(Ada_48h@meta.data))]



Ada_coculture <- Read10X(data.dir = "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Ada_coculture/raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
Ada_coculture <- CreateSeuratObject(counts = Ada_coculture, project = "Coculture", min.cells = 10, min.features = 500)
Ada_coculture
Ada_coculture[["percent.mt"]] <- PercentageFeatureSet(Ada_coculture, pattern = "^mt-")
Ada_coculture <- SCTransform(Ada_coculture, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(Ada_coculture, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results")
pdf("DF_Ada_coculture.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(Ada_coculture@meta.data))  ## Assuming 3.9% doublet formation rate - tailor for your dataset
nExp_poi
Ada_coculture <- doubletFinder_v3(Ada_coculture, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
Ada_coculture$Doublets <-  Ada_coculture@meta.data[, grep("DF.", colnames(Ada_coculture@meta.data))]
table(Ada_coculture$Doublets)

#DefaultAssay(Ada_0) <- "RNA"
#DefaultAssay(Ada_48h) <- "RNA"
#DefaultAssay(Ada_coculture) <- "RNA"


All_merged <- merge(Ada_0, y=c(Ada_48h, Ada_coculture), project="All_merged")
All_merged


#saveRDS(All_merged,"/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/All_merged.rds")
######################################################################### From here only baseline vs co-culture ###############################
All_merged <- readRDS("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results_all_samples/All_merged.rds")
All_merged <- subset(All_merged, subset = orig.ident != "Time_48h")

table(All_merged$orig.ident)
DirRes <- "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results"

pdf(file.path(DirRes, "QC.pdf"), width = 14, height = 7)
All_merged[["percent.mt"]] <- PercentageFeatureSet(All_merged, pattern = "^mt-")
VlnPlot(All_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
plot1 <- FeatureScatter(All_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(All_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf(file.path(DirRes, "QC_subset.pdf"), width = 14, height = 7)
All_merged <- subset(All_merged, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 7)
VlnPlot( All_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()
table(All_merged_subset$orig.ident)

pdf(file.path(DirRes ,"Initial_doublets.pdf"))
All_merged<- SCTransform(All_merged, vars.to.regress = c("percent.mt"), verbose = FALSE)
All_merged<- RunPCA(All_merged, features = VariableFeatures(object = All_merged))
VizDimLoadings(All_merged, dims = 1:2, reduction = "pca")
DimPlot(All_merged, reduction = "pca")
DimHeatmap(All_merged, dims = 1:20, cells = 500, balanced = TRUE)
All_merged<- JackStraw(All_merged, num.replicate = 100)
All_merged<- ScoreJackStraw(All_merged, dims = 1:20)
JackStrawPlot(All_merged, dims = 1:20)
ElbowPlot(All_merged)
All_merged<- FindNeighbors(All_merged, dims = 1:20)
All_merged <- FindClusters(All_merged, resolution = 0.4)
head(Idents(All_merged), 5)
All_merged<- RunUMAP(All_merged, dims=1:20)
DimPlot(All_merged, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(All_merged, reduction = "umap", group.by="Doublets")
DimPlot(All_merged, reduction = "umap", group.by="orig.ident")
All_merged
mouse_cell_cycle_genes <- readRDS("/Users/u0128760/Documents/PROJECTS/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
All_merged  <- CellCycleScoring(All_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(All_merged, reduction = "umap", group.by = 'Phase')
dev.off()

############### Subset for singlets ################
All_merged_subset <- subset(All_merged, subset = Doublets == "Singlet")
All_merged_subset <- SCTransform(All_merged_subset, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
All_merged_subset<- RunPCA(All_merged_subset, features = VariableFeatures(object = All_merged_subset))


All_merged_subset  <- FindNeighbors(All_merged_subset , dims = 1:20)
All_merged_subset  <- FindClusters(All_merged_subset , resolution = 0.2)
All_merged_subset  <- RunUMAP(All_merged_subset , dims=1:20)
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/UMAPs_including_EC.pdf")
DimPlot(All_merged_subset, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(All_merged_subset, reduction = "umap", group.by="orig.ident")
FeaturePlot(All_merged_subset, c("Pecam1", "Sox10", "S100a1"))
dev.off()

All_merged_subset <- subset(All_merged_subset, subset = Pecam1 <0.0001)

DimPlot(All_merged_subset)
####################################################
################### SCT with CELL CYCLE ############
####################################################
All_merged_subset <- SCTransform(All_merged_subset, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score', 'orig.ident'))
All_merged_subset<- RunPCA(All_merged_subset, features = VariableFeatures(object = All_merged_subset))
All_merged_subset  <- FindNeighbors(All_merged_subset , dims = 1:20)
All_merged_subset  <- FindClusters(All_merged_subset , resolution = 0.2)
All_merged_subset  <- RunUMAP(All_merged_subset , dims=1:20)

pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/UMAPs.pdf")
DimPlot(All_merged_subset, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(All_merged_subset, reduction = "umap", group.by="orig.ident")
dev.off()
#saveRDS(All_merged_subset, "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/All_merged_subset_no_pecam1.rds")

markers <- FindAllMarkers(All_merged_subset, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3) # 0.3, 0.4 before
write.table(markers, "/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Heatmap.pdf", width = 10, height = 13)
DoHeatmap(All_merged_subset, features = top20$gene) + NoLegend()
dev.off()


######################### ON/OFF all states
#################################### Karras et all ON OFF
################ Signatures from Landscape top 50 ###########################################
##########################        Neural_like        ##########################
All_merged_subset$Neural_like <- NULL
genes<-c("Mest","Wnt4","Fn1","Akap12","Sema5a","Trf","Slc29a1","Csn3","Kctd12","Emilin1","Postn","Sema3d","Dhh","Igf1","Moxd1","Lmcd1","Cd200","Fibin","Igfbp4","Aqp1","Qpct","Thsd7a","Mgp","Lbhd2","Timp1","Cavin2","Plvap","Enpp2","Serpina3n","Spry4","Gja1","Ltbp1","Tmem37","Tmem158","Tm4sf1","Gsn","Egfl8","Sulf2","Fxyd5","Col11a1","Chl1","Ephx1","Bpgm","Spon1","Abcg2","Sdc4","Fth1","Itgb5","Cxxc4","Egflam")
geneSets <- GeneSet(genes, setName="Neural_like")
geneSets
#cells_rankings <- AUCell_buildRankings(All_merged_subset@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Neural_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Neural_like<-getAUC(cells_AUC)
Neural_like<-t(Neural_like)
All_merged_subset@meta.data<-cbind(All_merged_subset@meta.data, Neural_like)
FeaturePlot(All_merged_subset, features = "Neural_like", label = T) +NoAxes()
VlnPlot(All_merged_subset, features = "Neural_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(All_merged_subset, features = "Neural_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(All_merged_subset, features = "Neural_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
dev.off()

##########################        Melanocytic_OXPHOS        ##########################
All_merged_subset$Melanocytic_OXPHOS <- NULL
genes<-c("Mlana","Ptgds","Car2","Dct","Pmel","Slc45a2","Car6","Syngr1","Gstp1","Cst6","Gpnmb","Chchd10","Lgals3","Car14","Sparcl1","Hpse","Gsta4","Mcoln3","Aebp1","Fxyd3","Cox6c","Uba52","Cox7c","Lgals1","Cox6a1","Gm2115","Ndufb2","Cox7a2","Apoe","Enho","Gjb6","Cyb5a","Cited1","Hpgds","Ndufa1","Rps28","Scrg1","Cox4i1","2010107E04Rik","Cox7b","Cd63","Kcnj10","Pla2g2e","Uqcc2","Cox6b1","S100b","Cryab","Ppia","Rpl35","Uqcr11")
geneSets <- GeneSet(genes, setName="Melanocytic_OXPHOS")
geneSets
#cells_rankings <- AUCell_buildRankings(All_merged_subset@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Melanocytic_OXPHOS.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Melanocytic_OXPHOS<-getAUC(cells_AUC)
Melanocytic_OXPHOS<-t(Melanocytic_OXPHOS)
All_merged_subset@meta.data<-cbind(All_merged_subset@meta.data, Melanocytic_OXPHOS)
FeaturePlot(All_merged_subset, features = "Melanocytic_OXPHOS", label = T)
VlnPlot(All_merged_subset, features = "Melanocytic_OXPHOS", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(All_merged_subset, features = "Melanocytic_OXPHOS", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(All_merged_subset, features = "Melanocytic_OXPHOS", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
dev.off()



##########################        Stem_like         ##########################
All_merged_subset$Stem_like <- NULL

genes<-c("Pltp","Serpine2","Vcan","Nes","Gpx3","Mdm2","Siva1","Btc","Carhsp1","1700007K13Rik","Cald1","Dcxr","Itga6","Celf5","Pam","Sdc1","Man2b1","Espn","Pmm1","Cdkn1a","Ak1","Kdm7a","Zfp385a","Mybl1","Ahnak2","Notch3","Ccng1","Plxdc2","Dnajc9","Cep170b","Il11","Crip2","Palm3","Ctxn1","Pde3b","Angptl2","Ptprt","Agpat4","Ift27","Rnase4","Gtse1","Fas","Creb3l1","Axl","Ahi1","Llph","Akr1b10","Ercc5","Dusp15","Fosl1")
geneSets <- GeneSet(genes, setName="Stem_like")
geneSets
#cells_rankings <- AUCell_buildRankings(All_merged_subset@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Stem_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stem_like<-getAUC(cells_AUC)
Stem_like<-t(Stem_like)
All_merged_subset@meta.data<-cbind(All_merged_subset@meta.data, Stem_like)
FeaturePlot(All_merged_subset, features = "Stem_like", label = T) +NoAxes()
VlnPlot(All_merged_subset, features = "Stem_like", pt.size = 0) + theme(legend.position = 'none')
VlnPlot(All_merged_subset, features = "Stem_like", pt.size = 0, group.by = "orig.ident") + theme(legend.position = 'none')
dev.off()




######################## stem on off
All_merged_subset$STEM_ON_OFF <- ifelse(All_merged_subset$Stem_like > 0.08,  "ON", "OFF")
DimPlot(All_merged_subset, group.by = "STEM_ON_OFF")
All_merged_subset@meta.data$"Sample" <- plyr::revalue(as.character(All_merged_subset$orig.ident),
                                                      c(  "Coculture" = "Coculture",
                                                          "Time_0h" = "Control" ))

TEST <- All_merged_subset@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(Sample, STEM_ON_OFF,  sep="_"))) %>%
  #mutate(Mutation = as.factor(Mutation)) %>%
  group_by(sample_id, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("Sample", "STEM_ON_OFF")) 
cell_num
total_cells<- TEST %>%
  group_by(Sample) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage <- subset(cell_percentage, subset = STEM_ON_OFF =="ON")
cell_percentage
cell_percentage$Sample <- factor(cell_percentage$Sample, levels = c("Control", "Coculture"))
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Stem_signature_ON_OFF.pdf", width = 3, height = 3)
ggbarplot(cell_percentage, x = "Sample", y = "percentage",shape = "STEM_ON_OFF", fill = "STEM_ON_OFF", title = "preEMT signature ON") +NoLegend()
dev.off()


################## neural signature

All_merged_subset$Neural_like_ON_OFF <- ifelse(All_merged_subset$Neural_like > 0.13,  "ON", "OFF")
DimPlot(All_merged_subset, group.by = "Neural_like_ON_OFF")
TEST <- All_merged_subset@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(Sample, Neural_like_ON_OFF,  sep="_"))) %>%
  #mutate(Mutation = as.factor(Mutation)) %>%
  group_by(sample_id, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("Sample", "Neural_like_ON_OFF")) 
cell_num
total_cells<- TEST %>%
  group_by(Sample) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage <- subset(cell_percentage, subset = Neural_like_ON_OFF =="ON")
cell_percentage
cell_percentage$Sample <- factor(cell_percentage$Sample, levels = c("Control", "Coculture"))

pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Neural_signature_ON_OFF.pdf", width = 3, height = 3)
ggbarplot(cell_percentage, x = "Sample", y = "percentage",shape = "Neural_like_ON_OFF", fill = "Neural_like_ON_OFF", title = "Neural-like signature ON") +NoLegend()
dev.off()

#Melanocytic

All_merged_subset$Melanocytic_ON_OFF <- ifelse(All_merged_subset$Melanocytic_OXPHOS > 0.3,  "ON", "OFF")
DimPlot(All_merged_subset, group.by = "Melanocytic_ON_OFF")

TEST <- All_merged_subset@meta.data
cell_num <- TEST %>%
  mutate(sample_id = as.factor(paste(Sample, Melanocytic_ON_OFF,  sep="_"))) %>%
  #mutate(Mutation = as.factor(Mutation)) %>%
  group_by(sample_id, .drop=FALSE) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::separate(sample_id, c("Sample", "Melanocytic_ON_OFF")) 
cell_num
total_cells<- TEST %>%
  group_by(Sample) %>%
  dplyr::summarise(total = n())
total_cells
cell_percentage<- left_join(cell_num, total_cells) %>%
  mutate(percentage = n/total*100)
cell_percentage
cell_percentage <- subset(cell_percentage, subset = Melanocytic_ON_OFF =="ON")
cell_percentage
cell_percentage$Sample <- factor(cell_percentage$Sample, levels = c("Control", "Coculture"))
pdf("/Users/u0128760/Documents/PROJECTS/Endothelial_expr/Results/Melanocytic_signature_ON_OFF.pdf", width = 3, height = 3)
ggbarplot(cell_percentage, x = "Sample", y = "percentage",shape = "Melanocytic_ON_OFF", fill = "Melanocytic_ON_OFF", title = "Melanocytic signature ON") +NoLegend()
dev.off()


