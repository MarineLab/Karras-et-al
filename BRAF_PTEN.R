# By Joanna Pozniak

library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(GSEABase)
library(GSA)
library(dplyr)
library(MAST)
library(future)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

########### read the output from CellRanger ################
sc5rCMA298_mm10 <- Read10X(data.dir = "/Volumes/Samsung_1/BRAF_PTEN_heterogeneity/sc5rCMA298_mm10/mapped_tomato/sc5rCMA298_mm10/outs/raw_feature_bc_matrix")
sc5rCMA298_mm10 <- CreateSeuratObject(counts = sc5rCMA298_mm10, project = "sc5rCMA298_mm10", min.cells = 10, min.features = 500)
sc5rCMA298_mm10
sc5rCMA298_mm10[["percent.mt"]] <- PercentageFeatureSet(sc5rCMA298_mm10, pattern = "^mt-")
sc5rCMA298_mm10 <- SCTransform(sc5rCMA298_mm10, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(sc5rCMA298_mm10, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500")
pdf("DF_sc5rCMA298_mm10.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(sc5rCMA298_mm10@meta.data))  
nExp_poi
sc5rCMA298_mm10 <- doubletFinder_v3(sc5rCMA298_mm10, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
sc5rCMA298_mm10$Doublets <-  sc5rCMA298_mm10@meta.data[, grep("DF.", colnames(sc5rCMA298_mm10@meta.data))]


sc5rCMA299_mm10 <- Read10X(data.dir = "/Volumes/Samsung_1/BRAF_PTEN_heterogeneity/sc5rCMA299_mm10/mapped_tomato/sc5rCMA299_mm10/outs/raw_feature_bc_matrix")
sc5rCMA299_mm10 <- CreateSeuratObject(counts = sc5rCMA299_mm10, project = "Primary", min.cells = 10, min.features = 500) #primary is an arbitrary name 
sc5rCMA299_mm10
sc5rCMA299_mm10[["percent.mt"]] <- PercentageFeatureSet(sc5rCMA299_mm10, pattern = "^mt-")
sc5rCMA299_mm10 <- SCTransform(sc5rCMA299_mm10, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(sc5rCMA299_mm10, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500")
pdf("DF_sc5rCMA299_mm10.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.048*nrow(sc5rCMA299_mm10@meta.data))  
nExp_poi
sc5rCMA299_mm10 <- doubletFinder_v3(sc5rCMA299_mm10, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
sc5rCMA299_mm10$Doublets <-  sc5rCMA299_mm10@meta.data[, grep("DF.", colnames(sc5rCMA299_mm10@meta.data))]




CMA332 <- Read10X(data.dir = "/Volumes/Samsung_1/BRAF_PTEN_heterogeneity/CMA332/raw_feature_bc_matrix")
CMA332 <- CreateSeuratObject(counts = CMA332, project = "CMA332", min.cells = 10, min.features = 500)
CMA332
CMA332[["percent.mt"]] <- PercentageFeatureSet(CMA332, pattern = "^mt-")
CMA332 <- SCTransform(CMA332, verbose = TRUE, vars.to.regress = c("percent.mt"), return.only.var.genes = F)
sweep.res.list <- paramSweep_v3(CMA332, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
setwd("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500")
pdf("DF_CMA332.pdf")
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
dev.off()
nExp_poi <- round(0.079*nrow(CMA332@meta.data))  
nExp_poi
CMA332 <- doubletFinder_v3(CMA332, PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
CMA332$Doublets <-  CMA332@meta.data[, grep("DF.", colnames(CMA332@meta.data))]


Merged_1 <- merge(sc5rCMA298_mm10, y=c(CMA332), project="BRAF_PTEN")
Merged_1
saveRDS(Merged_1,"/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Merged1.rds")

Merged_1@meta.data$"orig.ident" <- plyr::revalue(as.character(Merged_1$orig.ident),
                                                 c("Primary" = "sc5rCMA298"
                                                 ))

table(Merged_1$orig.ident)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/QC.pdf", width = 14, height = 7)
VlnPlot(Merged_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
plot1 <- FeatureScatter(Merged_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Merged_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

Merged <- subset(Merged_1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)


mouse_cell_cycle_genes <- readRDS("/Users/u0128760/Documents/PROJECTS/mouse_cell_cycle_genes.rds")

pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Initial_doublets.pdf")
Merged<- SCTransform(Merged, vars.to.regress = c("percent.mt", "orig.ident"), verbose = FALSE)
Merged<- RunPCA(Merged, features = VariableFeatures(object = Merged))
VizDimLoadings(Merged, dims = 1:2, reduction = "pca")
DimPlot(Merged, reduction = "pca")
DimHeatmap(Merged, dims = 1:20, cells = 500, balanced = TRUE)
Merged<- FindNeighbors(Merged, dims = 1:20)
Merged <- FindClusters(Merged, resolution = 0.4)
head(Idents(Merged), 5)
Merged<- RunUMAP(Merged, dims=1:20)
DimPlot(Merged, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(Merged, reduction = "umap", group.by="Doublets")
DimPlot(Merged, reduction = "umap", group.by="orig.ident")
Merged
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
Merged  <- CellCycleScoring(Merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(Merged, reduction = "umap", group.by = 'Phase')
dev.off()
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Initial_singlet.pdf")
############### Subset for singlets ################
Merged_subset <- subset(Merged, subset = Doublets == "Singlet")
Merged_subset
table(Merged$orig.ident)
####################################################
################### SCT with CELL CYCLE ############
####################################################
Merged_subset <- SCTransform(Merged_subset, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
Merged_subset<- RunPCA(Merged_subset, features = VariableFeatures(object = Merged))
VizDimLoadings(Merged_subset, dims = 1:2, reduction = "pca")
DimPlot(Merged_subset, reduction = "pca")
DimHeatmap(Merged_subset, dims = 1:20, cells = 500, balanced = TRUE)
Merged_subset<- FindNeighbors(Merged_subset, dims = 1:20)
Merged_subset <- FindClusters(Merged_subset, resolution = 0.4)
head(Idents(Merged_subset), 5)
Merged_subset<- RunUMAP(Merged_subset, dims=1:20)
DimPlot(Merged_subset, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(Merged_subset, reduction = "umap", group.by="Doublets")
DimPlot(Merged_subset, reduction = "umap", group.by="orig.ident")
dev.off()
merged_markers <- FindAllMarkers(Merged_subset, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(merged_markers, "/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/markers.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- merged_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/heatmap_markers.pdf", width = 15, height = 30)
DoHeatmap(Merged_subset, features = top20$gene, group.by='seurat_clusters')
dev.off()
table(Merged_subset$orig.ident)
################################harmony
Merged_subset <- SCTransform(Merged_subset, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
Merged_subset <- RunPCA(Merged_subset , verbose = TRUE)
Merged_subset <- RunHarmony(Merged_subset, group.by.vars = "orig.ident", assay.use="SCT")

#after harmony, plot harmony Embeddings on a heatmap to asses after which number the variance drops
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/harmony_heatmap.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(Merged_subset, 'harmony')
harmony_embeddings[1:5, 1:5]
#col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = TRUE,
        show_row_names = FALSE,
        name = "Hramony_embeedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
        #col = col_fun,
        row_title_rot = 0)
dev.off()
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Initial_singlet_harmony.pdf")

Merged_subset  <- FindNeighbors(Merged_subset , dims = 1:8, reduction = "harmony")
Merged_subset  <- FindClusters(Merged_subset , resolution = 0.3, reduction = "harmony")
Merged_subset  <- RunUMAP(Merged_subset , dims=1:8, reduction = "harmony") 

DimPlot(Merged_subset, reduction = "umap", label = T, group.by = "seurat_clusters")
DimPlot(Merged_subset, reduction = "umap", group.by="Doublets")
DimPlot(Merged_subset, reduction = "umap", group.by="orig.ident")
dev.off()
#saveRDS(Merged_subset, "/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Braf_Pten.rds")

Merged_subset <- readRDS("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Braf_Pten.rds")

pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Braf_Pten_Tirosh_supermel.pdf")
###################Tirosh
genes<-c("Mia","Tyr","Slc45A2","Cdh19","Pmel","Slc24A5","Magea6","Gjb1","Plp1","Prame","Capn3","Erbb3","Gpm6B","S100B","Fxyd3","Pax3","S100A1","Mlana","Slc26A2","Gpr143","Cspg4","Sox10","Mlph","Loxl4","Plekhb1","Rab38","Qpct","Birc7","Mfi2","Linc00473","Sema3B","Serpina3","Pir","Mitf","St6Galnac2","Ropn1B","Cdh1","Abcb5","Qdpr","Serpine2","Atp1A1","St3Gal4","Cdk2","Acsl3","Nt5Dc3","Igsf8","Mbp")
geneSets <- GeneSet(genes, setName="Tirosh_malignant")
geneSets
cells_rankings <- AUCell_buildRankings(Merged_subset@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tirosh_malignant<-getAUC(cells_AUC)
Tirosh_malignant<-t(Tirosh_malignant)
Merged_subset@meta.data<-cbind(Merged_subset@meta.data, Tirosh_malignant)
FeaturePlot(Merged_subset, features = "Tirosh_malignant")
###################SuperMEL

genes<-c("Copg2","Cd59a","Nceh1","Cdh19","Gjc3","Cers4","Sort1","Plekhb1","Pax3","Sox10","Rapgef4","Kcnn4","Akr1b7","Syngr1")
geneSets <- GeneSet(genes, setName="SuperMEL")
geneSets
#cells_rankings <- AUCell_buildRankings(Merged_subset@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
SuperMEL<-getAUC(cells_AUC)
SuperMEL<-t(SuperMEL)
Merged_subset@meta.data<-cbind(Merged_subset@meta.data, SuperMEL)
FeaturePlot(Merged_subset, features = c("SuperMEL","Tirosh_malignant", "tdTomato"))
VlnPlot(Merged_subset, features = c("SuperMEL","Tirosh_malignant", "tdTomato"))
FeatureScatter(Merged_subset, "SuperMEL", "Notch3")
dev.off()



#################################Malignant subset
Mailg_braf_pten <- subset(Merged_subset, subset = SuperMEL >0.1 | tdTomato > 2)
DimPlot(Mailg_braf_pten)
Mailg_braf_pten <- SCTransform(Mailg_braf_pten, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
Mailg_braf_pten<- RunPCA(Mailg_braf_pten) 
Mailg_braf_pten <- RunHarmony(Mailg_braf_pten, group.by.vars = "orig.ident", assay.use="SCT")

#after harmony plot harmony Embeddings on a heatmap to asses after which number the variance drops
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/harmony_heatmap_malig.pdf", width = 7, height = 7)
harmony_embeddings <- Embeddings(Mailg_braf_pten, 'harmony')
harmony_embeddings[1:5, 1:5]
#col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = TRUE,
        show_row_names = FALSE,
        name = "Hramony_embeedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
        #col = col_fun,
        row_title_rot = 0)
dev.off()
Mailg_braf_pten  <- FindNeighbors(Mailg_braf_pten , dims = 1:7, reduction = "harmony")
Mailg_braf_pten  <- FindClusters(Mailg_braf_pten , resolution = 0.4, reduction = "harmony")
Mailg_braf_pten  <- RunUMAP(Mailg_braf_pten , dims=1:7, reduction = "harmony") 
DimPlot(Mailg_braf_pten, reduction = "umap", group.by = "seurat_clusters") +NoAxes()

pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/UMAP_clusters.pdf", width = 4, height = 4)
DimPlot(Mailg_braf_pten, reduction = "umap", group.by = "seurat_clusters") +NoAxes()
DimPlot(Mailg_braf_pten, reduction = "umap", group.by="orig.ident")
dev.off()

pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/MITF_SOX10_MALANA_Prrx1_clusters.pdf", width = 4, height = 4)
FeaturePlot(Mailg_braf_pten, features = c("Mitf"))
FeaturePlot(Mailg_braf_pten, features = c("Sox10"))
FeaturePlot(Mailg_braf_pten, features = c("Mlana"))
FeaturePlot(Mailg_braf_pten, features = c("Prrx1"))
dev.off()

Idents(Mailg_braf_pten) <- "seurat_clusters"
merged_markers <- FindAllMarkers(Mailg_braf_pten, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.table(merged_markers, "/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/markers_malig.txt", sep='\t', quote = FALSE, col.names = T, row.names = F)
top20 <- merged_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/heatmap_markers_malig.pdf", width = 15, height = 15)
DoHeatmap(Mailg_braf_pten, features = top20$gene, group.by='seurat_clusters')
dev.off()

##########################        Verfaillie_PRO        ##########################
genes<-c("Gdpd5",  "Slc22A23",  "Rp11-3L8.3",  "Pdzrn3",  "Nr4A1",  "Pros1",  "Ss18L1",  "Ap1S2",  "Tspan10",  "Bhlhe41",  "Ccdc171",  "Il17D",  "Pdk4",  "C2Orf88",  "Prr5",  "Cdk2",  "Ctsh",  "Mreg",  "Fam124A",  "Rp11-558F24.4",  "Asah1",  "Hey2",  "Sdc3",  "Rap1Gap",  "Capg",  "Rab6B",  "Sort1",  "Mxi1",  "Pcsk9",  "St3Gal4",  "Apoc1",  "Adamts17",  "Epha5",  "Hey1",  "Slc12A7",  "Sh3Tc1",  "Pag1",  "Rab38",  "Rp11-1055B8.2",  "Scin",  "Linc00340",  "Gab2",  "Ces3",  "Dusp15",  "Cyp27A1",  "Atrnl1",  "Faxc",  "Hes6",  "Klf15",  "Dgcr5",  "Ac005786.5",  "Inpp5F",  "Sgk1",  "Npas1",  "Cygb",  "Pir",  "Cers1",  "Aatk",  "Sesn3",  "Cdk5R1",  "Scube2",  "Slc17A9",  "Slc16A10",  "Tmtc2",  "Kiaa1598",  "Rassf2",  "Mast1",  "Rp11-80F22.9",  "Tbc1D7",  "Prkch",  "Pde3B",  "C19Orf71",  "Slain1",  "Gja3",  "Fam53B",  "Rab11Fip4",  "Bai1",  "C11Orf96",  "Slc27A3",  "Fgf13",  "Brsk2",  "Egln3",  "Gnal",  "Cables1",  "Gpr137B",  "Cxadr",  "Shc2",  "St6Galnac1",  "Fbxl16",  "Z83851.1",  "Asrgl1",  "Tnfrsf19",  "Ceacam1",  "Sorl1",  "Ankrd6",  "Isg20",  "Rims4",  "Myom2",  "Lad1",  "Adrbk2",  "Lzts1",  "Rnf125",  "Tmem255A",  "Trpm8",  "Fam20A",  "Lonrf3",  "Ppm1H",  "Sptbn2",  "Tex41",  "St6Gal1",  "Pou3F3",  "Adam23",  "Ano4",  "Mfsd12",  "Rp11-390P2.4",  "St6Galnac2",  "Rp11-137H2.6",  "Olfm2",  "Tmcc2",  "Greb1",  "Ttc39A",  "Fam213A",  "Kcns1",  "Tnfrsf14",  "Stxbp6",  "Itga7",  "Aldh1A1",  "Znf704",  "Bambi",  "Pgbd5",  "Prkcz",  "Il6R",  "Plcl1",  "Egr3",  "Itpkb",  "Nat16",  "Lrrc4",  "Stox2",  "Ogdhl",  "Pik3Ap1",  "Pnliprp3",  "Cntn3",  "Baat",  "Col25A1",  "Celf2",  "Rasip1",  "Tmem229B",  "Plekhg1",  "Pknox2",  "Krtap19-1",  "Slc7A4",  "Slc24A4",  "Asb4",  "St3Gal6-As1",  "Efr3B",  "Nkain1",  "Caskin1",  "Ldlrad4",  "Wnk2",  "Tktl1",  "Rab17",  "Kndc1",  "Tesk2",  "Chn2",  "Slc7A8",  "Fgd4",  "Crtac1",  "Ppargc1A",  "Rp11-527H14.2",  "Fam134B",  "Rasef",  "Acan",  "Chst6",  "Tfcp2L1",  "Hmcn1",  "Rp3-395M20.8",  "Pnmal1",  "Rp3-527G5.1",  "Pip5K1B",  "Il12Rb2",  "Tenm1",  "Rps6Ka2",  "Cecr2",  "Vgf",  "Bcan",  "Adcy1",  "Rab3C",  "Cldn14",  "Rp11-481A20.11",  "Mertk",  "Linc00937",  "Linc00504",  "Ccl18",  "Rxrg",  "Phactr1",  "Card14",  "Qpct", 
         "Grasp",  "Dll3",  "Pklr",  "Lrguk",  "C1Orf51",  "Tex15",  "B4Galnt3",  "Kiaa1211",  "Peli2",  "Pou3F2",  "Prkcb",  "Hcg20",  "Rp11-98L5.2",  "Disc1Fp1",  "Rp11-317M11.1",  "Kbtbd11",  "Lamc3",  "Rp11-557H15.4",  "Mpz",  "Ac011294.3",  "Renbp",  "Prune2",  "Lama1",  "B3Gat1",  "Tincr",  "Ndn",  "Ttyh2",  "Ctb-151G24.1",  "Tc2N",  "Slc16A6",  "Rp11-143A12.3",  "Linc00518",  "Mfi2",  "Prodh",  "Golga7B",  "Pou3F1",  "Cdh3",  "Gpm6A",  "Nr4A3",  "Mapt",  "Large",  "Lrp2",  "Lpl",  "Tmprss5",  "Rlbp1",  "Gng7",  "Acp5",  "Rp11-93B14.5",  "Fam167B",  "Ednrb",  "Tmc6",  "Linc00426",  "Ca8",  "Myo16",  "St3Gal6",  "Cited1",  "Rasgef1A",  "Ac009784.3",  "Dapk1",  "Kit",  "Gyg2",  "Plekhb1",  "Ap000479.1",  "Tmprss13",  "Nup210",  "Aldh1A2",  "Znf536",  "Fam19A5",  "Rgs1",  "Nkx2-5",  "Oplah",  "Tspan7",  "Kdr",  "Plekhh1",  "Sbk1",  "Fam155B",  "Itga9",  "Best1",  "Kcnab2",  "Rp11-2E17.1",  "Mageb2",  "Rp13-735L24.1",  "Tubb8P7",  "Cpvl",  "Dennd1C",  "Mob3B",  "St8Sia6",  "Wipf3",  "Prss33",  "Cntn1",  "Ccdc64",  "Apod",  "Nrg3",  "Rp3-395M20.7",  "Ac009499.1",  "Ppp1R14C",  "Capn3",  "Sox6",  "Mbp",  "Ism1",  "Mitf",  "Cacna1H",  "Rp4-718J7.4",  "Mgat4A",  "Cpb2-As1",  "Lrrc4B",  "Extl1",  "Sox8",  "Cryab",  "Ropn1B",  "Syt3",  "Sgca",  "Nat8L",  "Rp11-509E16.1",  "Pmp2",  "Nmrk2",  "Ac004988.1",  "Gstt1",  "Linc00589",  "Rp11-189B4.6",  "Tyrp1",  "Fam189A2",  "Aldh3B2",  "Rp11-161M6.2",  "Rab33A",  "Zdhhc11B",  "Robo2",  "Slc35F1",  "Rragd",  "Lingo1",  "Rp11-290F20.3",  "Hspb8",  "Rp11-669N7.2",  "Rp3-332B22.1",  "Cacna1D",  "Samd5",  "Gfpt2",  "Itgax",  "Cpn1",  "Rp11-726G1.1",  "Msi1",  "Tnrc6C-As1",  "Lgi3",  "Mlip",  "Tubb4A",  "Cobl",  "Lhfpl3-As1",  "Linc00488",  "Rp11-557H15.2",  "Galnt3",  "S100B",  "Ac002511.1",  "Sorbs1",  "Igf1",  "Adcy2",  "Dlgap1",  "Sftpc",  "Hils1",  "Plxnc1",  "Myh14",  "Ac145110.1",  "Pla1A",  "Glb1L2",  "Chl1",  "Gapdhs",  "Ctnna2",  "Rp11-599J14.2",  "She",  "Ctd-2207A17.1",  "Rp11-1055B8.3",  "Rp11-104E19.1",  "Ac096559.1",  "Fxyd3",  "Hpgd",  "Mmp8",  "Irx6",  "Rp11-347E10.1",  "Gjb1",  "Gpr143",  "Plp1",  "Cdh19",  "Il16",  "Mpped2",  "Rp11-252C15.1",  "Cdh1",  "C10Orf90",  "Apoe",  "Itih5",  "Oca2",  
         "Maf", "Fam69C", "Rp4-529N6.1",  "Sox10",  "Frmd4B",  "Erbb3",  "Sema6A",  "Mcf2L",  "Mapk4",  "Fam174B",  "Slc38A8",  "Myo1D",  "Linc00520",  "Nsg1",  "Igsf11",  "C20Orf26",  "Hrk",  "Paep",  "Slc45A2",  "Trim51",  "Atp10A",  "Ropn1",  "Col9A3",  "Scml4",  "Rp11-429E11.2",  "Gas7",  "Ca14",  "Linc00698",  "Ptprz1",  "Dct",  "Mid2",  "Pmel",  "Sorcs1",  "Abcb5",  "Prdm7",  "Gpm6B",  "Lcp2",  "Enthd1",  "Mlana",  "Tyr",  "Nkain4",  "Irf4",  "Sgcd",  "Slc24A5",  "Trim63",  "Birc7",  "Trpm1")
geneSets <- GeneSet(genes, setName="Verfaillie_PRO")
geneSets
cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Verfaillie_PRO<-getAUC(cells_AUC)
Verfaillie_PRO<-t(Verfaillie_PRO)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Verfaillie_PRO)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Verfaillie_PRO_vln.pdf", width = 5, height = 4)
FeaturePlot(Mailg_braf_pten, features = "Verfaillie_PRO")
VlnPlot(Mailg_braf_pten, features = "Verfaillie_PRO", pt.size = 0) + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Mitf.pdf", width = 5, height = 4)
FeaturePlot(Mailg_braf_pten, features = "Mitf")
dev.off()

##########################        Verfaillie_INV        ##########################
genes<-c("Mxra5",  "Gabre",  "Ddx3Y",  "Vcam1",  "Postn",  "Cd248",  "Abcc3",  "Anpep",  "Tmem200A",  "Adamts12",  "F2Rl2",  "Il1B",  "Rps4Y1",  "Ereg",  "Col1A1",  "Serpinb7",  "Adamtsl1",  "Kdm5D",  "Txlng2P",  "Gpr128",  "Alpk2",  "Slc14A1",  "Vipr1",  "Cdh13",  "Ptgfr",  "Ptx3",  "Cpa4",  "Serpine1",  "Ccl2",  "Ccdc68",  "Neto1",  "Edn1",  "Eltd1",  "Plau",  "Rp11-48O20.4",  "Cdh6",  "Eif1Ay",  "Rgs4",  "Col8A1",  "Gpr39",  "Fendrr",  "Esm1",  "Cta-392C11.1",  "Dpep3",  "Sstr1",  "Prky",  "Fpr3",  "Smek3P",  "Tnfrsf9",  "Pxdn",  "Fpr1",  "Magea4", 
         "Fam196B",  "Col1A2",  "Anxa10",  "Mkx",  "Wnt7B",  "Armc4",  "Fgf5",  "Muc13",  "Uty",  "Galnt5",  "Gdf5",  "Lrrc15",  "Kal1",  "Spanxd",  "Lrrc17",  "Smoc1",  "Ntng1",  "Bdnf",  "Ifi27",  "Zfy",  "Pappa",  "Loxl2",  "Il1A",  "Usp9Y",  "Tspyl5",  "Spanxc",  "Ccbe1",  "Lypd6B",  "Pappa2",  "Adamts6",  "C7Orf69",  "Cpa3",  "Igfn1",  "C3",  "G0S2",  "Trhde",  "Itga11",  "Cda",  "Ctd-2171N6.1",  "Abcc9",  "Apbb1Ip",  "Il6",  "Fam43B",  "Itgbl1",  "Creb3L1",  "Ttty15",  "Foxf1",  "Foxr2",  "Psg5",  "Vgll3",  "Foxg1",  "Tnfsf18",  "Rab27B",  "Ido1",  
         "Fam180A",  "Parm1",  "Rab3B",  "Krt81",  "Meox2",  "Grem1",  "Ecscr",  "Galnt6",  "Il32",  "Cd163L1",  "Trhde-As1",  "C6Orf141",  "Mmp19",  "Thbd",  "Ldb2",  "Pbx1",  "Gpr68",  "Bdkrb2",  "Tenm2",  "Cfi",  "Edil3",  "Linc00707",  "Mir137Hg",  "Gbp1",  "F2Rl1",  "Pde1C",  "Rtn1",  "Axl",  "Abi3Bp",  "Triml2",  "Nexn",  "Ca9",  "Tcf4",  "Nid2",  "Myct1",  "Hs3St3A1",  "Clmp",  "Slc15A3",  "Gstm1",  "Sh2D4A",  "Col13A1",  "Glis3",  "Bnc1",  "Pou2F2",  "Spata18",  "Col5A1",  "Aox1",  "S1Pr1",  "Scn9A",  "Mt1E",  "Efemp1",  "Glipr1",  "Gbp4",  "Rp11-371I1.2",  
         "Oas2",  "Nmnat2",  "Rp11-224O19.2",  "Rsad2",  "Chmp4C",  "Fam155A",  "Pitx1",  "Cfh",  "Gabrq",  "Trim58",  "Prdm8",  "Loxl1",  "Apol3",  "Arhgdib",  "Jph2",  "Krt80",  "Emilin1",  "Tnfrsf11B",  "Htr1F",  "Rrad",  "Il18R1",  "Il7R",  "Amigo2",  "Atp8B1",  "Fbn2",  "Cxcl2",  "Nrg1",  "Isg15",  "Hspb6",  "Inhba",  "Rasgrf2",  "Ntn4",  "Gucy1B3",  "Il4I1",  "Col12A1",  "Mt1A",  "Il11",  "Tmem158",  "Stk33",  "Arntl2",  "Pkp2",  "Lamc2",  "Tor4A",  "Sod3",  "Tmem119",  "Srgn",  "Igfbp6",  "Scg2",  "Il1Rl1",  "Tox2",  "Nlrp3",  "Birc3",  "Hspb7",  "Fbn1",  "Sulf1", 
         "Kcnma1",  "Proser2-As1",  "Ssc5D",  "Dennd2A",  "Pdcd1Lg2",  "S1Pr3",  "Col6A3",  "Lpar1",  "Vit",  "Col3A1",  "Ptprn",  "Dnah2",  "Dkk3",  "Fgf1",  "Xaf1",  "Hecw2",  "Rp11-166D19.1",  "Cfb",  "Flnc",  "Ace",  "Cyp1B1",  "Galnt13",  "Slc6A10P",  "Col5A2",  "Thbs1",  "Serpinb2",  "Clic2",  "Mpp4",  "Klhl4",  "Ca12",  "Rarres3",  "Vegfc",  "Col4A6",  "Wnt5A",  "Hhipl2",  "Oasl",  "Cldn11",  "Csf1R",  "Mypn",  "S100A4",  "Tgm2",  "Slit2",  "Tnfaip2",  "Krt18",  "Ccl5",  "Ifi44L",  "Magea10",  "Mx2",  "Nrp1",  "Ccdc80",  "Magea11",  "Gabra2",  "Hrh1",  "Ifitm1",  "Nr2F1",  
         "Trpc4",  "Tgfb2",  "Cdk15",  "Znf385D",  "Krt8",  "Rp5-1011O1.3",  "Oxtr",  "Ngef",  "Lox",  "Fmod",  "Pparg",  "Rp11-513O17.2",  "Srpx2",  "Irak2",  "Mgst1",  "Crispld2",  "Ar",  "Ephb2",  "Lama3",  "Cdh12",  "Trim22",  "Rp11-400K9.4",  "Ptplad2",  "Arhgap22",  "Chrdl1",  "Slfn11",  "Adam19",  "Colec12",  "Stac",  "Myl9",  "Ntsr1",  "Tslp",  "Pdgfrb",  "Stc2",  "Pcdhga7",  "Gpr1",  "Slc16A9",  "Egfr",  "Spocd1",  "Frmpd4",  "Akr1C3",  "Prtfdc1",  "Clu",  "Fam101A",  "Plagl1",  "Il8",  "Cacna1C",  "Dse",  "Fam84A",  "Slc16A2",  "Col6A2",  "Fzd2",  "Efnb2",  "Loxl1-As1",  
         "Artn",  "Stc1",  "Ltbp1",  "Mecom",  "Znf804A",  "Rarres1",  "Reln",  "Tmem200B",  "Laptm5",  "Mmp2",  "Camk1D",  "Pdgfc",  "Podxl",  "Gcnt1",  "Fgf2",  "Cd96",  "Cyr61",  "Hic1",  "Mpp7",  "A4Galt",  "Linc00960",  "Layn",  "Wisp1",  "Mx1",  "Nptx1",  "Wnt5A-As1",  "Wnt5B",  "Ac005789.11",  "Fam65C",  "Irf1",  "Samd9L",  "Kcnt2",  "Foxe1",  "Nhs",  "Ly6K",  "Foxc2",  "Crim1",  "Dock3",  "Tnfrsf10A",  "Dock2",  "Cyp27C1",  "Amz1",  "Sectm1",  "Sox9",  "Fibcd1",  "Sh3Rf2",  "Dcn",  "Scara3",  "Blk",  "Tagln",  "Pcolce",  "Thsd4",  "Kiaa1462",  "Dkk1",  "Fam167A",  "Pcdhgb4",  
         "Samd12",  "B3Gnt9",  "Ifi44",  "Dpyd",  "Tle4",  "Il20Rb",  "Fam20C",  "Pdgfb",  "Syt1",  "C12Orf75",  "Scg5",  "Epha2",  "Adam8",  "Ampd3",  "Ltbp2",  "Magec2",  "Eps8L2",  "Ryr2",  "Pcdhb3",  "Nnmt",  "Vasn",  "Speg",  "Lphn2",  "Pcdhga11",  "Plaur",  "Itgb4",  "Itga2",  "Cxcl3",  "Adamts4",  "Prdm1",  "Steap2",  "Basp1",  "Spred3",  "S100A16",  "Tdrd9",  "Mylk",  "Col6A1",  "Angptl4",  "Pcdh10",  "Tfpi",  "Pdzd2",  "Rp11-575F12.2",  "Susd1",  "Htr7",  "Tpm2",  "Cpe",  "Flnb",  "Znf415",  "Gpx1P1",  "Sybu",  "C1R",  "Errfi1",  "Igfbp3",  "Dpf3",  "Znf185",  "Cxcl1",  "Pip5Kl1",  "Gpr176",  "Trabd2A",  "Nr2F1-As1",  "Ifi6",  "Vcan",  "Shc3",  "Scn2A",  "Nuak2",  "Itga5",  "Spock1",  "Ikzf2",  "Ddah1",  "Bcl3",  "Tfpi2",  "Inhbe",  "Rras",  "Osmr",  "Lamb3",  "Efemp2",  "Nmi",  "Ifit1")
geneSets <- GeneSet(genes, setName="Verfaillie_INV")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Verfaillie_INV<-getAUC(cells_AUC)
Verfaillie_INV<-t(Verfaillie_INV)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Verfaillie_INV)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Verfaillie_INV_vln.pdf", width = 5, height = 4)
FeaturePlot(Mailg_braf_pten, features = "Verfaillie_INV")
VlnPlot(Mailg_braf_pten, features = "Verfaillie_INV", pt.size = 0) + theme(plot.margin=unit(c(1,1,1.5,1.5),"cm"), legend.position = 'none')
dev.off()

Mailg_braf_pten <- readRDS("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Mailg_braf_pten.rds")

################ Signatures from NRAS/INK4A mouse  ###########################################
##########################        Neural_like        ##########################
Mailg_braf_pten$Neural_like <- NULL
genes<-c("Tfap2b","Prnp","Mef2c","Gfra1","Cd200","Syt11","Thsd7a","Cxxc4","Sema5a","Tbx3","Kif26b","Efhd1","Neto2","Hmcn1","Igsf10","Olfml3","Rgs2","Nt5e","Morc4","Aqp1","Gfra2","Wnt4","Abcg2","Elovl5","Emilin1","Fibin")
geneSets <- GeneSet(genes, setName="Neural_like")
geneSets
cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Neural_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Neural_like<-getAUC(cells_AUC)
Neural_like<-t(Neural_like)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Neural_like)
FeaturePlot(Mailg_braf_pten, features = "Neural_like", label = T) +NoAxes()
VlnPlot(Mailg_braf_pten, features = "Neural_like", pt.size = 0) + theme(legend.position = 'none')
dev.off()
Mailg_braf_pten$Neural_like_ON_OFF <- ifelse(Mailg_braf_pten$Neural_like > 0.13,  "ON", "OFF")
DimPlot(Mailg_braf_pten, group.by = "Neural_like_ON_OFF", cols = c('grey', 'blue'), pt.size = 1.5) +NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Neural_like_ON_OFF.pdf",useDingbats=FALSE,family="sans")


##########################        Melanocytic_OXPHOS        ##########################
genes<-c("Mlph","Mitf","Mlana","Pmel","Slc45a2","Apoe","Dct","Gpnmb","Tyr","Cited1","Bace2","Cox17","Cox7a2","Ndufb2","Ndufb4","Uqcr10","Uqcrb")
geneSets <- GeneSet(genes, setName="Melanocytic_OXPHOS")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Melanocytic_OXPHOS.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Melanocytic_OXPHOS<-getAUC(cells_AUC)
Melanocytic_OXPHOS<-t(Melanocytic_OXPHOS)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Melanocytic_OXPHOS)
FeaturePlot(Mailg_braf_pten, features = "Melanocytic_OXPHOS", label = T)
VlnPlot(Mailg_braf_pten, features = "Melanocytic_OXPHOS") + theme(legend.position = 'none')
dev.off()
Mailg_braf_pten$Melanocytic_OXPHOS_ON_OFF <- ifelse(Mailg_braf_pten$Melanocytic_OXPHOS > 0.34,  "ON", "OFF")
table(Mailg_braf_pten$Melanocytic_OXPHOS_ON_OFF )
DimPlot(Mailg_braf_pten, group.by = "Melanocytic_OXPHOS_ON_OFF", cols = c('grey', 'blue'), pt.size = 1.5) +NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Melanocytic_OXPHOS_ON_OFF.pdf",useDingbats=FALSE,family="sans")


##########################        Trans_reg        ##########################
genes<-c("Pabpn1","Ccnl1","Pnisr","Sfpq","Ccnl2","Bclaf1","Hnrnpu","Rbm25","Luc7l2","Ddx17","Cpsf6","Snrnp70","Srek1","Prpf4b","Srrm2","Rbm39","Son","Ddx5","Prpf38b","Pan3","Tra2a","Hnrnph1","Fus")
geneSets <- GeneSet(genes, setName="Trans_reg")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Trans_reg.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Trans_reg<-getAUC(cells_AUC)
Trans_reg<-t(Trans_reg)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Trans_reg)
FeaturePlot(Mailg_braf_pten, features = "Trans_reg", label = T)
VlnPlot(Mailg_braf_pten, features = "Trans_reg") + theme(legend.position = 'none')
dev.off()


##########################        Stress_hypoxia        ##########################
genes<-c("Bnip3","Tpi1","Slc2a1","Mif","Vldlr","Hk2","Vegfa","Ldha","Pfkl","P4ha1","Fam162a","Bhlhe40","Pgk1","Aldoa","Pfkp","Pdk1","Hspa9","Pgam1","Pkm","Atf4")
geneSets <- GeneSet(genes, setName="Stress_hypoxia")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Stress_hypoxia.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stress_hypoxia<-getAUC(cells_AUC)
Stress_hypoxia<-t(Stress_hypoxia)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Stress_hypoxia)
FeaturePlot(Mailg_braf_pten, features = "Stress_hypoxia", label = T)
VlnPlot(Mailg_braf_pten, features = "Stress_hypoxia") + theme(legend.position = 'none')
dev.off()


##########################        Stem_like         ##########################
genes<-c("Gtse1","Cep170b","Tap1","Ercc5","Kank3","Rap2a","Ei24","Zmat3","Fas","Igdcc4","Notch3","Vcan","Mybl1","Dusp15","Dgka","Nrp1","Slc27a3","Slc19a2","Dgki","Siva1","Rbm38","Cad","Mcm7","Nme4","Thyn1","Gpx3","Arl4c","Ass1","Efnb2","Frrs1","Sdc1","Itga6","Icam1","Chst3","Nes","Sox4","Fat1","Angptl2")
geneSets <- GeneSet(genes, setName="Stem_like")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Stem_like.pdf", width = 4, height = 4)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Stem_like<-getAUC(cells_AUC)
Stem_like<-t(Stem_like)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Stem_like)
FeaturePlot(Mailg_braf_pten, features = "Stem_like", label = T) +NoAxes()
FeaturePlot(Mailg_braf_pten, features = "Nes", label = T)
FeaturePlot(Mailg_braf_pten, features = "Thra", label = T)
VlnPlot(Mailg_braf_pten, features = "Stem_like", pt.size = 0) + theme(legend.position = 'none')
dev.off()

Mailg_braf_pten$Stem_like_ON_OFF <- ifelse(Mailg_braf_pten$Stem_like > 0.04,  "ON", "OFF")
DimPlot(Mailg_braf_pten, group.by = "Stem_like_ON_OFF", cols = c('grey', 'blue'), pt.size = 1.5) +NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Stem_like_ON_OFF.pdf",useDingbats=FALSE,family="sans")


##########################        Immune        ##########################
genes<-c("Gbp4","Gbp6","Irf1","Nlrc5","Stat2","Cd74","Il12rb1","B2m","Gbp2","Ifit3","Gbp3","Psmb8","Stat1","Ifit1","Isg15","Psmb9","Gbp7","Xaf1","Tap1","Gbp5","Usp18","Ciita","Irf7","Ube2l6","Tap2","Psme1","Ifitm3","Psme2","Tapbp","Ifi35","Eif2ak2","Irf9","Ddx58","Trim25","Ifit2","Irf8","H2-DMa","H2-Aa","H2-DMb1","H2-Eb1","H2-Ab1","H2-T23")
library(nichenetr)
#genes1 <- convert_mouse_to_human_symbols(genes)
#genes1 <- unique(genes1)
geneSets <- GeneSet(genes, setName="Immune")
geneSets
cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Immune.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Immune)
FeaturePlot(Mailg_braf_pten, features = "Immune", label = T)
VlnPlot(Mailg_braf_pten, features = "Immune") + theme(legend.position = 'none')
dev.off()


##########################        Mesenchymal        ##########################
genes<-c("Fap","Lama2","Prrx1","Slit3","Abi3bp","Loxl1","Cdh11","Col6a3","Col6a2","Loxl2","Mfap5","Fbn1","Col4a2","Pcolce","Lum","Col4a1","Col5a2","Thy1","Fbn2","Bgn","Tgfbi","Pdgfrb","Sulf1","Inhba","Col1a1","Col3a1","Col1a2","Ctsk","Dcn","Serpinf1","Sparc","Fstl1")
geneSets <- GeneSet(genes, setName="Mesenchymal")
geneSets
#cells_rankings <- AUCell_buildRankings(Mailg_braf_pten@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
pdf("/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Mesenchymal.pdf", width = 7.08, height = 5.8)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Mesenchymal<-getAUC(cells_AUC)
Mesenchymal<-t(Mesenchymal)
Mailg_braf_pten@meta.data<-cbind(Mailg_braf_pten@meta.data, Mesenchymal)
FeaturePlot(Mailg_braf_pten, features = "Mesenchymal", label = T)
VlnPlot(Mailg_braf_pten, features = "Mesenchymal") + theme(legend.position = 'none')
dev.off()
hist(Mailg_braf_pten$Mesenchymal)
Mailg_braf_pten$Mesenchymal_ON_OFF <- ifelse(Mailg_braf_pten$Mesenchymal > 0.067,  "ON", "OFF")
table(Mailg_braf_pten$Mesenchymal_ON_OFF)
DimPlot(Mailg_braf_pten, group.by = "Mesenchymal_ON_OFF", cols = c('grey', 'blue'), pt.size = 1.5) +NoAxes()
dev.copy2pdf(file="/Users/u0128760/Documents/PROJECTS/BRAF_PTEN_heterogeneity/Results_500/Mesenchymal_ON_OFF.pdf",useDingbats=FALSE,family="sans")


