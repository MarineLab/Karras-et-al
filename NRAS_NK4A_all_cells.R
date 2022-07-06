# Florian Rambow
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
library(data.table)
library(ggraph)
library(clustree)

options(future.globals.maxSize = 10000 * 1024^2)


NRAS__1 <- readRDS("NRAS__1.rds")
NRAS__3 <- readRDS("NRAS__3.rds")
#NRAS__1 = CMA001
#NRAS__3 = CMA079 + CMA080 + CMA081 + CMA082 + CMA083 + CMA084 + CMA085 + CMA086

########### names for seurat objects ################
file_list1 <- c("NRAS__1",
                "NRAS__3")

########### This is Seurat objects list ################
file_list2 <- c(NRAS__1,
                NRAS__3)

# Create Seurat objects
for (i in 1:length(list_data)) {
assign(file_list1[i], CreateSeuratObject(counts = list_data[[i]], project = file_list1[[i]], min.cells = 10, min.features = 1000))
}

# Filter Seurat objects
for (i in 1:length(file_list2)) {
  file_list2[[i]][["percent.mt"]] <- PercentageFeatureSet(file_list2[[i]], pattern = "^mt-")
  assign(file_list1[i], subset(file_list2[[i]], subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 7.5))
}


# SCT normalization for all Seurat objects
for (i in 1:length(file_list2)) {
  assign(file_list1[i], SCTransform(file_list2[[i]], verbose = TRUE, vars.to.regress = c("percent.mt")))
}

# Doublet Finder for all Seurat objects_target8k
for (i in 1:length(file_list2)) {
  sweep.res.list <- paramSweep_v3(file_list2[[i]], PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  setwd("/Documents/")
  pdf(paste(file_list1[i], ".pdf", sep=""))
  plot(x = pK, y = BCmetric, pch = 16,type="b",
       col = "blue",lty=1)
  abline(v=pK_choose,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 3,col = "red")
  dev.off()
  nExp_poi <- round(0.061*nrow(file_list2[[i]]@meta.data))  
  nExp_poi
  assign(file_list1[i], doubletFinder_v3(file_list2[[i]], PCs = 1:10, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE))
}

##### file_list2 ################
file_list2 <- c(NRAS__1,
                NRAS__3)
#Add Doublet finder column called Doublets in each seurat object so that when you merge the objects you have one doublets column in meta data. 
#The old DF column per sample will stay in meta data with NAs for other samples
for (i in 1:length(file_list2)) {
  file_list2[[i]]$Doublets <-  file_list2[[i]]@meta.data[, grep("DF.", colnames(file_list2[[i]]@meta.data))]
  assign(file_list1[i], file_list2[[i]])
}

NRAS13merge <- merge(NRAS__1, y=c(NRAS__3), project="NRAS13merge")

slotNames(NRAS13merge)
NRAS13merge
DirRes <- "/Documents/"

pdf(file.path(DirRes, "QC.pdf"), width = 14, height = 7)
NRAS13merge[["percent.mt"]] <- PercentageFeatureSet(NRAS13merge, pattern = "^mt-")
VlnPlot(NRAS13merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(NRAS13merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NRAS13merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
NRAS13merge <- subset(NRAS13merge, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 10)
VlnPlot(NRAS13merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# run SC transformation and harmony on the merged data
NRAS13mergeH <- SCTransform(NRAS13merge, verbose = FALSE, vars.to.regress = c("percent.mt"))
NRAS13mergeH <- RunPCA(NRAS13mergeH , verbose = TRUE)
NRAS13mergeH <- RunHarmony(NRAS13mergeH, group.by.vars = "orig.ident", assay.use="SCT")
#NRAS13mergeH
#saveRDS(NRAS13mergeH , "NRAS123mergeH8k.rds")
#NRAS13mergeH<-readRDS("NRAS123mergeH8k.rds")
Idents(NRAS13mergeH ) <- 'orig.ident'

VlnPlot(NRAS13mergeH, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0.0)

#after harmony plot harmony Embeddings on a heatmap to asses after which number the variance drops
pdf(file.path(DirRes,"Harmony_Heatmap_all.pdf"), width = 7, height = 7)
harmony_embeddings <- Embeddings(NRAS13mergeH, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = TRUE,
        show_row_names = FALSE,
        name = "Harmony_embedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
        #col = col_fun,
        row_title_rot = 0)
dev.off()

# Rename identity classes
new.cluster.ids <- c()
names(new.cluster.ids) <- levels(NRAS13mergeH)
NRAS13mergeH <- RenameIdents(NRAS13mergeH, new.cluster.ids)
DimPlot(NRAS13mergeH, reduction = "umap", label = TRUE, pt.size = 0.5)
NRAS13mergeH@meta.data$seurat_clusters
WhichCells(NRAS13mergeH, idents = "stromal")
WhichCells(NRAS13mergeH, idents = "immune")

NRAS13mergeH  <- FindNeighbors(NRAS13mergeH , dims = 1:10, reduction = "harmony")
NRAS13mergeH  <- FindClusters(NRAS13mergeH , resolution = 0.5, reduction = "harmony")
NRAS13mergeH  <- RunUMAP(NRAS13mergeH , dims=1:10, reduction = "harmony")

monoc <- read.delim("43k_tumorID.txt")
monoc <- read.delim("43k_compartment.txt")
monoc <- read.delim("43K_celltype.txt")

NRAS13mergeH <- AddMetaData(NRAS13mergeH, monoc$cell_type, col.name = "cell_type")
head(monoc)


Idents(NRAS13mergeH) <- "cell_type"
### cell cycle scoring
CellCycleMouse <- readRDS('mouse_cell_cycle_genes.rds')
s.genes <- CellCycleMouse$s.genes
g2m.genes <- CellCycleMouse$g2m.genes
NRAS13mergeH  <- CellCycleScoring(NRAS13mergeH, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(NRAS13mergeH , reduction = "umap", group.by = 'Phase')

############### Subset for singlets ################
table(NRAS13mergeH$orig.ident)
NRAS13mergeH <- subset(NRAS13mergeH, subset = Doublets == "Singlet")
NRAS13mergeH
pdf(file.path(DirRes,"Fig2_CC.pdf"), width = 7, height = 7)
####################################################
################### SCT with CELL CYCLE ############
####################################################
NRAS13mergeH <- SCTransform(NRAS13mergeH, verbose = FALSE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
NRAS13mergeH <- RunPCA(NRAS13mergeH , verbose = TRUE)
NRAS13mergeH <- RunHarmony(NRAS13mergeH, group.by.vars = "orig.ident", assay.use="SCT")

pdf(file.path(DirRes,"Harmony_Heatmap_all_after_CC_reg.pdf"), width = 7, height = 7)
harmony_embeddings <- Embeddings(NRAS13mergeH, 'harmony')
harmony_embeddings[1:5, 1:5]
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(harmony_embeddings, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE,  
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        show_column_names = TRUE,
        show_row_names = FALSE,
        name = "Harmony_embedding",
        #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
        #col = col_fun,
        row_title_rot = 0)
dev.off()
NRAS13mergeH  <- FindNeighbors(NRAS13mergeH , dims = 1:20, reduction = "harmony")
NRAS13mergeH  <- FindClusters(NRAS13mergeH , resolution = 0.6, reduction = "harmony")
head(Idents(NRAS13mergeH ), 5)
NRAS13mergeH  <- RunUMAP(NRAS13mergeH , dims=1:20, reduction = "harmony")
DimPlot(NRAS13mergeH, group.by = c('seurat_clusters'), pt.size = 0.02, label=T)
DimPlot(NRAS13mergeH, group.by = c('orig.ident'), pt.size = 0.02)
DimPlot(NRAS13mergeH , reduction = "umap", group.by = 'Phase')

#markers.to.plot <- c("S100a8","Il1b","Lyz","S100a9","Cst3","Vcan","Cxcl8","Timp1","Ier3","Fth1","Nfkbia","Sod2","Ereg","Fcn1")
#DotPlot(NRAS13mergeH, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(NRAS13mergeH, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.5)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(pbmc.markers, file="NRAS1_3_all_0_22clusters.csv")

getwd()
setwd("/Documents/")
saveRDS(NRAS13mergeH,"NRAS13_all_43k.rds")

NRAS13mergeH<-readRDS("NRAS13_all_43k.rds")
############################add_CNV########################
cnv_honey <- fread("NRAS13_HoneyBadger_results.txt")
dim(cnv_honey)
cnv_honey <- as.data.frame(cnv_honey)
cnv_honey[1:5, 1:5]
cnv_honey$V1 <-  NULL
dim(cnv_honey)
class(cnv_honey)
abs_honey <- apply(cnv_honey, 2, abs)
abs_honey[1:5, 1:5]
mean_cnv <- apply(abs_honey, 2, mean)
mean_cnv <- as.data.frame(mean_cnv)

head(NRAS13mergeH@meta.data)
NRAS13mergeH@meta.data$x_merge <- rownames(NRAS13mergeH@meta.data)
mean_cnv$x_merge <- rownames(mean_cnv)
NRAS13mergeH@meta.data <- NRAS13mergeH@meta.data %>% inner_join(mean_cnv, by="x_merge")
rownames(NRAS13mergeH@meta.data) <- NRAS13mergeH@meta.data$x_merge
VlnPlot(NRAS13mergeH, features = c("mean_cnv"),pt.size = 0.0)
saveRDS(NRAS13mergeH,"NRAS13_all_43k.rds")

