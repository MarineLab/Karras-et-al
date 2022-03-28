# By Joanna Pozniak
library("HoneyBADGER")
library("Rsamtools")
library("GenomicFeatures")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("edgeR")
library("DESeq2")
library("ComplexHeatmap")
#xfun::pkg_load2(c('base64enc', 'htmltools', 'mime'))
library(dplyr)
library(devtools)
library(Seurat)
library(Matrix)
library(AUCell)
library(dplyr)
library(stringr)
library(GSEABase)
library(GSA)
library(Seurat)
library(circlize)
Takis10x <- readRDS("/Users/u0128760/Documents/PROJECTS/Takis_10x_NRAS_INK4A/NRAS13_all_43k.rds")
Takis10x

################################################################# for HoneyBADGER ##############################
genes4 <- c("Acap1","Akna","Alox5Ap","Ankrd44","Apobec3G","Arhgap15","Arhgap25","Arhgap30","Arhgap4","Arhgap9","Arhgdib","Atp2A3","Bin2","C16Orf54","Ccdc88B","Cd37","Cd48","Cd52","Cd53","Cd69","Cd84","Cdc42Se2","Celf2","Cntrl","Coro1A","Csk","Cxcr4","Cyth4","Cytip","Def6","Dennd1C","Dock2","Dock8","Dusp2","Evi2B","Fermt3","Fgd3","Fnbp1","Gbp5","Gpr65","Gpsm3","Hcls1","Hmha1","Ikzf1","Il10Ra","Il16","Il2Rg","Inpp5D","Itga4","Itgal","Itgb2","Lair1","Laptm5","Lcp1","Lilrb3","Limd2","Lpxn","Lsp1","Ly9","Map4K1","Myo1G","Nckap1L","Nr4A2","Parp8","Parvg","Pik3Cd","Pim2","Plcb2","Plekha2","Prkcb","Psd4","Pstpip2","Ptk2B","Ptpn22","Ptpn6","Ptpn7","Ptprc","Rac2","Rassf5","Rcsd1","Rgs1","Rhoh","Rps6Ka1","Samsn1","Sash3","Sla","Snx20","Sp140","Stk17B","Tagap","Tbc1D10C","Tmc6","Tmc8","Tmsb4X","Traf3Ip3","Tsc22D3","Tstd1","Ucp2","Vav1","Wipf1")
geneSets <- GeneSet(genes4, setName="Immune")
geneSets
cells_rankings <- AUCell_buildRankings(Takis10x@assays[["RNA"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Immune<-getAUC(cells_AUC)
Immune<-t(Immune)
Takis10x@meta.data<-cbind(Takis10x@meta.data, Immune)

FeaturePlot(Takis10x , features = c ("Immune"))
### annotation ###
origin <- Takis10x@meta.data[["orig.ident"]]
immune <- Takis10x@meta.data[["Immune"]]
malignant <- Takis10x@meta.data[["Lineage"]]
ha = HeatmapAnnotation(origin = origin,
                       immune = immune,
                       malignant = malignant,
                       annotation_name_side = "left")

### subsetting for immune only to use as normal for HoneyBADGER ####
#Takis10x_normal <- subset(x=Takis10x, subset = malign_call =='non_malig')
#Takis10x_normal@assays[["RNA"]]@counts@Dimnames[[1]] <- str_to_title(Takis10x_normal@assays[["RNA"]]@counts@Dimnames[[1]])
#Takis10x_normal <- as.matrix(Takis10x_normal@assays[["RNA"]]@counts)
### subsetting for immune only to use as normal for HoneyBADGER ####
Immune_sub <- subset(Takis10x, subset = Immune > 0.1)
Immune_sub
Normal_immune <- as.matrix(GetAssayData(Immune_sub@assays[["RNA"]]))


# Load mart object for later annotation
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'mmusculus_gene_ensembl')

# Load the expression data as a gene x cell matrix
#Takis10x@assays[["RNA"]]@counts@Dimnames[[1]] <- str_to_title(Takis10x@assays[["RNA"]]@counts@Dimnames[[1]])

exprs	<- as.matrix(Takis10x@assays[["RNA"]]@counts)
##colnames(exprs) <- gsub(".141218.*", "", colnames(exprs))

# Clean up... and normalize
##rownames(exprs)	<- exprs$ENSEMBL_ID
##exprs		<- as.matrix(exprs[,-1])
exprs 		<- cpm(exprs, log=TRUE)

# Load a control expression matrix of more or less normal expression. This one was found in the GTex database... Make sure the IDs in both datasets are the same.
exprsNorm <- Normal_immune

# Clean up... and normalize
##rownames(exprsNorm) 	<- exprsNorm$ID_REF
##exprsNorm 		<- as.matrix(exprsNorm[,-1])
exprsNorm 		<- cpm(exprsNorm, log=TRUE)

# Only retain the genes shared in both datasets
sharedGenes 	<- intersect(rownames(exprsNorm), rownames(exprs))
exprs 		<- exprs[sharedGenes,]
exprsNorm 	<-exprsNorm[sharedGenes, ]

# Keep 5000 most expressed genes (other options are possible here...)
highExpr 	<- names(head(sort(rowMeans(exprs), decreasing=TRUE), 5000))
exprs 		<- exprs[highExpr,]
exprsNorm 	<-exprsNorm[highExpr, ]
# Take average ref
exprsNorm <- rowMeans(exprsNorm)

# Run honeybadger. Change to id to fit your data.
hb <- new('HoneyBADGER', name='Takis10x')

hb$setGexpMats(exprs, exprsNorm, mart.obj, filter = FALSE, scale=TRUE, verbose=TRUE, id = "mgi_symbol")
png("/Users/u0128760/Documents/PROJECTS/Takis_10x_NRAS_INK4A/HoneyBadger/InitViz7500Scaled_immune0_12.png")
hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = FALSE, returnPlot = FALSE)
gexpPlot <- hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = TRUE, returnPlot = FALSE)
dev.off()

# Run honeybadger again and create prettier plot with complexheatmap.
gexpPlot <- hb$plotGexpProfile(gexp.norm.sub = hb$gexp.norm, setOrder = FALSE, returnPlot = TRUE)
plotData <- do.call(rbind, gexpPlot)

write.table(plotData, file = "/Users/u0128760/Documents/PROJECTS/Takis_10x_NRAS_INK4A/HoneyBadger/HoneyBadger_results.txt", quote=F, sep="\t", row.names = T)
saveRDS(gexpPlot, file = "/Users/u0128760/Documents/PROJECTS/Takis_10x_NRAS_INK4A/HoneyBadger/HoneyBadger_gexplot.rds")


pdf("/Users/u0128760/Documents/PROJECTS/Takis_10x/HoneyBadger/Final_heatmap.pdf")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht<-Heatmap(plotData, 
            cluster_rows = FALSE, 
            cluster_columns = TRUE,  
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "complete", 
            show_column_names = FALSE,
            #heatmap_legend_param = list(legend_direction = "vertical", title_position = "leftcenter"),
            gap = unit(2, "mm"), 
            top_annotation = ha,
            col = col_fun,
            split = factor(rep(c(paste("Chr", c(1:19))), times = sapply(gexpPlot, nrow)), levels = c(paste("Chr", c(1:19)))), 
            row_title_rot = 0)
draw(ht)
dev.off()

