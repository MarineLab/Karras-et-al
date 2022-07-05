# Florian Rambow

library(devtools)
library(BiocManager)
library(Rtsne)
library(cowplot)
library(sctransform)
library(RcppEigen)
library(uwot)
library(spatstat)
library(rsvd)
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(harmony)
#read the seurat object with all TME cells
NRAS13mergeH <- readRDS(NRAS13mergeH,"NRAS13_all_43k.rds") 

#subset for malignant cells
NRAS13mergeH_Mstr <- subset(NRAS13mergeH_M, subset = M_TI > 0.11 | mean_cnv > 0.13) 
NRAS13mergeH_Mstr <- subset(NRAS13mergeH_Mstr , subset = superMEL > 0.16 | mean_cnv > 0.13) 
NRAS13mergeH_Mstr <- subset(NRAS13mergeH_Mstr, subset = Ptprc < 0.0001)
#run SC transformation and harmony
NRAS13mergeH_Mstr <- SCTransform(NRAS13mergeH_Mstr, verbose = TRUE, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'))
NRAS13mergeH_Mstr <- RunPCA(NRAS13mergeH_Mstr , verbose = TRUE)
NRAS13mergeH_Mstr <- RunHarmony(NRAS13mergeH_Mstr, group.by.vars = "orig.ident", assay.use="SCT")

#after harmony plot harmony Embeddings on a heatmap to asses after which number the variance drops
harmony_embeddings <- Embeddings(NRAS13mergeH_Mstr, 'harmony')
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

NRAS13mergeH_Mstr  <- FindNeighbors(NRAS13mergeH_Mstr , dims = 1:15 , reduction = "harmony")
NRAS13mergeH_Mstr  <- FindClusters(NRAS13mergeH_Mstr, save.snn=T , resolution = 0.3)
head(Idents(NRAS13mergeH_Mstr), 5)
NRAS13mergeH_Mstr  <- RunUMAP(NRAS13mergeH_Mstr , dims=1:15, reduction = "harmony")
DimPlot(NRAS13mergeH_Mstr, group.by = c ('seurat_clusters'),label=T)
VlnPlot(NRAS13mergeH_Mstr, features = c("mean_cnv"),pt.size = 0.0)


# find markers for every cluster compared to all remaining cells, report only the positive ones
NRAS13mergeH_Mstr.markers <- FindAllMarkers(NRAS13mergeH_Mstr, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
NRAS13mergeH_Mstr.markers %>% group_by(cluster) %>% top_n(n = 9, wt = avg_logFC)

write.csv(NRAS13mergeH_Mstr.markers, file="NRAS13mergeH_Mstr_7clusters.csv")

top10 <- NRAS13mergeH_Mstr.markers %>% group_by (cluster) %>% top_n(n = 9, wt = avg_logFC)
DoHeatmap(NRAS13mergeH_Mstr, features = top10$gene)



library(cluster)
#silhouette scores
cell_distance<- dist(NRAS13mergeH_Mstr@reductions$harmony@cell.embeddings[, 1:15])
cell_cluster<- as.numeric(as.character(Idents(NRAS13mergeH_Mstr)))
silhouette_score<- cluster::silhouette(cell_cluster, cell_distance)
head(silhouette_score)
silhouette <- silhouette(as.numeric(cell_cluster), dist = cell_distance)

NRAS13mergeH_Mstr@meta.data$silhouette_score <- silhouette[,3]

mean_silhouette_score <- mean(NRAS13mergeH_Mstr@meta.data$silhouette_score)



p <- NRAS13mergeH_Mstr@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  arrange(seurat_clusters,-silhouette_score) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  ggplot() +
  geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = TRUE) +
  geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
  scale_x_discrete(name = 'Cells') +
  scale_y_continuous(name = 'Silhouette score') +
  # scale_fill_manual(values = custom_colors$discrete) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p



