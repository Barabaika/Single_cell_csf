library(Matrix)
library(Seurat)
library(dplyr)

setwd ("/mnt/Data10tb/student9_data/GSE138266_RAW")
####loading file 
data.dir = "/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_3"

CSF_MS_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_1")
CSF_MS_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_2")
CSF_MS_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_3")
CSF_MS_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_4")
CSF_MS_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_5")
CSF_MS_6.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_6")
CSF_cont_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_1")
CSF_cont_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_2")
CSF_cont_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_3")
CSF_cont_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_4")
CSF_cont_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_5")
CSF_cont_6.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_6")

CSF_MS_1.data <- CreateSeuratObject(counts = CSF_MS_1.data, min.cells = 3, min.features = 200, project = "CSF_MS_1")
CSF_MS_2.data <- CreateSeuratObject(counts = CSF_MS_2.data, min.cells = 3, min.features = 200, project = "CSF_MS_2")
CSF_MS_3.data <- CreateSeuratObject(counts = CSF_MS_3.data, min.cells = 3, min.features = 200, project = "CSF_MS_3")
CSF_MS_4.data <- CreateSeuratObject(counts = CSF_MS_4.data, min.cells = 3, min.features = 200, project = "CSF_MS_4")
CSF_MS_5.data <- CreateSeuratObject(counts = CSF_MS_5.data, min.cells = 3, min.features = 200, project = "CSF_MS_5")
CSF_MS_6.data <- CreateSeuratObject(counts = CSF_MS_6.data, min.cells = 3, min.features = 200, project = "CSF_MS_6")

CSF_cont_1.data <- CreateSeuratObject(counts = CSF_cont_1.data, min.cells = 3, min.features = 200, project = "CSF_cont_1")
CSF_cont_2.data <- CreateSeuratObject(counts = CSF_cont_2.data, min.cells = 3, min.features = 200, project = "CSF_cont_2")
CSF_cont_3.data <- CreateSeuratObject(counts = CSF_cont_3.data, min.cells = 3, min.features = 200, project = "CSF_cont_3")
CSF_cont_4.data <- CreateSeuratObject(counts = CSF_cont_4.data, min.cells = 3, min.features = 200, project = "CSF_cont_4")
CSF_cont_5.data <- CreateSeuratObject(counts = CSF_cont_5.data, min.cells = 3, min.features = 200, project = "CSF_cont_5")
CSF_cont_6.data <- CreateSeuratObject(counts = CSF_cont_6.data, min.cells = 3, min.features = 200, project = "CSF_cont_6")

CSF <- merge(CSF_MS_1.data, c(CSF_MS_2.data, CSF_MS_3.data, CSF_MS_4.data, CSF_MS_5.data, CSF_MS_6.data, CSF_cont_1.data, CSF_cont_2.data, CSF_cont_3.data, CSF_cont_4.data, CSF_cont_5.data, CSF_cont_6.data), 
              add.cell.ids = c("CSF_MS_1", "CSF_MS_2", "CSF_MS_3", "CSF_MS_4", "CSF_MS_5", "CSF_MS_6.data", 
                               "CSF_cont_1", "CSF_cont_2", "CSF_cont_3", "CSF_cont_4", "CSF_cont_5", "CSF_cont_6.data"))


CSF.list <- SplitObject(CSF)

CSF.list$CSF_MS_1$stim <- "MS"
CSF.list$CSF_MS_2$stim <- "MS"
CSF.list$CSF_MS_3$stim <- "MS"
CSF.list$CSF_MS_4$stim <- "MS"
CSF.list$CSF_MS_5$stim <- "MS"
CSF.list$CSF_MS_6$stim <- "MS"

CSF.list$CSF_cont_1$stim <- "CONT"
CSF.list$CSF_cont_2$stim <- "CONT"
CSF.list$CSF_cont_3$stim <- "CONT"
CSF.list$CSF_cont_4$stim <- "CONT"
CSF.list$CSF_cont_5$stim <- "CONT"
CSF.list$CSF_cont_6$stim <- "CONT"

rm(CSF_MS_1.data, CSF_MS_2.data, CSF_MS_3.data, CSF_MS_4.data, 
   CSF_MS_5.data,CSF_MS_6.data, CSF_cont_1.data, CSF_cont_2.data, 
   CSF_cont_3.data, CSF_cont_4.data, CSF_cont_5.data, CSF_cont_6.data)
# run garbage collect to free up memory
gc()


for (i in 1:length(CSF.list)) {
  # filtering
  CSF.list[[i]][["percent.mt"]] <- PercentageFeatureSet(CSF.list[[i]], pattern = "^MT-")
  CSF.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(CSF.list[[i]], pattern = "^RP[SL]")
  CSF.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(CSF.list[[i]], pattern = "^HB[^(P)]")
  
  CSF.list[[i]] <- subset(CSF.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  CSF.list[[i]] <- SCTransform(CSF.list[[i]], verbose = TRUE)
}


options(future.globals.maxSize= 14000000000)

CSF.features <- SelectIntegrationFeatures(object.list = CSF.list, nfeatures = 3000)

CSF.list <- PrepSCTIntegration(object.list = CSF.list, anchor.features = CSF.features, 
                                verbose = FALSE)

CSF.anchors <- FindIntegrationAnchors(object.list = CSF.list, normalization.method = "SCT", 
                                       anchor.features = CSF.features, verbose = TRUE)
CSF.integrated <- IntegrateData(anchorset = CSF.anchors, normalization.method = "SCT", 
                                 verbose = FALSE) 

# ctrl$stim <- "CTRL"

#clastering
CSF.integrated <- RunPCA(CSF.integrated, verbose = FALSE)
CSF.integrated <- RunUMAP(CSF.integrated, dims = 1:30)

CSF.integrated <- FindNeighbors(CSF.integrated, reduction = "pca", dims = 1:30)
CSF.integrated <- FindClusters(CSF.integrated, resolution = 0.5)

p1 <- DimPlot(CSF.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(CSF.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# find markers for every cluster compared to all remaining cells, report only the positive ones
CSF.markers <- FindAllMarkers(CSF.integrated, min.pct = 0.4, logfc.threshold = 0.25, base = exp(1))

top_CSF.markers <- CSF.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(CSF.integrated, features = top_chf_ms.markers$gene)

write.csv(CSF.markers,file="/home/student9/projects/SC_annot_scSA/SCSA/csf_markers.csv",quote=FALSE)
#annot:через терминал
# python3 SCSA.py -d whole.db -s seurat -i seurat_GSE72056.csv -k All -E -g Human -p 0.01 -f 1.5

#annotatig cell types to clasters 

CSF.cluster.ids <- c('Hepatocyte','Regulatory T', 'Th17|NK cell', 'NK T cell', 'T cell', 'Regulatory T', 'Regulatory T|T cell', 'T helper1 (Th1) cell', 
                     'Macrophage|Monocyte', 'Germ cell|Primordial germ cell', 'Myeloid dendritic cell|Monocyte', 'B cell', 'NK cell', 'Astrocyte|NK T cell',
                     'B cell', 'Monocyte|Macrophage', 'B cell|Memory B cell', 'Plasmacytoid dendritic cell')
names(CSF.cluster.ids) <- levels(CSF.integrated)
CSF.integrated <- RenameIdents(CSF.integrated, CSF.cluster.ids)

DimPlot(CSF.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DimPlot(CSF.integrated, reduction = "umap", split.by = "stim")

#загружаем результат GO, полученный через scSA 
csf_go <- read.table(header = TRUE, '/mnt/Data10tb/student9_data/GSE138266_RAW/results_go_csf.go', sep = ";" )

# diff exp for each claster ms vs cont
CSF.integrated$celltype.stim <- paste(Idents(CSF.integrated), CSF.integrated$stim, sep = "_")
CSF.integrated$celltype <- Idents(CSF.integrated)
Idents(CSF.integrated) <- "celltype.stim"
# 1 cluster
Hepatocyte.markers <- FindMarkers(CSF.integrated, ident.1 = "Hepatocyte_MS", ident.2 = "Hepatocyte_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 2 cluster
Tregs.markers <- FindMarkers(CSF.integrated, ident.1 = "Regulatory T_MS", ident.2 = "Regulatory T_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 3 cluster
Th17_NK.markers <- FindMarkers(CSF.integrated, ident.1 = "Th17|NK cell_MS", ident.2 = "Th17|NK cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 4 cluster
Tcell.markers <- FindMarkers(CSF.integrated, ident.1 = "T cell_MS", ident.2 = "T cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 5 cluster
Tregs_Tcell.markers <- FindMarkers(CSF.integrated, ident.1 = "Regulatory T|T cell_MS", ident.2 = "Regulatory T|T cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 6 cluster
Th1.markers <- FindMarkers(CSF.integrated, ident.1 = "T helper1 (Th1) cell_MS", ident.2 = "T helper1 (Th1) cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 7 cluster
Macro_mono.markers <- FindMarkers(CSF.integrated, ident.1 = "Macrophage|Monocyte_MS", ident.2 = "Macrophage|Monocyte_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 8 cluster
Germ.markers <- FindMarkers(CSF.integrated, ident.1 = "Germ cell|Primordial germ cell_MS", ident.2 = "Germ cell|Primordial germ cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 9 cluster
Myeloid_dend_mono.markers <- FindMarkers(CSF.integrated, ident.1 = "Myeloid dendritic cell|Monocyte_MS", ident.2 = "Myeloid dendritic cell|Monocyte_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 10 cluster
Astrocyte_NK.markers <- FindMarkers(CSF.integrated, ident.1 = "Astrocyte|NK T cell_MS", ident.2 = "Astrocyte|NK T cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 11 cluster
Mono_Macro.markers <- FindMarkers(CSF.integrated, ident.1 = "Monocyte|Macrophage_MS", ident.2 = "Monocyte|Macrophage_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 12 cluster
Bcell.markers <- FindMarkers(CSF.integrated, ident.1 = "B cell_MS", ident.2 = "B cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 13 cluster
NK.markers <- FindMarkers(CSF.integrated, ident.1 = "NK cell_MS", ident.2 = "NK cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 14 cluster
Bcell_MemoryB.markers <- FindMarkers(CSF.integrated, ident.1 = "B cell|Memory B cell_MS", ident.2 = "B cell|Memory B cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 15 cluster
Plasma.markers <- FindMarkers(CSF.integrated, ident.1 = "Plasmacytoid dendritic cell_MS", ident.2 = "Plasmacytoid dendritic cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)


# markers = merge(Hepatocyte.markers, c(Tregs.markers, Th17_NK.markers, Tcell.markers, Tregs_Tcell.markers, Th1.markers, Macro_mono.markers, Bcell.markers, NK.markers,
#                                       Bcell_MemoryB.markers, Plasma.markers), add.cell.ids = c("Tregs",))


# fiter de gene markers  

top_pos_Hepatocyte.markers <- filter(Hepatocyte.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Hepatocyte.markers <- filter(Hepatocyte.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Tregs.markers <- filter(Tregs.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Tregs.markers <- filter(Tregs.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Th17_NK.markers <- filter(Th17_NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Th17_NK.markers <- filter(Th17_NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Tcell.markers <- filter(Tcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Tcell.markers <- filter(Tcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Tregs_Tcell.markers <- filter(Tregs_Tcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Tregs_Tcell.markers <- filter(Tregs_Tcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Th1.markers <- filter(Th1.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Th1.markers <- filter(Th1.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Macro_mono.markers <- filter(Macro_mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Macro_mono.markers <- filter(Macro_mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Germ.markers <- filter(Germ.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 2) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Germ.markers <- filter(Germ.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Myeloid_dend_mono.markers <- filter(Myeloid_dend_mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Myeloid_dend_mono.markers <- filter(Myeloid_dend_mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Astrocyte_NK.markers <- filter(Astrocyte_NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Astrocyte_NK.markers <- filter(Astrocyte_NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Mono_Macro.markers <- filter(Mono_Macro.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Mono_Macro.markers <- filter(Mono_Macro.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Bcell.markers <- filter(Bcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Bcell.markers <- filter(Bcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_NK.markers <- filter(NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_NK.markers <- filter(NK.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Bcell_MemoryB.markers <- filter(Bcell_MemoryB.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Bcell_MemoryB.markers <- filter(Bcell_MemoryB.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_Plasma.markers <- filter(Plasma.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 10, wt = abs(avg_log2FC))
top_neg_Plasma.markers <- filter(Plasma.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

# графики диф экспрессии

Idents(CSF.integrated) <- "celltype"

markers.to.plot <- c(row.names(top_neg_Plasma.markers))

# DotPlot(CSF.integrated, features = markers.to.plot, 
#         split.by = "stim", idents = "Hepatocyte") + RotatedAxis()

plots <- VlnPlot(CSF.integrated, features = markers.to.plot, split.by = "stim", idents = "Plasmacytoid dendritic cell",
                 pt.size = 0, combine = FALSE, split.plot = TRUE)
CombinePlots(plots = plots, ncol = 5)

# ###vision
# # devtools::install_github("YosefLab/VISION")
# library(VISION)
# signatures <- c("/mnt/Data10tb/student9_data/GSE138266_RAW/Signatures/c7.all.v7.2.symbols.gmt")
# 
# vision.obj_CSF <- Vision(CSF.integrated, signatures = signatures, projection_methods = "UMAP")
# 
# vision.obj_CSF <- analyze(vision.obj_CSF)
# 
# viewResults(vision.obj_CSF)
# 
# # vision.obj_csf_ms <- addProjection(vis, "UMAP", projection)


