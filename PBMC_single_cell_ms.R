# devtools::install_github('ZJUFanLab/scCATCH')
library(Matrix)
library(Seurat)
library(dplyr)
library(cowplot)

setwd ("/mnt/Data10tb/student9_data/GSE138266_RAW")
####loading file 
data.dir = "/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_3"

pbmc_ms_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_MS_1")
pbmc_ms_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_MS_2")
pbmc_ms_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_MS_3")
pbmc_ms_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_MS_4")
pbmc_ms_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_MS_5")
pbmc_cont_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_1")
pbmc_cont_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_2")
pbmc_cont_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_3")
pbmc_cont_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_4")
pbmc_cont_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_5")

pbmc_ms_1.data <- CreateSeuratObject(counts = pbmc_ms_1.data, min.cells = 3, min.features = 200, project = "pbmc_ms_1")
pbmc_ms_2.data <- CreateSeuratObject(counts = pbmc_ms_2.data, min.cells = 3, min.features = 200, project = "pbmc_ms_2")
pbmc_ms_3.data <- CreateSeuratObject(counts = pbmc_ms_3.data, min.cells = 3, min.features = 200, project = "pbmc_ms_3")
pbmc_ms_4.data <- CreateSeuratObject(counts = pbmc_ms_4.data, min.cells = 3, min.features = 200, project = "pbmc_ms_4")
pbmc_ms_5.data <- CreateSeuratObject(counts = pbmc_ms_5.data, min.cells = 3, min.features = 200, project = "pbmc_ms_5")
pbmc_cont_1.data <- CreateSeuratObject(counts = pbmc_cont_1.data, min.cells = 3, min.features = 200, project = "pbmc_cont_1")
pbmc_cont_2.data <- CreateSeuratObject(counts = pbmc_cont_2.data, min.cells = 3, min.features = 200, project = "pbmc_cont_2")
pbmc_cont_3.data <- CreateSeuratObject(counts = pbmc_cont_3.data, min.cells = 3, min.features = 200, project = "pbmc_cont_3")
pbmc_cont_4.data <- CreateSeuratObject(counts = pbmc_cont_4.data, min.cells = 3, min.features = 200, project = "pbmc_cont_4")
pbmc_cont_5.data <- CreateSeuratObject(counts = pbmc_cont_5.data, min.cells = 3, min.features = 200, project = "pbmc_cont_5")

pbmc <- merge(pbmc_ms_1.data, c(pbmc_ms_2.data, pbmc_ms_3.data, pbmc_ms_4.data, pbmc_ms_5.data, pbmc_cont_1.data, pbmc_cont_2.data, pbmc_cont_3.data, pbmc_cont_4.data, pbmc_cont_5.data), 
              add.cell.ids = c("pbmc_ms_1", "pbmc_ms_2", "pbmc_ms_3", "pbmc_ms_4", "pbmc_ms_5",
                               "pbmc_cont_1", "pbmc_cont_2", "pbmc_cont_3", "pbmc_cont_4", "pbmc_cont_5"))

rm(csf_ms_1.data, csf_ms_2.data, csf_ms_3.data, csf_ms_4.data, 
   csf_ms_5.data,csf_ms_6.data, csf_cont_1.data, csf_cont_2.data, 
   csf_cont_3.data, csf_cont_4.data, csf_cont_5.data, csf_cont_6.data, 
   pbmc_ms_1.data, pbmc_ms_2.data, pbmc_ms_3.data, pbmc_ms_4.data, 
   pbmc_ms_5.data, pbmc_cont_1.data, pbmc_cont_2.data, 
   pbmc_cont_3.data, pbmc_cont_4.data, pbmc_cont_5.data)
# run garbage collect to free up memory
gc()

pbmc.list <- SplitObject(pbmc)
rm(pbmc)

pbmc.list$pbmc_ms_1$stim <- "MS"
pbmc.list$pbmc_ms_2$stim <- "MS"
pbmc.list$pbmc_ms_3$stim <- "MS"
pbmc.list$pbmc_ms_4$stim <- "MS"
pbmc.list$pbmc_ms_5$stim <- "MS"
pbmc.list$pbmc_cont_1$stim <- "CONT"
pbmc.list$pbmc_cont_2$stim <- "CONT"
pbmc.list$pbmc_cont_3$stim <- "CONT"
pbmc.list$pbmc_cont_4$stim <- "CONT"
pbmc.list$pbmc_cont_5$stim <- "CONT"

for (i in 1:length(pbmc.list)) {
  # filtering
  pbmc.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pbmc.list[[i]], pattern = "^MT-")
  pbmc.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(pbmc.list[[i]], pattern = "^RP[SL]")
  pbmc.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(pbmc.list[[i]], pattern = "^HB[^(P)]")
  
  pbmc.list[[i]] <- subset(pbmc.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10 & percent.Ribo < 10)
  
  # normalization
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
}


options(future.globals.maxSize= 14000000000)

pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)

pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features, 
                                   verbose = FALSE)

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                          anchor.features = pbmc.features, verbose = FALSE)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", 
                                    verbose = FALSE)


#clastering
pbmc.integrated <- RunPCA(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30)

pbmc.integrated <- FindNeighbors(pbmc.integrated, reduction = "pca", dims = 1:30)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.5)

p1 <- DimPlot(pbmc.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(pbmc.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_pbmc.markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(pbmc.integrated, features = top_chf_ms.markers$gene)

###vision

library(VISION)
signatures <- c("/mnt/Data10tb/student9_data/GSE138266_RAW/Signatures/c7.all.v7.2.symbols.gmt")

vision.obj_pbmc <- Vision(pbmc.integrated, signatures = signatures, projection_methods = "UMAP")

vision.obj_pbmc <- analyze(vision.obj_pbmc)

viewResults(vision.obj_pbmc)

# vision.obj_csf_ms <- addProjection(vis, "UMAP", projection)


# finding cell types with scCATCH( only for pbmc becouse this there is no CSF tissue in scCATCH)
library(scCATCH)
pbmc.markers_annot<- scCATCH(object = pbmc.markers,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
colnames(pbmc.markers_annot)[3] <- 'scCATCH_cell_type'

# pbmc_cont.markers_annot<- scCATCH(object = pbmc_cont.markers,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
# colnames(pbmc_cont.markers_annot)[3] <- 'scCATCH_cell_type'

#annotatig cell types to clasters 

pbmc.cluster.ids <- pbmc.markers_annot$scCATCH_cell_type
names(pbmc.cluster.ids) <- levels(pbmc.integrated)
pbmc.integrated <- RenameIdents(pbmc.integrated, pbmc.cluster.ids)

DimPlot(pbmc.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DimPlot(pbmc.integrated, reduction = "umap", split.by = "stim")

#####diff exp analysis
DefaultAssay(pbmc.integrated) <- "RNA"
# 
# naive.t.cells <- subset(pbmc.integrated, idents = "Naive T Cell")
# Idents(naive.t.cells) <- "stim"
# avg.t.cells <- AverageExpression(naive.t.cells, verbose = FALSE)$RNA
# avg.t.cells$gene <- rownames(avg.t.cells)
# library(ggplot2)
# p3 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")

###########
# Find differentially expressed features between claster and all
# BiocManager::install('multtest')
# install.packages('metap')
# library(metap)

VSall_naiveTcell.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Naive T Cell", grouping.var = "stim", verbose = FALSE)
# head(naiveTcell.markers)
VSall_mono.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Monocyte", grouping.var = "stim", verbose = FALSE)

VSall_memoryTcell.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Memory T Cell", grouping.var = "stim", verbose = FALSE)

VSall_regulatoryTcell.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Regulatory T Cell", grouping.var = "stim", verbose = FALSE)

VSall_transitionalBcell.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Transitional B Cell", grouping.var = "stim", verbose = FALSE)

VSall_plasmacytoiddendriticcell.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Plasmacytoid Dendritic Cell", grouping.var = "stim", verbose = FALSE)

VSall_basophil.markers <- FindConservedMarkers(pbmc.integrated, ident.1 = "Basophil",  grouping.var = "stim", verbose = FALSE)





# dif exp between each claster by condition

pbmc.integrated$celltype.stim <- paste(Idents(pbmc.integrated), pbmc.integrated$stim, sep = "_")
pbmc.integrated$celltype <- Idents(pbmc.integrated)
Idents(pbmc.integrated) <- "celltype.stim"
# 1 cluster
naiveTcell.markers <- FindMarkers(pbmc.integrated, ident.1 = "Naive T Cell_MS", ident.2 = "Naive T Cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 2 cluster
mono.markers <- FindMarkers(pbmc.integrated, ident.1 = "Monocyte_MS", ident.2 = "Monocyte_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 3 cluster
memoryTcell.markers <- FindMarkers(pbmc.integrated, ident.1 = "Memory T Cell_MS", ident.2 = "Memory T Cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 4 
regulatoryTcell.markers <- FindMarkers(pbmc.integrated, ident.1 = "Regulatory T Cell_MS", ident.2 = "Regulatory T Cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 5 cluster
transitionalBcell.markers <- FindMarkers(pbmc.integrated, ident.1 = "Transitional B Cell_MS", ident.2 = "Transitional B Cell_CONT",min.pct = 0.4, logfc.threshold = 0.25, verbose = FALSE)
# 6 cluster
plasmacytoiddendriticcell.markers <- FindMarkers(pbmc.integrated, ident.1 = "Plasmacytoid Dendritic Cell_MS",  ident.2 = "Plasmacytoid Dendritic Cell_CONT",min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
# 7 cluster
basophil.markers <- FindMarkers(pbmc.integrated, ident.1 = "Basophil_MS", ident.2 = "Basophil_CONT",min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)



# fiter de gene markers  
# filter(mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, abs(avg_log2FC) > 3)
top_pos_mono.markers <- filter(mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 2) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_mono.markers <- filter(mono.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_memoryT.markers <- filter(memoryTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_memoryT.markers <- filter(memoryTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_naiveTcell.markers <- filter(naiveTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_naiveTcell.markers <- filter(naiveTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_plasm.markers <- filter(plasmacytoiddendriticcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_plasm.markers <- filter(plasmacytoiddendriticcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_reg.markers <- filter(regulatoryTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_reg.markers <- filter(regulatoryTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))

top_pos_transB.markers <- filter(transitionalBcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 0) %>% top_n(n = 30, wt = abs(avg_log2FC))
top_neg_transB.markers <- filter(transitionalBcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC < 0) %>% top_n(n = 10, wt = abs(avg_log2FC))


# plot variable features
markers.to.plot <- c(row.names(top_pos_plasm.markers), row.names(top_neg_plasm.markers))
# plots <- VlnPlot(pbmc.integrated, features = markers.to.plot, split.by = "stim", group.by = "celltype", 
#                  pt.size = 0, combine = FALSE, split.plot = TRUE)
# CombinePlots(plots = plots, ncol = 1)

DotPlot(pbmc.integrated, features = markers.to.plot, 
        split.by = "stim", idents = "Plasmacytoid Dendritic Cell") + RotatedAxis()

# DotPlot(pbmc.integrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, 
# ) + RotatedAxis()


# library(presto)
# 
# markers<- presto::wilcoxauc(pbmc.integrated, assay = 'data', group_by	= "celltype.stim")
# 
# markers<- filter(markers, padj <= 0.01, pct_in >= 0.4, pct_out >= 0.4, logFC > 0)
# 
# all_markers<- markers %>%
#   select(-rank) %>% 
#   unclass() %>% 
#   stack() %>%
#   pull(values) %>%
#   unique() %>%
#   .[!is.na(.)]
# DoHeatmap(pbmc, features = all_markers) + NoLegend()
# 
# mat<- pbmc.integrated[["RNA"]]@data[all_markers, ] %>% as.matrix()


# GSEA analisys##############
pos_memoryT.markers <- filter(memoryTcell.markers, p_val_adj <= 0.01, pct.1 >= 0.4, pct.2 >= 0.4, avg_log2FC > 1)

ranks <- pos_memoryT.markers$avg_log2FC
namrangs <- row.names(pos_memoryT.markers)

write.csv(namrangs, row.names = FALSE)

library(msigdbr)
m_df<- msigdbr(species = "Homo sapiens", category = "C7")

head(m_df)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

library(fgsea)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, scoreType = "pos")

# only plot the top 20 pathways
ggplot(fgseaRes %>% filter(padj < 0.01) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#GO analisis
# BiocManager::install("geneLenDataBase")

library(goseq)
library(geneLenDataBase)
supportedOrganisms() %>% filter(str_detect(Genome, "hsa"))

library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code')

sigGenes <- shrinkLvV$Entrez[ shrinkLvV$FDR < 0.01 & 
                                !is.na(shrinkLvV$FDR) &
                                abs(shrinkLvV$logFC) > 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'hsa')
head(kk, n=10)