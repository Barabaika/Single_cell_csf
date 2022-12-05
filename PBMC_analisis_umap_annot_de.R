#!/usr/bin/Rscript
# -*- coding: utf-8 -*-

if !require("Matrix") install.packages("Matrix")
if !require("Seurat") install.packages("Seurat")
if !require("dplyr") install.packages("dplyr")
if !require("harmony") install.packages("harmony")
if !require("SeuratDisk") install.packages("SeuratDisk")
if !require("ggplot2") install.packages("ggplot2")

# PARAMS

WORKDIR="/data/ashevtsov/MS_data/raw_data/"
setwd(WORKDIR)

min.cells = 3
min.features = 200

cell.ids <- c("PBMC_MS_1", "PBMC_MS_2", "PBMC_MS_3", "PBMC_MS_4", "PBMC_MS_5", 
              "PBMC_cont_1", "PBMC_cont_2", "PBMC_cont_3", "PBMC_cont_4", "PBMC_cont_5")

# LOAD DATA

PBMC_MS_1.data <- Read10X(file.path(WORKDIR, "GSM4104134_MS19270_PBMCs")
PBMC_MS_2.data <- Read10X(file.path(WORKDIR, "GSM4104135_MS71658_PBMCs")
PBMC_MS_3.data <- Read10X(file.path(WORKDIR, "GSM4104136_MS49131_PBMCs")
PBMC_MS_4.data <- Read10X(file.path(WORKDIR, "GSM4104137_MS60249_PBMCs")
PBMC_MS_5.data <- Read10X(file.path(WORKDIR, "GSM4104138_MS74594_PBMCs")

PBMC_cont_1.data <- Read10X(file.path(WORKDIR, "GSM4104139_PST83775_PBMCs")
PBMC_cont_2.data <- Read10X(file.path(WORKDIR, "GSM4104140_PTC32190_PBMCs")
PBMC_cont_3.data <- Read10X(file.path(WORKDIR, "GSM4104141_PST95809_PBMCs")
PBMC_cont_4.data <- Read10X(file.path(WORKDIR, "GSM4104142_PTC41540_PBMCs")
PBMC_cont_5.data <- Read10X(file.path(WORKDIR, "GSM4104143_PTC85037_PBMCs")

PBMC_MS_1.data <- CreateSeuratObject(counts = PBMC_MS_1.data, min.cells = min.cells, min.features = min.features, project = "PBMC_MS_1")
PBMC_MS_2.data <- CreateSeuratObject(counts = PBMC_MS_2.data, min.cells = min.cells, min.features = min.features, project = "PBMC_MS_2")
PBMC_MS_3.data <- CreateSeuratObject(counts = PBMC_MS_3.data, min.cells = min.cells, min.features = min.features, project = "PBMC_MS_3")
PBMC_MS_4.data <- CreateSeuratObject(counts = PBMC_MS_4.data, min.cells = min.cells, min.features = min.features, project = "PBMC_MS_4")
PBMC_MS_5.data <- CreateSeuratObject(counts = PBMC_MS_5.data, min.cells = min.cells, min.features = min.features, project = "PBMC_MS_5")

PBMC_cont_1.data <- CreateSeuratObject(counts = PBMC_cont_1.data, min.cells = min.cells, min.features = min.features, project = "PBMC_cont_1")
PBMC_cont_2.data <- CreateSeuratObject(counts = PBMC_cont_2.data, min.cells = min.cells, min.features = min.features, project = "PBMC_cont_2")
PBMC_cont_3.data <- CreateSeuratObject(counts = PBMC_cont_3.data, min.cells = min.cells, min.features = min.features, project = "PBMC_cont_3")
PBMC_cont_4.data <- CreateSeuratObject(counts = PBMC_cont_4.data, min.cells = min.cells, min.features = min.features, project = "PBMC_cont_4")
PBMC_cont_5.data <- CreateSeuratObject(counts = PBMC_cont_5.data, min.cells = min.cells, min.features = min.features, project = "PBMC_cont_5")

PBMC <- merge(PBMC_MS_1.data, c(PBMC_MS_2.data, PBMC_MS_3.data, PBMC_MS_4.data, PBMC_MS_5.data, PBMC_cont_1.data, PBMC_cont_2.data, PBMC_cont_3.data, PBMC_cont_4.data, PBMC_cont_5.data), 
              add.cell.ids = cell.ids)


PBMC.list <- SplitObject(PBMC)

PBMC.list$PBMC_MS_1$stim <- "MS"
PBMC.list$PBMC_MS_2$stim <- "MS"
PBMC.list$PBMC_MS_3$stim <- "MS"
PBMC.list$PBMC_MS_4$stim <- "MS"
PBMC.list$PBMC_MS_5$stim <- "MS"


PBMC.list$PBMC_cont_1$stim <- "CONT"
PBMC.list$PBMC_cont_2$stim <- "CONT"
PBMC.list$PBMC_cont_3$stim <- "CONT"
PBMC.list$PBMC_cont_4$stim <- "CONT"
PBMC.list$PBMC_cont_5$stim <- "CONT"


rm(PBMC_MS_1.data, PBMC_MS_2.data, PBMC_MS_3.data, PBMC_MS_4.data, 
   PBMC_MS_5.data, PBMC_cont_1.data, PBMC_cont_2.data, 
   PBMC_cont_3.data, PBMC_cont_4.data, PBMC_cont_5.data)
# run garbage collect to free up memory
gc()

for (i in 1:length(PBMC.list)) {
  # filtering
  PBMC.list[[i]][["percent.mt"]] <- PercentageFeatureSet(PBMC.list[[i]], pattern = "^MT-")
  PBMC.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(PBMC.list[[i]], pattern = "^RP[SL]")
  PBMC.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(PBMC.list[[i]], pattern = "^HB[^(P)]")
  
  PBMC.list[[i]] <- subset(PBMC.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  PBMC.list[[i]] <- SCTransform(PBMC.list[[i]], verbose = TRUE)
}


##### per batch for batch difference exploring

PBMC[["percent.mt"]] <- PercentageFeatureSet(PBMC, pattern = "^MT-")
PBMC[["percent.Ribo"]] <- PercentageFeatureSet(PBMC, pattern = "^RP[SL]")
PBMC[["percent.hemoglob"]] <- PercentageFeatureSet(PBMC, pattern = "^HB[^(P)]")
PBMC <- subset(PBMC, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
# normalization
PBMC <- SCTransform(PBMC, verbose = TRUE)
PBMC <- RunPCA(PBMC, verbose = FALSE)
PBMC <- RunUMAP(PBMC, dims = 1:30)
PBMC <- FindNeighbors(PBMC, reduction = "pca", dims = 1:30)
PBMC <- FindClusters(PBMC, resolution = 0.5)
PBMC <- SetIdent(PBMC, value = "orig.ident") # orig.ident
plot_name = paste('/data/ashevtsov/MS_data/PBMC/', 'ALL_Batches_merged_batches_split_umap.pdf', sep = '_')
pdf(plot_name)
DimPlot(PBMC, reduction = "umap", label = FALSE, repel = TRUE)
dev.off()
#####

options(future.globals.maxSize= 14000000000)

PBMC.features <- SelectIntegrationFeatures(object.list = PBMC.list, nfeatures = 3000)


PBMC.list <- PrepSCTIntegration(object.list = PBMC.list, anchor.features = PBMC.features, 
                                verbose = FALSE)

PBMC.anchors <- FindIntegrationAnchors(object.list = PBMC.list, normalization.method = "SCT", 
                                       anchor.features = PBMC.features, verbose = TRUE)
PBMC.integrated <- IntegrateData(anchorset = PBMC.anchors, normalization.method = "SCT", 
                                 verbose = FALSE) 


#clustering
PBMC.integrated <- RunPCA(PBMC.integrated, verbose = FALSE)
PBMC.integrated <- RunUMAP(PBMC.integrated, dims = 1:30)

PBMC.integrated <- FindNeighbors(PBMC.integrated, reduction = "pca", dims = 1:30)
PBMC.integrated <- FindClusters(PBMC.integrated, resolution = 0.5)

# SaveH5Seurat(PBMC.integrated,'/data/ashevtsov/MS_data/seurat_objects/PBMC_integrated_21_11_2022', overwrite = TRUE)
# PBMC.integrated <- LoadH5Seurat("/data/ashevtsov/MS_data/seurat_objects/PBMC_integrated_21_11_2022.h5seurat")

DefaultAssay(object = PBMC.integrated) <- "SCT"
PBMC.integrated <- SetIdent(PBMC.integrated, value = "seurat_clusters")
PBMC.markers <- FindAllMarkers(PBMC.integrated, min.pct = 0.3, logfc.threshold = 0.25)

plot_name = '/data/ashevtsov/MS_data/PBMC/PBMC_integrated_markers.csv'
write.csv(PBMC.markers, plot_name)

######## plots 
# p1 <- DimPlot(CSF.integrated, reduction = "umap", group.by = "stim")
# p2 <- DimPlot(CSF.integrated, reduction = "umap", label = TRUE, repel = TRUE)
# 
# plot_name = paste('/data/ashevtsov/MS_data/','batches_integrated_umap.pdf', sep = '_')
# pdf(plot_name)
# DimPlot(CSF.integrated, reduction = "umap", group.by = "orig.ident")
# dev.off()
# 
# DefaultAssay(object = CSF.integrated) <- "SCT"
# gene = 'IL7R'
# plot_name = paste('/data/ashevtsov/MS_data/', gene,'integrated_umap.pdf', sep = '_')
# pdf(plot_name)
# FeaturePlot(CSF.integrated, features = c(gene))
# dev.off()
#########




PBMC.integrated <- SetIdent(PBMC.integrated, value = "seurat_clusters")
cell_types_list <- c('na CD8+ T cell',
                     'monocyte' ,
                     'ab T cell',
                     'NK cell',
                     'a CD8+ T cell',
                     'na CD8+ T cell',
                     'ab T cell',
                     'naive B cell',
                     'naive B cell',
                     'a CD8+ T cell',
                     'ab T cell',
                     'monocyte',
                     'gd T cell',
                     'ab T cell',
                     'mDC',
                     'a CD8+ T cell',
                     'megakaryocyte',
                     'granulocyte',
                     'NK cell',
                     'megakaryocyte',
                     'GrB+ pDC',
                     'mDC2')

names(cell_types_list) <- levels(PBMC.integrated)
PBMC.integrated <- RenameIdents(PBMC.integrated, cell_types_list)
PBMC.integrated$CellType <- Idents(PBMC.integrated)
PBMC.integrated$celltype.stim <- paste(Idents(PBMC.integrated), PBMC.integrated$stim, sep = "_")

# SaveH5Seurat(PBMC.integrated,'/data/ashevtsov/MS_data/seurat_objects/PBMC_integrated_21_11_2022', overwrite = TRUE)

plot_name = paste('/data/ashevtsov/MS_data/PBMC/','PBMC_batches_integrated_clust_names_umap.pdf', sep = '')
pdf(plot_name, width=9, height=6)
DimPlot(PBMC.integrated, reduction = "umap", label = FALSE, split.by = 'stim') + 
  ggtitle('PBMC cell types') +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c(
    'na CD8+ T cell' = '#CC6633',
    'ab T cell' = '#990000',
    'gd T cell' = '#996666',
    'a CD8+ T cell' = '#FF6633',
    'NK cell' = '#CC3399',
    'granulocyte' = '#99FF33',
    'megakaryocyte' = '#006600',
    'naive B cell' = '#CC9900',
    'monocyte' = '#000099',
    'mDC' = '#0099FF',
    'mDC2' = '#66FFFF',
    'GrB+ pDC' = '#FF99FF'
  ))

dev.off()

plot_name = paste('/data/ashevtsov/MS_data/PBMC/','PBMC_batches_integrated_clust_names_umap.pdf', sep = '')
pdf(plot_name, width=9, height=6)
DimPlot(PBMC.integrated, reduction = "umap", label = TRUE)
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
PBMC.markers <- FindAllMarkers(PBMC.integrated, min.pct = 0.3, logfc.threshold = 0.25)

plot_name = '/data/ashevtsov/MS_data/PBMC/PBMC_integrated_markers.csv'
write.csv(PBMC.markers, plot_name)



## plot number of cells among condition
cells_counts <- data.frame(Idents(PBMC.integrated), PBMC.integrated$stim)
file_name = paste('/data/ashevtsov/MS_data/PBMC/','PBMC_cell_counts.csv', sep = '')
colnames(cells_counts) <- c('Cell_Type', 'Condition')
write_csv(cells_counts, file_name)

library(dplyr)
cells_counts <- cells_counts %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

library(ggplot2)
plot_name = paste('/data/ashevtsov/MS_data/PBMC/','PBMC_cell_percistange.pdf', sep = '')
pdf(plot_name, width=3, height=6)

ggplot(cells_counts, aes(x = Condition, y = percent, fill = Cell_Type))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(
    'na CD8+ T cell' = '#CC6633',
    'ab T cell' = '#990000',
    'gd T cell' = '#996666',
    'a CD8+ T cell' = '#FF6633',
    'NK cell' = '#CC3399',
    'granulocyte' = '#99FF33',
    'megakaryocyte' = '#006600',
    'naive B cell' = '#CC9900',
    'monocyte' = '#000099',
    'mDC' = '#0099FF',
    'mDC2' = '#66FFFF',
    'GrB+ pDC' = '#FF99FF'
  )) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
  ggtitle('PBMC') +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()