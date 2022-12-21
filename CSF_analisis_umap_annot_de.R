library(Matrix)
library(Seurat)
library(dplyr)
# install.packages("harmony")
library(harmony)
library(SeuratDisk)



setwd ("/mnt/Data10tb/student9_data/GSE138266_RAW")
####loading file 
data.dir = "/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_3"

CSF_MS_1.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104122_MS19270_CSF")
CSF_MS_2.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104123_MS58637_CSF")
CSF_MS_3.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104124_MS71658_CSF")
CSF_MS_4.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104125_MS49131_CSF")
CSF_MS_5.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104126_MS60249_CSF")
CSF_MS_6.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104127_MS74594_CSF")

CSF_cont_1.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104128_PST83775_CSF")
CSF_cont_2.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104129_PTC32190_CSF")
CSF_cont_3.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104130_PST95809_CSF")
CSF_cont_4.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104131_PTC41540_CSF")
CSF_cont_5.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104132_PST45044_CSF")
CSF_cont_6.data <- Read10X("/data/ashevtsov/MS_data/raw_data/GSM4104133_PTC85037_CSF")

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


# CSF <- SetIdent(CSF, value = "orig.ident") # orig.ident
# plot_name = paste('/data/ashevtsov/MS_data/umap_batches/', 'ALL_Batches_merged_batches_split_umap.pdf', sep = '_')
# pdf(plot_name)
# DimPlot(CSF, reduction = "umap", label = FALSE, repel = TRUE)
# dev.off()
# CSF.list_stim <- SplitObject(CSF, split.by = "stim")
# 
# plot_name = paste('/data/ashevtsov/MS_data/umap_batches/', 'MS_Batches_merged_batches_split_umap.pdf', sep = '_')
# pdf(plot_name)
# DimPlot(CSF.list_stim$MS, reduction = "umap", label = FALSE, repel = TRUE)
# dev.off()
# 
# gene = 'CD8A'
# plot_name = paste('/data/ashevtsov/MS_data/umap_batches/', 'MS', gene,'Batches_merged_batches_split_umap.pdf', sep = '_')
# pdf(plot_name)
# FeaturePlot(CSF.list_stim$MS, features = c(gene))
# dev.off()


options(future.globals.maxSize= 14000000000)

CSF.features <- SelectIntegrationFeatures(object.list = CSF.list, nfeatures = 3000)


CSF.list <- PrepSCTIntegration(object.list = CSF.list, anchor.features = CSF.features, 
                                verbose = FALSE)

CSF.anchors <- FindIntegrationAnchors(object.list = CSF.list, normalization.method = "SCT", 
                                       anchor.features = CSF.features, verbose = TRUE)
CSF.integrated <- IntegrateData(anchorset = CSF.anchors, normalization.method = "SCT", 
                                 verbose = FALSE) 


#clustering
CSF.integrated <- RunPCA(CSF.integrated, verbose = FALSE)
CSF.integrated <- RunUMAP(CSF.integrated, dims = 1:30)

CSF.integrated <- FindNeighbors(CSF.integrated, reduction = "pca", dims = 1:30)
CSF.integrated <- FindClusters(CSF.integrated, resolution = 0.5)

######## plots 
plot_name = paste('/data/ashevtsov/MS_data/','integrated_umap_mas_and_cont.pdf', sep = '_')
pdf(plot_name)
DimPlot(CSF.integrated, reduction = "umap", group.by = c("stim", "orig.ident"))
dev.off()


CSF.integrated <- SetIdent(CSF.integrated, value = "seurat_clusters")
cell_types_list <- c('T cell',
                     'ab T cell',
                     'a CD8+ T cell',
                     'a CD8+ T cell',
                     'na CD8+ T cell',
                     'T cell',
                     'ab T cell',
                     'ab T cell',
                     'monocyte',
                     'T cell',
                     'mDC2',
                     'mDC',
                     'NK cell',
                     'gd T cell',
                     'naive B cell',
                     'granulocyte',
                     'Plasmablast',
                     'GrB+ pDC',
                     'mDC1')

names(cell_types_list) <- levels(CSF.integrated)
CSF.integrated <- RenameIdents(CSF.integrated, cell_types_list)
CSF.integrated$CellType <- Idents(CSF.integrated)
CSF.integrated$celltype.stim <- paste(Idents(CSF.integrated), CSF.integrated$stim, sep = "_")

# SaveH5Seurat(CSF.integrated,'/data/ashevtsov/MS_data/seurat_objects/CSF_integrated_21_11_2022', overwrite = TRUE)
# CSF.integrated <- LoadH5Seurat("/data/ashevtsov/MS_data/seurat_objects/CSF_integrated_21_11_2022.h5seurat")

plot_name = paste('/data/ashevtsov/MS_data/CSF/','CSF_batches_integrated_clust_names_umap.pdf', sep = '')
pdf(plot_name, width=9, height=6)
DimPlot(CSF.integrated, reduction = "umap", label = FALSE, split.by = 'stim') + 
  ggtitle('CSF cell types') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c(
  'na CD8+ T cell' = '#CC6633',
  'a CD8+ T cell' = '#FF6633',
  'ab T cell' = '#990000',
  'T cell' = '#330000',
  'gd T cell' = '#996666',
  'CD4+ T cell' = '#FF9999',
  'NK cell' = '#CC3399',
  'granulocyte' = '#99FF33',
  'Plasmablast' = '#006600',
  'naive B cell' = '#CC9900',
  'monocyte' = '#000099',
  'mDC' = '#0099FF',
  'mDC2' = '#66FFFF',
  'mDC1' = '#33FFFF',
  'GrB+ pDC' = '#FF99FF'
  ))

dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
CSF.markers <- FindAllMarkers(CSF.integrated, min.pct = 0.3, logfc.threshold = 0.25)

plot_name = '/data/ashevtsov/MS_data/integrated_markers.csv'
write.csv(CSF.markers, plot_name)

file_name = paste('/data/ashevtsov/MS_data/CSF/','CSF_cell_counts.csv', sep = '')
cells_counts <- as.data.frame.matrix(table(Idents(CSF.integrated), CSF.integrated$stim))
write.csv(cells_counts, file_name)

## plot number of cells among condition
cells_counts <- data.frame(Idents(CSF.integrated), CSF.integrated$stim)
file_name = paste('/data/ashevtsov/MS_data/CSF/','CSF_cell_counts.csv', sep = '')
colnames(cells_counts) <- c('Cell_Type', 'Condition')
write_csv(cells_counts, file_name)

library(dplyr)
cells_counts <- cells_counts %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

library(ggplot2)
plot_name = paste('/data/ashevtsov/MS_data/CSF/','CSF_cell_percistange.pdf', sep = '')
pdf(plot_name, width=3, height=6)

ggplot(cells_counts, aes(x = Condition, y = percent, fill = Cell_Type))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(
    'na CD8+ T cell' = '#CC6633',
    'a CD8+ T cell' = '#FF6633',
    'ab T cell' = '#990000',
    'T cell' = '#330000',
    'gd T cell' = '#996666',
    'CD4+ T cell' = '#FF9999',
    'NK cell' = '#CC3399',
    'granulocyte' = '#99FF33',
    'Plasmablast' = '#006600',
    'naive B cell' = '#CC9900',
    'monocyte' = '#000099',
    'mDC' = '#0099FF',
    'mDC2' = '#66FFFF',
    'mDC1' = '#33FFFF',
    'GrB+ pDC' = '#FF99FF'
  )) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())+
  ggtitle('CSF') +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()