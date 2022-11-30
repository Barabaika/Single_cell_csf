remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
suppressMessages(require(DoubletFinder))

library(Matrix)
library(Seurat)
library(dplyr)


setwd ("/mnt/Data10tb/student9_data/GSE138266_RAW")
####loading file 
data.dir = "/mnt/Data10tb/student9_data/GSE138266_RAW/PBMC_cont_3"

csf_ms_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_1")
csf_ms_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_2")
csf_ms_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_3")
csf_ms_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_4")
csf_ms_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_5")
csf_ms_6.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_MS_6")
csf_cont_1.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_1")
csf_cont_2.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_2")
csf_cont_3.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_3")
csf_cont_4.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_4")
csf_cont_5.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_5")
csf_cont_6.data <- Read10X("/mnt/Data10tb/student9_data/GSE138266_RAW/CSF_cont_6")

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



########################## creating seurat object + QC FILTERING
csf_ms_1.data <- CreateSeuratObject(counts = csf_ms_1.data, min.cells = 3, min.features = 200, project = "csf_ms_1")
csf_ms_2.data <- CreateSeuratObject(counts = csf_ms_2.data, min.cells = 3, min.features = 200, project = "csf_ms_2")
csf_ms_3.data <- CreateSeuratObject(counts = csf_ms_3.data, min.cells = 3, min.features = 200, project = "csf_ms_3")
csf_ms_4.data <- CreateSeuratObject(counts = csf_ms_4.data, min.cells = 3, min.features = 200, project = "csf_ms_4")
csf_ms_5.data <- CreateSeuratObject(counts = csf_ms_5.data, min.cells = 3, min.features = 200, project = "csf_ms_5")
csf_ms_6.data <- CreateSeuratObject(counts = csf_ms_6.data, min.cells = 3, min.features = 200, project = "csf_ms_6")
csf_cont_1.data <- CreateSeuratObject(counts = csf_cont_1.data, min.cells = 3, min.features = 200, project = "csf_cont_1")
csf_cont_2.data <- CreateSeuratObject(counts = csf_cont_2.data, min.cells = 3, min.features = 200, project = "csf_cont_2")
csf_cont_3.data <- CreateSeuratObject(counts = csf_cont_3.data, min.cells = 3, min.features = 200, project = "csf_cont_3")
csf_cont_4.data <- CreateSeuratObject(counts = csf_cont_4.data, min.cells = 3, min.features = 200, project = "csf_cont_4")
csf_cont_5.data <- CreateSeuratObject(counts = csf_cont_5.data, min.cells = 3, min.features = 200, project = "csf_cont_5")
csf_cont_6.data <- CreateSeuratObject(counts = csf_cont_6.data, min.cells = 3, min.features = 200, project = "csf_cont_6")

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

#совместить все в одно
csf_ms <- merge(csf_ms_1.data, c(csf_ms_2.data, csf_ms_3.data, csf_ms_4.data, csf_ms_5.data,csf_ms_6.data), add.cell.ids = c("csf_ms_1", "csf_ms_2", "csf_ms_3", "csf_ms_4", "csf_ms_5", "csf_ms_6"))
csf_cont <- merge(csf_cont_1.data, c(csf_cont_2.data, csf_cont_3.data, csf_cont_4.data, csf_cont_5.data, csf_cont_6.data),add.cell.ids = c("csf_cont_1", "csf_cont_2", "csf_cont_3", "csf_cont_4", "csf_cont_5", "csf_cont_6" ))
pbmc_ms <- merge(pbmc_ms_1.data, c(pbmc_ms_2.data, pbmc_ms_3.data, pbmc_ms_4.data, pbmc_ms_5.data),
                 add.cell.ids = c("pbmc_ms_1", "pbmc_ms_2", "pbmc_ms_3", "pbmc_ms_4", "pbmc_ms_5"))
pbmc_cont <- merge(pbmc_cont_1.data, c(pbmc_cont_2.data, pbmc_cont_3.data, pbmc_cont_4.data, pbmc_cont_5.data), 
                   add.cell.ids = c("pbmc_cont_1", "pbmc_cont_2", "pbmc_cont_3", "pbmc_cont_4", "pbmc_cont_5" ))

# remove all objects that will not be used.
rm(csf_ms_1.data, csf_ms_2.data, csf_ms_3.data, csf_ms_4.data, 
   csf_ms_5.data,csf_ms_6.data, csf_cont_1.data, csf_cont_2.data, 
   csf_cont_3.data, csf_cont_4.data, csf_cont_5.data, csf_cont_6.data, 
   pbmc_ms_1.data, pbmc_ms_2.data, pbmc_ms_3.data, pbmc_ms_4.data, 
   pbmc_ms_5.data, pbmc_cont_1.data, pbmc_cont_2.data, 
   pbmc_cont_3.data, pbmc_cont_4.data, pbmc_cont_5.data)
# run garbage collect to free up memory
gc()

#  сделать список (обратное действие,но както по другому у меня не вышло=))
csf_ms.list <- SplitObject(csf_ms)
csf_cont.list <- SplitObject(csf_cont)
pbmc_ms.list <- SplitObject(pbmc_ms)
pbmc_cont.list <- SplitObject(pbmc_cont)

for (i in 1:length(csf_ms.list)) {
  # filtering
  csf_ms.list[[i]][["percent.mt"]] <- PercentageFeatureSet(csf_ms.list[[i]], pattern = "^MT-")
  csf_ms.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(csf_ms.list[[i]], pattern = "^RP[SL]")
  csf_ms.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(csf_ms.list[[i]], pattern = "^HB[^(P)]")
  
  csf_ms.list[[i]] <- subset(csf_ms.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  csf_ms.list[[i]] <- SCTransform(csf_ms.list[[i]], verbose = FALSE)
}

for (i in 1:length(csf_cont.list)) {
  # filtering
  csf_cont.list[[i]][["percent.mt"]] <- PercentageFeatureSet(csf_cont.list[[i]], pattern = "^MT-")
  csf_cont.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(csf_cont.list[[i]], pattern = "^RP[SL]")
  csf_cont.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(csf_cont.list[[i]], pattern = "^HB[^(P)]")
  
  csf_cont.list[[i]] <- subset(csf_cont.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  csf_cont.list[[i]] <- SCTransform(csf_cont.list[[i]], verbose = FALSE)
}

for (i in 1:length(pbmc_ms.list)) {
  # filtering
  pbmc_ms.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pbmc_ms.list[[i]], pattern = "^MT-")
  pbmc_ms.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(pbmc_ms.list[[i]], pattern = "^RP[SL]")
  pbmc_ms.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(pbmc_ms.list[[i]], pattern = "^HB[^(P)]")
  
  pbmc_ms.list[[i]] <- subset(pbmc_ms.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  pbmc_ms.list[[i]] <- SCTransform(pbmc_ms.list[[i]], verbose = FALSE)
}

for (i in 1:length(pbmc_cont.list)) {
  # filtering
  pbmc_cont.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pbmc_cont.list[[i]], pattern = "^MT-")
  pbmc_cont.list[[i]][["percent.Ribo"]] <- PercentageFeatureSet(pbmc_cont.list[[i]], pattern = "^RP[SL]")
  pbmc_cont.list[[i]][["percent.hemoglob"]] <- PercentageFeatureSet(pbmc_cont.list[[i]], pattern = "^HB[^(P)]")
  
  pbmc_cont.list[[i]] <- subset(pbmc_cont.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10 & percent.hemoglob <10)
  
  # normalization
  pbmc_cont.list[[i]] <- SCTransform(pbmc_cont.list[[i]], verbose = FALSE)
}


# change the max size of global objects, for defiting eror in csf_ms.list making (PrepSCTIntegration)
options(future.globals.maxSize= 14000000000)
# selecting features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
csf_ms.features <- SelectIntegrationFeatures(object.list = csf_ms.list, nfeatures = 3000)
csf_ms.list <- PrepSCTIntegration(object.list = csf_ms.list, anchor.features = csf_ms.features, 
                                    verbose = FALSE)
# identify anchors and integrate the datasets
csf_ms.anchors <- FindIntegrationAnchors(object.list = csf_ms.list, normalization.method = "SCT", 
                                           anchor.features = csf_ms.features, verbose = FALSE)
csf_ms.integrated <- IntegrateData(anchorset = csf_ms.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

#clastering
csf_ms.integrated <- RunPCA(csf_ms.integrated, verbose = FALSE)
csf_ms.integrated <- RunUMAP(csf_ms.integrated, dims = 1:30)

csf_ms.integrated <- FindNeighbors(csf_ms.integrated, reduction = "pca", dims = 1:30)
csf_ms.integrated <- FindClusters(csf_ms.integrated, resolution = 0.5)

p1 <- DimPlot(csf_ms.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(csf_ms.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# find markers for every cluster compared to all remaining cells, report only the positive ones
chf_ms.markers <- FindAllMarkers(csf_ms.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_chf_ms.markers <- chf.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(csf_ms.integrated, features = top_chf_ms.markers$gene)
###vision
###vision



library(VISION)
signatures <- c("/mnt/Data10tb/student9_data/GSE138266_RAW/Signatures/c7.all.v7.2.symbols.gmt")

vision.obj_pbmc_cont <- Vision(pbmc_cont.integrated, signatures = signatures, projection_methods = "UMAP")

vision.obj_pbmc_cont <- analyze(vision.obj_pbmc_cont)

viewResults(vision.obj_pbmc_cont)

vision.obj_csf_ms <- addProjection(vis, "UMAP", projection)


# finding cell types with scCATCH( only for pbmc becouse this there is no CSF tissue in scCATCH)
library(scCATCH)
pbmc_ms.markers_annot<- scCATCH(object = pbmc_ms.markers,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
colnames(pbmc_ms.markers_annot)[3] <- 'scCATCH_cell_type'

pbmc_cont.markers_annot<- scCATCH(object = pbmc_cont.markers,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
colnames(pbmc_cont.markers_annot)[3] <- 'scCATCH_cell_type'

#annotatig cell types to clasters 

pbmc_ms.cluster.ids <- pbmc_ms.markers_annot$scCATCH_cell_type
names(pbmc_ms.cluster.ids) <- levels(pbmc_ms.integrated)
pbmc_ms.integrated <- RenameIdents(pbmc_ms.integrated, pbmc_ms.cluster.ids)

DimPlot(pbmc_ms.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc_cont.cluster.ids <- pbmc_cont.markers_annot$scCATCH_cell_type
names(pbmc_cont.cluster.ids) <- levels(pbmc_cont.integrated)
pbmc_cont.integrated <- RenameIdents(pbmc_cont.integrated, pbmc_cont.cluster.ids)

DimPlot(pbmc_cont.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#####diff exp analysis
pbmc.merge <- merge(pbmc_ms.integrated, pbmc_cont.integrated, add.cell.ids = c("pbmc_ms.integrated", "pbmc_cont.integrated"))
pbmcall.list <- SplitObject(pbmc.merge)

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
monocyte.de.markers <- FindMarkers(pbmcall.list,  cells.1 = "pbmc_ms",cells.1 = "pbmc_cont")
# view results
head(monocyte.de.markers)


list_objects_batches = SplitObject(CSF.integrated, split.by <- 'orig.ident')
batch1 <- list_objects_batches[['CSF_MS_1']]
