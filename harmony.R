library(Matrix)
library(Seurat)
library(dplyr)
# install.packages("harmony")
library(harmony)


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
}

CSF <-  merge(CSF.list[[1]], y = CSF.list[2:12],
              add.cell.ids = c("CSF_MS_1", "CSF_MS_2", "CSF_MS_3", "CSF_MS_4", "CSF_MS_5", "CSF_MS_6.data",
                               "CSF_cont_1", "CSF_cont_2", "CSF_cont_3", "CSF_cont_4", "CSF_cont_5", "CSF_cont_6.data"))
# normalization
CSF <- SCTransform(CSF, verbose = TRUE)
CSF <- RunPCA(CSF, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
CSF <- CSF %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

CSF <- RunUMAP(CSF, reduction = "harmony", dims = 1:30)
CSF <- FindNeighbors(CSF, reduction = "harmony", dims = 1:30)
CSF <- FindClusters(CSF, resolution = 0.5)

p1 <- DimPlot(CSF.integrated, reduction = "umap", group.by = "stim")
p2 <- DimPlot(CSF.integrated, reduction = "umap", label = TRUE, repel = TRUE)

plot_name = paste('/data/ashevtsov/MS_data/','harmony_cluster_ind_integrated_umap.pdf', sep = '')
pdf(plot_name)
DimPlot(CSF, reduction = "umap", group.by = "seurat_clusters")
dev.off()
