set.seed(12345)
library(Seurat)
library(hdf5r)
library(patchwork)
library(harmony)
library(cowplot)
library(dplyr)
library(limma)
library(reticulate)
use_condaenv(condaenv = "py39", conda = "/home/anaconda3/bin/conda")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

bm1_1<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp1/Run4xSI-GA-H11_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM1_1", min.cells = 3, min.features = 200)
bm1_1[["percent.mt"]] <- PercentageFeatureSet(bm1_1, pattern = "^MT-")
plot1 <- FeatureScatter(bm1_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm1_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm1_1 <- subset(bm1_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#LayerData(bm1_1,layer = 'data')<- LayerData(bm1_1,layer = 'counts')
bm1_1 <- NormalizeData(bm1_1)
bm1_1 <- CellCycleScoring(bm1_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)


#bm1_1 <- SCTransform(bm1_1, method = "glmGamPoi",vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#------------
bm1_2<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp1/Run5xSI-GA-D10_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM1_2", min.cells = 3, min.features = 200)
bm1_2[["percent.mt"]] <- PercentageFeatureSet(bm1_2, pattern = "^MT-")
plot1 <- FeatureScatter(bm1_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm1_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm1_2 <- subset(bm1_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10)
bm1_2 <- NormalizeData(bm1_2)
bm1_2 <- CellCycleScoring(bm1_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm1_2 <- SCTransform(bm1_2, method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#---------------
bm2_1<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp2/BM1_IGO_07465_1_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM2_1", min.cells = 3, min.features = 200)
bm2_1[["percent.mt"]] <- PercentageFeatureSet(bm2_1, pattern = "^MT-")
plot1 <- FeatureScatter(bm2_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm2_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm2_1 <- subset(bm2_1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
bm2_1 <- NormalizeData(bm2_1)
bm2_1 <- CellCycleScoring(bm2_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm2_1 <- SCTransform(bm2_1, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#-----------------
bm2_2<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp2/BM1_IGO_07465_1_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM2_2", min.cells = 3, min.features = 200)
bm2_2[["percent.mt"]] <- PercentageFeatureSet(bm2_2, pattern = "^MT-")
plot1 <- FeatureScatter(bm2_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm2_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm2_2 <- subset(bm2_2, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
bm2_2 <- NormalizeData(bm2_2)
bm2_2 <- CellCycleScoring(bm2_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm2_2 <- SCTransform(bm2_2, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#------------------
bm2_3<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp2/BM2_10xSI-GA-C12_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM2_3", min.cells = 3, min.features = 200)
bm2_3[["percent.mt"]] <- PercentageFeatureSet(bm2_3, pattern = "^MT-")
plot1 <- FeatureScatter(bm2_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm2_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm2_3 <- subset(bm2_3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
bm2_3 <- NormalizeData(bm2_3)
bm2_3 <- CellCycleScoring(bm2_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm2_3 <- SCTransform(bm2_3, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#-----------------
bm2_4<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp2/BM2_10xSI-GA-D12_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM2_4", min.cells = 3, min.features = 200)
bm2_4[["percent.mt"]] <- PercentageFeatureSet(bm2_4, pattern = "^MT-")
plot1 <- FeatureScatter(bm2_4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm2_4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm2_4 <- subset(bm2_4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
bm2_4 <- NormalizeData(bm2_4)
bm2_4 <- CellCycleScoring(bm2_4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm2_4 <- SCTransform(bm2_4, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#------------------
bm2_5<- Read10X_h5( "/scratch/CD34cells_Palantir/bmp2/BM2_10xSI-GA-F12_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM2_5", min.cells = 3, min.features = 200)
bm2_5[["percent.mt"]] <- PercentageFeatureSet(bm2_5, pattern = "^MT-")
plot1 <- FeatureScatter(bm2_5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm2_5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm2_5 <- subset(bm2_5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
bm2_5 <- NormalizeData(bm2_5)
bm2_5 <- CellCycleScoring(bm2_5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm2_5 <- SCTransform(bm2_5, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#-----------------
bm3_1<- Read10X_h5("/scratch/CD34cells_Palantir/bmp3/BoneMarrow_CD34_1_IGO_07861_1_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM3_1", min.cells = 3, min.features = 200)
bm3_1[["percent.mt"]] <- PercentageFeatureSet(bm3_1, pattern = "^MT-")
plot1 <- FeatureScatter(bm3_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm3_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm3_1 <- subset(bm3_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
bm3_1 <- NormalizeData(bm3_1)
bm3_1 <- CellCycleScoring(bm3_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm3_1 <- SCTransform(bm3_1, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#--------------------
bm3_2<- Read10X_h5("/scratch/CD34cells_Palantir/bmp3/BoneMarrow_CD34_2_IGO_07861_2_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM3_2", min.cells = 3, min.features = 200)
bm3_2[["percent.mt"]] <- PercentageFeatureSet(bm3_2, pattern = "^MT-")
plot1 <- FeatureScatter(bm3_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm3_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm3_2 <- subset(bm3_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10)
bm3_2 <- NormalizeData(bm3_2)
bm3_2 <- CellCycleScoring(bm3_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm3_2 <- SCTransform(bm3_2, method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#--------------------
bm3_3<- Read10X_h5("/scratch/CD34cells_Palantir/bmp3/BoneMarrow_CD34_3_IGO_07861_3_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM3_3", min.cells = 3, min.features = 200)
bm3_3[["percent.mt"]] <- PercentageFeatureSet(bm3_3, pattern = "^MT-")
plot1 <- FeatureScatter(bm3_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm3_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm3_3 <- subset(bm3_3, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10)
bm3_3 <- NormalizeData(bm3_3)
bm3_3 <- CellCycleScoring(bm3_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm3_3 <- SCTransform(bm3_3, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)

#------------------
bm3_4<- Read10X_h5("/scratch/CD34cells_Palantir/bmp3/BoneMarrow_CD34_4_IGO_07861_4_out/outs/filtered_feature_bc_matrix.h5") %>%
  CreateSeuratObject( project = "BM3_4", min.cells = 3, min.features = 200)
bm3_4[["percent.mt"]] <- PercentageFeatureSet(bm3_4, pattern = "^MT-")
plot1 <- FeatureScatter(bm3_4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm3_4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm3_4 <- subset(bm3_4, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10)
bm3_4 <- NormalizeData(bm3_4)
bm3_4 <- CellCycleScoring(bm3_4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#bm3_4 <- SCTransform(bm3_4, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)


#-------------------
#intFeater<- SelectIntegrationFeatures(object.list =c(bm1_1,bm1_2,bm2_1,bm2_2,bm2_3,bm2_4,bm2_5,bm3_1,bm3_2,bm3_3,bm3_4), nfeatures = 3000, verbose = FALSE )
#comib_prep<- PrepSCTIntegration(object.list =c(bm1_1,bm1_2,bm2_1,bm2_2,bm2_3,bm2_4,bm2_5,bm3_1,bm3_2,bm3_3,bm3_4), assay="SCT", anchor.features = intFeater)
#intAnchor<- FindIntegrationAnchors(comib_prep,normalization.method = "SCT", anchor.features = intFeater, reduction = "cca", verbose = FALSE)#object.list =c(Soluplus,SR_1,UM171))
#hspc_comb <- IntegrateData(intAnchor, normalization.method = "SCT",new.assay.name = "Combine", verbose = FALSE)
hspc_comb <- merge(x= bm1_1,y=list(bm1_2,bm2_1,bm2_2,bm2_3,bm2_4,bm2_5,bm3_1,bm3_2,bm3_3,bm3_4))
bm1_1=c()
bm1_2 = c()
bm2_1 = c()
bm2_2 = c()
bm2_3 = c()
bm2_4 = c()
bm2_5 = c()
bm3_1 = c()
bm3_2 = c()
bm3_3 = c()
bm3_4 = c()
#hspc_comd <- NormalizeData(hspc_comd)
#hspc_comd <-

hspc_comb<- SCTransform(hspc_comb, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE)
DefaultAssay(hspc_comb)<- 'SCT'
hspc_comb <- RunPCA(hspc_comb, npcs=50, assay = 'SCT')
hspc_comb <- IntegrateLayers(
  object = hspc_comb, method = CCAIntegration,
  layers = c('data.SoluP', 'data.SR1', 'data.UM171'),
  normalization.method = 'SCT',
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
#hspc_comb<- JoinLayers(hspc_comb,layers = c('counts.SoluP, counts.SR1, counts.UM171))
saveRDS(hspc_comb, file = paste0("IntGrat_hspc_ManuReproduce_combined_CellCycle_",Sys.Date(),".rds"))
#quit(save = "no",status = 0)
#hspc_comb <- CellCycleScoring(hspc_comb, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, verbose = FALSE)

#hspc_comb <- ScaleData(hspc_comb, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(hspc_comb), verbose = FALSE)

#hspc_comb <- ScaleData(hspc_comb, features = rownames(hspc_comb), verbose = FALSE)
#hspc_comb <- RunPCA(hspc_comb,npcs = 40, verbose = FALSE)

hspc_comb<- FindNeighbors(hspc_comb,  dims = 1:50, verbose = FALSE,reduction = "integrated.cca")#, graph.name = "Combine_snn")
hspc_comb<- FindClusters(hspc_comb, resolution = 0.45, algorithm = 4,method = "igraph")#, graph.name = "Combine_snn")
hspc_comb<- RunUMAP(hspc_comb,  dims = 1:50, verbose = FALSE, reduction = 'integrated.cca', metric = "euclidean", return.model = TRUE)
DimPlot(hspc_comb, label = T)
#hspc_comb<- JoinLayers(hspc_comb,assay = 'SCT')
hspc_comb<- PrepSCTFindMarkers(object = hspc_comb)
hspc_comb.markers<- FindAllMarkers(hspc_comb, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

hspc_comb.markers %>% group_by(cluster) %>% top_n(2,avg_log2FC) -> top2
saveRDS(hspc_comb, file = paste0("hspc_ManuReproduce_combined_CellCycle_",Sys.Date(),".rds"))
saveRDS(hspc_comb.markers,  file = paste0("hspc_ManuReporduce_SCT_allmarkers_CellCycle_",Sys.Date(),".rds"))
#hspc_comb<- JoinLayers(hspc_comb)
DefaultAssay(hspc_comb)<- "RNA"
hspc_comb<- JoinLayers(hspc_comb, assay = 'RNA')
hspc_comb.RNAmarkers<- FindAllMarkers(hspc_comb, assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
hspc_comb.RNAmarkers %>% group_by(cluster) %>% top_n(2,avg_log2FC) -> rnatop2
saveRDS(hspc_comb.RNAmarkers,  file = paste0("hspc_ManuReporduce_RNA_allmarkers_CellCycle_",Sys.Date(),".rds"))
