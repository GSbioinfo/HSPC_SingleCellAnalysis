set.seed(12345)
library(Seurat)
library(hdf5r)
library(patchwork)
#library(harmony)
library(cowplot)
library(dplyr)
library(limma)
library(reticulate)
use_virtualenv('/root/.virtualenvs/r-reticulate')

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

bmcd34_markers<- readRDS("/scratch/CD34cells_Palantir/PalantirCD34_RNA_allmarkers_CellCycle_2025-01-06.rds")
bmcd34_markers %>% group_by(cluster) %>% top_n(5,avg_log2FC) -> top5
bmcd34_markers %>% group_by(cluster) %>% top_n(2,avg_log2FC) -> top2
bmcd34_markers %>% group_by(cluster) %>% top_n(10,avg_log2FC) -> top10 

bmcd34_combined<- readRDS("/scratch/CD34cells_Palantir/PalantirCD34_combined_CellCycle_2025-01-06.rds")
# cluste -> cellType
#  
# 1 -> HSC/MPP
# 2 -> eryth
# 3 -> ProB
# 4 -> GMP
# 5 -> LT-HSC
# 6 -> ProMyelo
# 7 -> MegkP
# 8 -> LMPP
# 9 -> pDC
# 10 -> CLP
# 11 -> PreProB
# 12 -> DCp
# 13 -> NMP (Neutrophil-myeloid progenitor)
# 14 -> preB
# 15 -> Mast cell

bmcd34_combined <- RenameIdents(bmcd34_combined, # genes
                                `1` = 'HSC/MPP', # SPINK2
                                `2` = 'eryth',   # HBB
                                `3` = 'ProB',    # ARPP21, VPREB3, CD79B
                                `4`= 'GMP',      # ELANE, MPO
                                `5` = 'LT-HSC',  # AVP, HLF
                                `6` = 'ProMyelo',# LYZ, MPO
                                `7` = 'MegkP',   # PLXDC2, GP1BB
                                `8` = 'LMPP',    # NEGR1, CNTNAP2, CD99
                                `9` = 'pDC',     # IRF8, TCF4
                                `10` = 'CLP',    # EBF1, DNTT, LTB
                                `11` = 'PreProB',# HIST1H4C, TUBA1B
                                `12` = 'DCp',    # LYZ, CST3, S100A10, S100A9
                                `13` = 'NMP',    # CLC, MS4A3
                                `14` = 'preB',   # CD24, IGHM
                                `15` = 'MastCell')# FTL, FTH1, C1QA

ggsave2("Annotated_dimplot.pdf", DimPlot(bmcd34_combined,label = T), height = 10, width = 10)
ggsave2("dimplot_split_by_donor_cellcycle.pdf", DimPlot(bmcd34_combined,group.by = "orig.ident", split.by = "Phase"), width = 12, height = 4)
ggsave2("dimplot_cellcycle.pdf",DimPlot(bmcd34_combined,group.by = "Phase"), width = 12, height = 4)
#bmcd34_markers<- readRDS("PalantirCD34_SCT_allmarkers_CellCycle_2023-07-10.rds")
#bmcd34_markers %>% group_by(cluster) %>% top_n(5,avg_log2FC) -> top5
#bmcd34_markers %>% group_by(cluster) %>% top_n(2,avg_log2FC) -> top2
#bmcd34_markers %>% group_by(cluster) %>% top_n(10,avg_log2FC) -> top10
#DotPlot(bmcd34_combined,assay = "SCT",features = unique(top2$gene))+RotatedAxis()

saveRDS(bmcd34_combined, file = paste0("/scratch/CD34cells_Palantir/Anno_PalantirCD34_combined_CellCycle_",Sys.Date(),".rds"))
stem_gene<- c("HLF", "AVP", "CD34","PROM1","TOP2A","MKI67","S100A8","S100A9",
              "MPO","PRTN3", "ELANE", "LTB", "JCHAIN", "GATA2",
              "ITGA2B", "KLF1", "TFRC", "HBG2", "AHSP", "PF4", 
              "ITGB3", "KIT", "KRT1")

DefaultAssay(bmcd34_combined)<- "Combine"
dtplt<- DotPlot(bmcd34_combined,assay = "SCT",features = stem_gene)+RotatedAxis()
ggsave2("Dotplot_select_stem_genes.pdf",dtplt, height = 12, width = 10) 
#bmcd34_combined<- PrepSCTFindMarkers(object = bmcd34_combined)
DefaultAssay(bmcd34_combined)<- "SCT"
feat_plt<- FeaturePlot(bmcd34_combined, features = c("S100A10","MPO","LTB","HLF","IRF8","ITGA2B","GATA1","CD79A"), label = T)
ggsave2("Feature_plot_topCluster_genes.pdf",feat_plt, height = 12, width = 10) 