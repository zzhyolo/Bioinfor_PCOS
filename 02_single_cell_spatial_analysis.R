# ==============================================================================
# Script 02: Single-Cell RNA-Seq, Virtual Knockout, CellChat, and Spatial Transcriptomics
# Manuscript: Integrated network toxicology and multi-dimensional transcriptomics
#             reveal DEHP-associated pathogenic mechanisms in Polycystic Ovary Syndrome
# ==============================================================================

# ==============================================================================
# Part 1. Single-Cell RNA Sequencing (scRNA-seq) Processing & Harmony Integration
# ==============================================================================
library(Seurat)
library(harmony)
library(tidyverse)
library(RColorBrewer)

# 1.1 载入卵巢单细胞表达对象 (GSE255690, 6 normal ovarian tissue samples)
# sce.all <- readRDS("GSE255690_raw_seurat.rds")

# 标准质控、标准化、找高变基因与 PCA
# sce.all <- sce.all %>%
#   NormalizeData() %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(features = rownames(.)) %>%
#   RunPCA(npcs = 30)

# 1.2 Harmony 去除样本间批次效应 (dims = 1:20)
# sce.all <- sce.all %>%
#   RunHarmony(group.by.vars = "orig.ident", dims.use = 1:20) %>%
#   FindNeighbors(dims = 1:20, reduction = "harmony") %>%
#   FindClusters(resolution = 0.5) %>%
#   RunUMAP(dims = 1:20, reduction = "harmony")

# 1.3 细胞类型注释 (基于经典卵巢细胞标志基因列表: Zhang 2026 Aging Cell + Wu 2024 Nat Aging)
ovary_marker_list <- list(
  T_S    = c("DCN", "CLDN1", "FYB2", "STAR", "GSTA1", "CYP17A1"), # 膜/基质细胞
  GCs    = c("CYP19A1", "INHA", "AMH", "ACOXL"),                   # 颗粒细胞
  ECs    = c("VWF", "CLDN5", "TFF3"),                             # 内皮细胞
  SMCs   = c("RGS5", "TAGLN", "ACTA2", "MUSTN1"),                 # 平滑肌细胞
  MONO   = c("TYROBP", "IFI30"),                                  # 单核细胞
  NK     = c("CCL5", "NKG7"),                                     # NK细胞
  T_cell = c("IL7R", "KLRB1"),                                    # T细胞
  Oocyte = c("TUBB8", "ZP3", "FIGLA")                             # 卵母细胞
)

# 绘制 Marker 基因气泡图 (Figure 6B, 6D)
# DotPlot(sce.all, features = ovary_marker_list, group.by = "seurat_clusters") + RotatedAxis()

# Cluster 到 8 大细胞类型的映射字典
cluster_annot <- data.frame(
  cluster = as.character(0:17),
  celltype = c(
    "T_S", "T_S", "T_S", "ECs", "SMCs", "ECs", "T_S", "NK",
    "MONO", "SMCs", "T_cells", "T_cells", "GCs", "T_S", "GCs", "OO",
    "ECs", "T_S"
  )
)

# sce.all$celltype <- cluster_annot$celltype[match(sce.all$seurat_clusters, cluster_annot$cluster)]
# Idents(sce.all) <- "celltype"

# 绘制 UMAP 分群图 (Figure 6A, 6C) 与核心基因表达提琴图 (Figure 6E, 6F)
# DimPlot(sce.all, reduction = "umap", group.by = "celltype", label = TRUE)
# FeaturePlot(sce.all, features = c("APOE", "DHCR7"), reduction = "umap")


# ==============================================================================
# Part 2. Single-Cell Virtual Gene Knockout via scTenifoldKnk
# ==============================================================================
# library(scTenifoldNet)

# 针对颗粒细胞 (GCs) 进行 APOE 与 DHCR7 的计算机虚拟单细胞基因敲除模拟
run_vko <- function(count_matrix, target_gene) {
  result <- scTenifoldKnk(
    countMatrix    = count_matrix,
    qc             = TRUE,
    gKO            = target_gene,    # 敲除目标基因: "APOE" 或 "DHCR7"
    qc_mtThreshold = 0.1,
    qc_minLSize    = 6792,
    nc_lambda      = 0,
    nc_nNet        = 10,
    nc_nCells      = 6792,
    nc_nComp       = 3,
    nc_scaleScores = TRUE,
    nc_symmetric   = FALSE,
    nc_q           = 0.9,
    td_K           = 3,
    td_maxIter     = 1000,
    td_maxError    = 1e-05,
    td_nDecimal    = 3,
    ma_nDim        = 2
  )
  return(result)
}

# 示例:
# count_gc <- GetAssayData(subset(sce.all, idents = "GCs"), layer = "counts")
# vko_apoe  <- run_vko(count_gc, target_gene = "APOE")
# vko_dhcr7 <- run_vko(count_gc, target_gene = "DHCR7")


# ==============================================================================
# Part 3. Cell–Cell Communication Analysis (CellChat)
# ==============================================================================
library(CellChat)
library(patchwork)

# 为保证计算效率与代表性，随机抽取 25,000 个细胞
# set.seed(123)
# selected_cells <- sample(colnames(sce.all), size = min(25000, ncol(sce.all)), replace = FALSE)
# scRNA <- sce.all[, selected_cells]

# 准备表达矩阵与元数据
# data.input <- GetAssayData(scRNA, layer = "data")
# meta <- scRNA@meta.data[, c("orig.ident", "celltype")]
# colnames(meta) <- c("group", "celltypes")

# 创建 CellChat 对象
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltypes")
# cellchat@DB <- CellChatDB.human # 使用人全配受体库

# 3.1 识别过表达基因及相互作用
# cellchat <- subsetData(cellchat)
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)

# 3.2 计算细胞间通讯概率 (truncatedMean 模式, min.cells = 1)
# cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE)
# cellchat <- filterCommunication(cellchat, min.cells = 1)

# 3.3 计算信号通路水平通讯概率与整合网络 (Figure 7A, 7B)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)

# 3.4 提取关注的核心相互作用轴 (APOE-TREM2, desmosterol-NR1H2, cholesterol-RORA)
# df.net <- subsetCommunication(cellchat)
# core_targets <- c("DHCR7", "APOE", "TREM2", "NR1H2", "RORA")
# df.net_core <- df.net %>% filter(str_detect(interaction_name, paste(core_targets, collapse = "|")))

# 绘制气泡图展示核心配受体对通信概率与显著性 (Figure 7D)
# netVisual_bubble(cellchat, pairLR.use = df.net_core, remove.isolate = TRUE)


# ==============================================================================
# Part 4. Human Ovarian Spatial Transcriptomics (Stereo-seq) & semla Co-Localization
# ==============================================================================
library(semla)

# 4.1 载入并处理人类正常卵巢空间转录组对象 (GSE267315)
# seu_st <- qs2::qs_read("st_ova_sct_12.qs2")

# SCTransform 标准化并降维聚类 (dims = 1:30, resolution = 0.2)
# seu_st <- SCTransform(seu_st, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
# seu_st <- RunPCA(seu_st, assay = "SCT", verbose = FALSE)
# seu_st <- FindNeighbors(seu_st, reduction = "pca", dims = 1:30)
# seu_st <- FindClusters(seu_st, resolution = 0.2, verbose = FALSE)
# seu_st <- RunUMAP(seu_st, reduction = "pca", dims = 1:30)

# 4.2 空间特征分布图 (Spatial Feature Plots, Figure 8A, 8C, 8E)
# SpatialPlot(seu_st, label = FALSE)

# 4.3 转换至 semla 绘制高分辨率多基因空间共定位图 (Figure 8C, 8E)
# colnames(seu_st@images$image@coordinates) <- c("x", "y")
# SCC <- UpdateSeuratForSemla(seu_st)

# (1) APOE 与 TREM2 的空间共定位分布图
# MapMultipleFeatures(SCC,
#                     features = c("APOE", "TREM2"),
#                     pt_size = 0.9, max_cutoff = 0.99,
#                     override_plot_dims = TRUE,
#                     colors = c("#ccebc5", "#e41a1c"),
#                     scale = "free")

# (2) DHCR7, RORA 与 NR1H2 的空间共定位分布图
# MapMultipleFeatures(SCC,
#                     features = c("DHCR7", "RORA", "NR1H2"),
#                     pt_size = 0.8, max_cutoff = 0.99, min_cutoff = 0.2,
#                     override_plot_dims = TRUE,
#                     colors = c("#e41a1c", "#984ea3", "#a65628"))
