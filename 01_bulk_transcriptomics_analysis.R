# ==============================================================================
# Script 01: Bulk Transcriptomics and Network Toxicology Analysis Pipeline
# Manuscript: Integrated network toxicology and multi-dimensional transcriptomics
#             reveal DEHP-associated pathogenic mechanisms in Polycystic Ovary Syndrome
# ==============================================================================

# ==============================================================================
# Part 1. Multi-Cohort Integration & Batch Effect Correction (ComBat)
# ==============================================================================
library(sva)
library(FactoMineR)
library(factoextra)
library(tidyverse)

# 载入原始未去批次的表达矩阵及表型数据 (含 GSE138518, GSE193123, GSE168404)
# load("pcos_data.rda") # data_all contains exp_all and pd_data

exp_all <- data_all$exp_all
pd_all  <- data_all$pd_data

# 1.1 绘制去批次前的 PCA 图
dat.pca_before <- PCA(as.data.frame(t(exp_all)), graph = FALSE)
pca_before_plot <- fviz_pca_ind(dat.pca_before,
                                geom.ind = "points",
                                col.ind = pd_all$batch,
                                palette = c("#00AFBB", "#E7B800", "navy"),
                                addEllipses = TRUE,
                                legend.title = "Batch")
print(pca_before_plot)

# 1.2 使用 sva::ComBat 去除批次效应
pd_all$batch <- as.factor(pd_all$batch)
exp_all_combat <- ComBat(dat = exp_all, batch = pd_all$batch)

# 1.3 绘制去批次后的 PCA 诊断图 (Supplementary Figure S1B)
dat.pca_after <- PCA(as.data.frame(t(exp_all_combat)), graph = FALSE)
pca_after_plot <- fviz_pca_ind(dat.pca_after,
                               geom.ind = "points",
                               col.ind = pd_all$batch,
                               palette = c("#00AFBB", "#E7B800", "navy"),
                               addEllipses = TRUE,
                               legend.title = "Batch")
print(pca_after_plot)

data_merge <- list(exp_all = exp_all_combat, pd_all = pd_all)
save(data_merge, file = "pcos_merge.rda")


# ==============================================================================
# Part 2. Differential Expression Analysis (limma) & Volcano Plot
# ==============================================================================
library(limma)
library(data.table)

load("pcos_merge.rda")

exp <- data_merge$exp_all
exprSet <- as.data.frame(exp)
group_list <- factor(data_merge$pd_all$group, levels = c("Normal", "Disease"))

cat("Sample counts per group:\n")
print(table(group_list))

# 2.1 构建线性模型与对比矩阵
design <- model.matrix(~ 0 + group_list)
colnames(design) <- c("Normal", "Disease")
rownames(design) <- colnames(exprSet)

contrast.matrix <- makeContrasts(Disease - Normal, levels = design)

# 2.2 运行经验贝叶斯检验
fit <- lmFit(exprSet, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
allDiff <- topTable(fit, number = Inf, adjust.method = "BH")

write.csv(allDiff, file = "pcos_allDiff.csv")

# 2.3 依据预设阈值划分差异表达基因 (|log2FC| > 0.2, p-value < 0.05)
logFC_cutoff <- 0.2
pvalue_cutoff <- 0.05

allDiff$change <- "stable"
allDiff$change[allDiff$P.Value < pvalue_cutoff & allDiff$logFC > logFC_cutoff]  <- "up"
allDiff$change[allDiff$P.Value < pvalue_cutoff & allDiff$logFC < -logFC_cutoff] <- "down"
allDiff$symbol <- rownames(allDiff)

cat("\nDEG counts (|logFC| > 0.2 & P < 0.05):\n")
print(table(allDiff$change))

DEG_limma <- allDiff
save(DEG_limma, file = "DEG_limma.rda")

# 2.4 绘制火山图 (Figure 2A)
p_volcano <- ggplot(data = allDiff, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.5, size = 2.5, aes(color = change)) +
  scale_color_manual(values = c("down" = "#2166AC", "stable" = "grey75", "up" = "#B2182B")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4, col = "black", lwd = 0.6) +
  geom_hline(yintercept = -log10(pvalue_cutoff), lty = 4, col = "black", lwd = 0.6) +
  labs(x = "log2(Fold Change)", y = "-log10(P-value)", title = "Volcano Plot (PCOS vs Normal)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_volcano)


# ==============================================================================
# Part 3. Target Intersection & Functional Enrichment Analysis (GO / KEGG)
# ==============================================================================
library(org.Hs.eg.db)
library(clusterProfiler)

# 3.1 提取显著差异基因与 DEHP 预测靶点交集
load("DEG_limma.rda")
degs <- DEG_limma %>% filter(change != "stable")

# 读取公共数据库 (CTD, STITCH, SwissTargetPrediction, ChEMBL) 汇总靶点
if (file.exists("toxic_targets.csv")) {
  toxic_targets <- fread("toxic_targets.csv", header = TRUE)
} else if (file.exists("toxic_targets_merge/toxic_targets.csv")) {
  toxic_targets <- fread("toxic_targets_merge/toxic_targets.csv", header = TRUE)
} else {
  # 若从补充材料读取
  library(openxlsx)
  toxic_targets <- read.xlsx("C:/Users/Administrator/Desktop/revise_bbrc/supp_materials.xlsx",
                             sheet = "Supplementary Table 2")
  colnames(toxic_targets) <- c("Index", "x")
}

target_col <- ifelse("x" %in% colnames(toxic_targets), "x", colnames(toxic_targets)[2])
dehp_targets <- unique(as.character(toxic_targets[[target_col]]))

toxic_pcos_genes <- intersect(degs$symbol, dehp_targets)
cat("Overlapping candidate targets (Figure 2B):", length(toxic_pcos_genes), "\n")

save(toxic_pcos_genes, degs, file = "deg_intersect.rda")
write.csv(data.frame(Gene = toxic_pcos_genes), file = "DEHP_PCOS_overlapping_genes.csv", row.names = FALSE)

# 3.2 转换 Entrez ID
gene_ids <- bitr(toxic_pcos_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 3.3 GO 富集分析 (BP, CC, MF)
GO_all <- enrichGO(gene          = gene_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "ALL",
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   readable      = TRUE)

GO_result <- as.data.frame(GO_all)
write.csv(GO_result, file = "GO_result.csv", row.names = FALSE)

# 3.4 KEGG 通路富集分析
kk <- enrichKEGG(gene         = gene_ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 minGSSize    = 5,
                 maxGSSize    = 500)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
kegg_result <- as.data.frame(kk)
write.csv(kegg_result, file = "kegg_result.csv", row.names = FALSE)


# ==============================================================================
# Part 4. Core Target Prioritization via Machine Learning (LASSO & Random Forest)
# ==============================================================================
library(glmnet)
library(randomForest)
library(caret)

load("pcos_merge.rda")
load("deg_intersect.rda")

exp_matrix <- data_merge$exp_all
common_genes <- intersect(rownames(exp_matrix), toxic_pcos_genes)

# 构建样本-特征矩阵 (转置后: 样本 x 基因)
mat_features <- as.data.frame(t(exp_matrix[common_genes, ]))
group_vec <- factor(data_merge$pd_all$group, levels = c("Normal", "Disease"))
y_numeric <- ifelse(group_vec == "Disease", 1, 0)
x_mat <- as.matrix(mat_features)

# 4.1 LASSO 回归筛选 (10-fold CV, 1-SE rule)
set.seed(12345)
cv_lasso <- cv.glmnet(x = x_mat, y = y_numeric, family = "binomial", alpha = 1,
                      nfolds = 10, type.measure = "deviance")

best_lambda <- cv_lasso$lambda.1se
coef_mat <- coef(cv_lasso, s = "lambda.1se")
lasso_genes <- rownames(coef_mat)[which(coef_mat != 0)]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]
cat("LASSO prioritized genes (lambda.1se):", length(lasso_genes), "\n")
print(lasso_genes)

# 4.2 随机森林特征重要性筛选 (ntree = 100, mtry = 4, seed = 123)
set.seed(123)
rf_data <- mat_features
rf_data$group <- group_vec

rf_model <- randomForest(group ~ ., data = rf_data, ntree = 100, mtry = 4, importance = TRUE)

imp_matrix <- importance(rf_model)
imp_df <- data.frame(Gene = rownames(imp_matrix), MeanDecreaseGini = imp_matrix[, "MeanDecreaseGini"])
imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), ]
rf_genes <- head(imp_df$Gene, 15)
cat("Random Forest Top 15 genes:\n")
print(rf_genes)

# 4.3 取两算法并集获得 22 个核心靶点
core_candidate_genes <- union(lasso_genes, rf_genes)
cat("Core candidate genes prioritized by ML (Union = 22):", length(core_candidate_genes), "\n")
print(core_candidate_genes)

save(lasso_genes, rf_genes, core_candidate_genes, file = "machine_learning_core_genes.rda")
write.csv(data.frame(CoreGene = core_candidate_genes), file = "core_candidate_genes_22.csv", row.names = FALSE)


# ==============================================================================
# Part 5. Consensus Clustering and GSEA Analysis
# ==============================================================================
library(ConsensusClusterPlus)
library(msigdbr)
library(pheatmap)

load("pcos_merge.rda")
load("machine_learning_core_genes.rda")

exp_all <- data_merge$exp_all
common_core <- intersect(rownames(exp_all), core_candidate_genes)

expr_core <- exp_all[common_core, ]
expr_centered <- sweep(expr_core, 1, apply(expr_core, 1, median, na.rm = TRUE))

# 5.1 一致性聚类 (k = 2 to 6, 500 reps)
out_dir <- "ConsensusCluster"
dir.create(out_dir, showWarnings = FALSE)

results <- ConsensusClusterPlus(
  d             = expr_centered,
  maxK          = 6,
  reps          = 500,
  pItem         = 0.8,
  pFeature      = 1,
  title         = out_dir,
  clusterAlg    = "hc",
  distance      = "pearson",
  seed          = 123456,
  plot          = "pdf"
)

consensus_class_k2 <- results[[2]]$consensusClass
cat("Sample distribution in Consensus Clusters (k = 2):\n")
print(table(consensus_class_k2))

# 5.2 绘制双向层次聚类热图 (Figure 5B)
annotation_col <- data.frame(
  Group   = data_merge$pd_all$group,
  Cluster = paste0("Cluster_", consensus_class_k2),
  row.names = colnames(expr_core)
)

pheatmap(expr_centered,
         annotation_col = annotation_col,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_colnames = FALSE,
         fontsize_row = 9,
         filename = "Figure5B_Heatmap.pdf",
         width = 7, height = 6)

# 5.3 GSEA 分析 (基于 MSigDB Hallmark 基因集)
load("DEG_limma.rda")
gene_list <- DEG_limma$logFC
names(gene_list) <- rownames(DEG_limma)
gene_list <- sort(gene_list, decreasing = TRUE)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

set.seed(12345)
gsea_res <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = hallmark_sets,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 5,
  maxGSSize     = 500
)

gsea_df <- as.data.frame(gsea_res)
cat("Significant Hallmark pathways identified by GSEA:", nrow(gsea_df), "\n")
write.csv(gsea_df, file = "GSEA_Hallmark_results.csv", row.names = FALSE)
save(gsea_res, file = "GSEA_results.rda")
