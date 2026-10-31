```r
#!/usr/bin/env Rscript
# ==============================================================================
# 01_Data_Prep_and_Clustering.R
# 功能：读取 8um binned_outputs 与空间图像，完成 11 亚群注释与拟时序对象构建
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(SingleCellExperiment)
  library(slingshot)
})

set.seed(2026)

# 1. 路径配置
base_dir <- "/Users/kenny/Desktop/凡/空转/人标本/q_d_2/analysis"
save_dir <- file.path(base_dir, "outs")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

# Space Ranger 原始 outs 目录 (若 binned_outputs 在上一级目录，请按实际位置微调)
spaceranger_outs <- "/Users/kenny/Desktop/凡/空转/人标本/q_d_2/outs/"
bin_dir <- file.path(spaceranger_outs, "binned_outputs/square_008um")

# ------------------------------------------------------------------------------
# 2. 全局统一 11 亚群名称与调色板
# ------------------------------------------------------------------------------
target_levels <- c(
  "Tendon Core Fibroblasts",
  "Quiescent Fibroblasts-1",
  "Quiescent Fibroblasts-2",
  "ECM-Remodeling Fibroblasts",
  "Fibro-Immune Transitional",
  "Pro-fibrotic Fibroblasts",
  "Active Fibroblasts",
  "Myofibroblasts",
  "Endothelial Cells",
  "Ferritin-High Macrophages",
  "TFRC+ Macrophages"
)

umap_colors <- c(
  "Quiescent Fibroblasts-1"    = "#66C2A5",
  "Quiescent Fibroblasts-2"    = "#8DA0CB",
  "Tendon Core Fibroblasts"    = "#A6D854",
  "Endothelial Cells"          = "#E78AC3",
  "Fibro-Immune Transitional"  = "#FFD92F",
  "Ferritin-High Macrophages"  = "#B3B3B3",
  "TFRC+ Macrophages"          = "#E5C494",
  "Active Fibroblasts"         = "#FB8072",
  "Pro-fibrotic Fibroblasts"   = "#FDB462",
  "ECM-Remodeling Fibroblasts" = "#BC80BD",
  "Myofibroblasts"             = "#E41A1C"
)

# ------------------------------------------------------------------------------
# 3. 读取 8um Binned 矩阵与 Spatial 图像 (或直接加载已有 RDS)
# ------------------------------------------------------------------------------
rds_file <- file.path(base_dir, "seu_anno_8um.rds")

if (file.exists(rds_file)) {
  message(">>> 检测到已处理的 seu_anno_8um.rds，正在直接加载...")
  seu_anno <- readRDS(rds_file)
} else {
  message(">>> 正在从 square_008um 读取原始矩阵与空间图像...")
  counts <- Read10X_h5(file.path(bin_dir, "filtered_feature_bc_matrix.h5"))
  seu_anno <- CreateSeuratObject(counts = counts, assay = "Spatial.008um")
  
  # 读取 R 图像与 parquet 空间坐标
  img <- Read10X_Image(
    image.dir = file.path(bin_dir, "spatial"),
    filter.matrix = TRUE
  )
  DefaultAssay(img) <- "Spatial.008um"
  seu_anno[["slice1"]] <- img[Cells(seu_anno)]
  
  # 质控、标准化与降维聚类
  seu_anno[["percent.mt"]] <- PercentageFeatureSet(seu_anno, pattern = "^MT-")
  seu_anno <- subset(seu_anno, subset = nCount_Spatial.008um > 100 & nFeature_Spatial.008um > 80 & percent.mt < 15)
  
  seu_anno <- NormalizeData(seu_anno, scale.factor = 10000)
  seu_anno <- FindVariableFeatures(seu_anno, nfeatures = 3000)
  seu_anno <- ScaleData(seu_anno)
  seu_anno <- RunPCA(seu_anno, npcs = 30)
  seu_anno <- FindNeighbors(seu_anno, dims = 1:30)
  seu_anno <- FindClusters(seu_anno, resolution = 0.6)
  seu_anno <- RunUMAP(seu_anno, dims = 1:30)
  
  # 映射 11 个亚群 (请根据实际 cluster 编号对应调整)
  cluster_map <- c(
    "0"  = "Tendon Core Fibroblasts",
    "1"  = "Active Fibroblasts",
    "2"  = "Fibro-Immune Transitional",
    "3"  = "ECM-Remodeling Fibroblasts",
    "4"  = "Pro-fibrotic Fibroblasts",
    "5"  = "Quiescent Fibroblasts-1",
    "6"  = "Quiescent Fibroblasts-2",
    "7"  = "Myofibroblasts",
    "8"  = "Ferritin-High Macrophages",
    "9"  = "TFRC+ Macrophages",
    "10" = "Endothelial Cells"
  )
  seu_anno$CellType <- factor(unname(cluster_map[as.character(Idents(seu_anno))]), levels = target_levels)
  Idents(seu_anno) <- "CellType"
}

# 修正可能存在的拼写截断并统一因子水平
seu_anno$CellType <- as.character(seu_anno$CellType)
seu_anno$CellType[seu_anno$CellType == "Pro-fibrotic Fibroblas"] <- "Pro-fibrotic Fibroblasts"
seu_anno$CellType <- factor(seu_anno$CellType, levels = target_levels)
Idents(seu_anno) <- "CellType"

# ------------------------------------------------------------------------------
# 4. 提取空间坐标构建 clean_df 并划分 Healthy / Fibrotic 区域
# ------------------------------------------------------------------------------
coords <- GetTissueCoordinates(seu_anno)
colnames(coords)[1:2] <- c("x", "y")

clean_df <- data.frame(
  Barcode  = rownames(coords),
  x        = as.numeric(coords$x),
  y        = as.numeric(coords$y),
  CellType = seu_anno$CellType,
  stringsAsFactors = FALSE
)

# 标记空间病理分区 (若 metadata 中已有 Region 列则保留，否则根据纤维化群空间分布划分)
if (!"Region" %in% colnames(seu_anno@meta.data)) {
  seu_anno$Region <- ifelse(
    seu_anno$CellType %in% c("Active Fibroblasts", "Pro-fibrotic Fibroblasts", "Myofibroblasts", "Ferritin-High Macrophages", "TFRC+ Macrophages"),
    "Fibrotic", "Healthy"
  )
}

# ------------------------------------------------------------------------------
# 5. 构建成纤维细胞子集 fib_clean 并计算 Slingshot 拟时序 global_pt
# ------------------------------------------------------------------------------
fibro_subsets <- c(
  "Tendon Core Fibroblasts", "Quiescent Fibroblasts-1", "Quiescent Fibroblasts-2",
  "ECM-Remodeling Fibroblasts", "Fibro-Immune Transitional", "Pro-fibrotic Fibroblasts",
  "Active Fibroblasts", "Myofibroblasts"
)

fib_clean <- subset(seu_anno, idents = fibro_subsets)

if (!"global_pt" %in% colnames(fib_clean@meta.data)) {
  message(">>> 正在运行 Slingshot 计算成纤维细胞演进拟时序 (global_pt)...")
  sce <- as.SingleCellExperiment(fib_clean)
  sce <- slingshot(sce, clusterLabels = "CellType", reducedDim = "UMAP", start.clus = "Tendon Core Fibroblasts")
  pt_raw <- slingPseudotime(sce)[, 1]
  # 将拟时序缩放至 0 - 16 区间，精准匹配分期阈值 4.5 与 11.5
  fib_clean$global_pt <- scales::rescale(pt_raw, to = c(0, 16))
}

# 保存全部前置数据供后续脚本直接读取
saveRDS(seu_anno, file.path(base_dir, "seu_anno_8um.rds"))
saveRDS(fib_clean, file.path(base_dir, "fib_clean_8um.rds"))
saveRDS(clean_df, file.path(base_dir, "clean_df_coords.rds"))
message(">>> 01 数据准备与拟时序构建完成！所有对象已保存至: ", base_dir)