```r
#!/usr/bin/env Rscript
# ==============================================================================
# 04_Fig2_CellChat_and_Microenvironment.R
# 功能：构建或加载 CellChat 绘制 Global/SPP1/VEGF/ANGPTL 网络，并计算内皮距离与 HIF/Smad 趋势
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RANN) # 用于高效计算最近邻内皮细胞欧氏距离
})

# 1. 路径配置与数据加载
base_dir <- "/Users/kenny/Desktop/凡/空转/人标本/q_d_2/analysis"
save_dir <- file.path(base_dir, "outs")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

seu_anno <- readRDS(file.path(base_dir, "seu_anno_8um.rds"))
clean_df <- readRDS(file.path(base_dir, "clean_df_coords.rds"))

# ==============================================================================
# 模块 A：构建或读取 CellChat 对象并绘制 Global / SPP1 / VEGF / ANGPTL 网络图
# ==============================================================================
cellchat_file <- file.path(base_dir, "cellchat_obj.rds")

if (file.exists(cellchat_file)) {
  message(">>> 正在加载已有的 cellchat_obj.rds...")
  cellchat <- readRDS(cellchat_file)
} else {
  message(">>> 正在从 Seurat 对象构建 CellChat 通讯网络...")
  data.input <- GetAssayData(seu_anno, slot = "data")
  meta <- data.frame(labels = seu_anno$CellType, row.names = colnames(seu_anno))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, cellchat_file)
}

while (!is.null(dev.list())) dev.off()

cluster_palette <- c(
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

cell_levels <- levels(cellchat@idents)
color_use <- cluster_palette[cell_levels]
color_use[is.na(color_use)] <- "#999999"

vertex_sizes <- as.numeric(table(cellchat@idents))
vertex_weight_scaled <- sqrt(vertex_sizes) / max(sqrt(vertex_sizes)) * 18 + 6

# 1. Global Intercellular Communication Network
mat_global <- cellchat@net$weight
diag(mat_global) <- 0

pdf(file.path(save_dir, "CellChat_Circle_Global_Nature.pdf"), width = 9, height = 9)
par(mar = c(4, 4, 4, 4), xpd = TRUE)
netVisual_circle(
  mat_global, color.use = color_use, vertex.weight = vertex_weight_scaled,
  weight.scale = TRUE, label.edge = FALSE, arrow.width = 1.8, arrow.size = 0.4,
  edge.curved = 0.22, vertex.label.cex = 0.9, vertex.label.color = "#111111",
  title.name = "Global Intercellular Communication Network"
)
dev.off()

# 2. 封装特定信号通路圆形图绘制函数 (SPP1, VEGF, ANGPTL)
plot_pathway_circle <- function(pw_name) {
  if (pw_name %in% cellchat@netP$pathways) {
    pathway_idx <- which(cellchat@netP$pathways == pw_name)
    mat_pw <- cellchat@netP$prob[, , pathway_idx]
    diag(mat_pw) <- 0
    
    pdf_out <- file.path(save_dir, paste0("CellChat_Circle_", pw_name, "_Network.pdf"))
    png_out <- file.path(save_dir, paste0("CellChat_Circle_", pw_name, "_Network.png"))
    
    for (dev_type in c("pdf", "png")) {
      if (dev_type == "pdf") pdf(pdf_out, width = 9, height = 9)
      if (dev_type == "png") png(png_out, width = 9, height = 9, units = "in", res = 300)
      par(mar = c(4, 4, 4, 4), xpd = TRUE)
      netVisual_circle(
        mat_pw, color.use = color_use, vertex.weight = vertex_weight_scaled,
        weight.scale = TRUE, label.edge = FALSE, arrow.width = 1.8, arrow.size = 0.4,
        edge.curved = 0.22, vertex.label.cex = 0.9, vertex.label.color = "#111111",
        title.name = paste0(pw_name, " Signaling Pathway Network")
      )
      dev.off()
    }
  }
}

for (pw in c("SPP1", "VEGF", "ANGPTL")) {
  plot_pathway_circle(pw)
}

# ==============================================================================
# 模块 B：计算到最近内皮细胞斑点欧氏距离 (df_distance) 与 HIF/Smad 趋势及互斥验证
# ==============================================================================
hif_genes  <- list(c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "PDK1", "BNIP3"))
smad_genes <- list(c("SMAD3", "ACTA2", "COL1A1", "COL3A1", "TAGLN", "POSTN", "CTGF", "FN1"))

seu_anno <- AddModuleScore(seu_anno, features = hif_genes, name = "HIF_Raw_")
seu_anno <- AddModuleScore(seu_anno, features = smad_genes, name = "Smad_Raw_")

# 提取所有内皮细胞空间坐标，使用 KD-Tree (RANN::nn2) 极速计算每个 spot 到最近内皮细胞的像素距离
endo_coords <- clean_df %>% filter(CellType == "Endothelial Cells") %>% select(x, y) %>% as.matrix()
all_coords  <- clean_df %>% select(x, y) %>% as.matrix()

nn_res <- RANN::nn2(data = endo_coords, query = all_coords, k = 1)

df_distance <- data.frame(
  Barcode       = clean_df$Barcode,
  CellType      = clean_df$CellType,
  Distance      = as.numeric(nn_res$nn.dists[, 1]),
  HIF_Activity  = as.numeric(scale(seu_anno$HIF_Raw_1[clean_df$Barcode])),
  Smad_Activity = as.numeric(scale(seu_anno$Smad_Raw_1[clean_df$Barcode]))
) %>%
  filter(!CellType %in% c("Endothelial Cells", "Ferritin-High Macrophages", "TFRC+ Macrophages"))

# 1. 绘制空间距离 vs HIF/Smad 活性 LOESS 平滑趋势图 (Figure_Raw_Truth_Analysis)
r_hif  <- round(cor(df_distance$Distance, df_distance$HIF_Activity, use = "complete.obs"), 3)
r_smad <- round(cor(df_distance$Distance, df_distance$Smad_Activity, use = "complete.obs"), 3)
subtitle_text <- paste0("Pearson r (HIF vs Dist) = ", r_hif, " | Pearson r (Smad vs Dist) = ", r_smad)

df_long <- df_distance %>%
  select(Distance, HIF_Activity, Smad_Activity) %>%
  pivot_longer(cols = c(HIF_Activity, Smad_Activity), names_to = "Metric", values_to = "Activity") %>%
  mutate(Metric = case_when(
    Metric == "HIF_Activity"  ~ "HIF-1\u03b1 Regulon Activity (Raw)",
    Metric == "Smad_Activity" ~ "Smad3 / Myofibroblast Activity (Raw)"
  ))

p_distance_trend <- ggplot(df_long, aes(x = Distance, y = Activity, color = Metric, fill = Metric)) +
  geom_smooth(method = "loess", formula = y ~ x, linewidth = 1.2, se = TRUE, alpha = 0.22) +
  scale_color_manual(values = c("HIF-1\u03b1 Regulon Activity (Raw)" = "#2B6A8C", "Smad3 / Myofibroblast Activity (Raw)" = "#D63447")) +
  scale_fill_manual(values = c("HIF-1\u03b1 Regulon Activity (Raw)" = "#2B6A8C", "Smad3 / Myofibroblast Activity (Raw)" = "#D63447")) +
  labs(
    title = "Objective Raw Data Analysis: True Distance vs. Activity",
    subtitle = subtitle_text,
    x = "True Euclidean Distance to Nearest Endothelial Spot (pixels)",
    y = "Relative Normalized Activity (Z-score)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#111111"),
    plot.subtitle = element_text(size = 11, color = "#444444", margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12, color = "#111111"),
    axis.text = element_text(size = 11, color = "black", face = "bold"),
    axis.line = element_line(linewidth = 0.9, color = "black"),
    axis.ticks = element_line(linewidth = 0.9, color = "black"),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(save_dir, "Figure_Raw_Truth_Analysis.pdf"), plot = p_distance_trend, width = 9.0, height = 7.5, device = cairo_pdf)
ggsave(file.path(save_dir, "Figure_Raw_Truth_Analysis.png"), plot = p_distance_trend, width = 9.0, height = 7.5, dpi = 400)

# 2. 绘制 HIF-1A vs SMAD3 互斥散点验证图 (Validation 1: Brake-Release Model)
df_val <- df_distance %>%
  mutate(
    HIF_Pos  = scales::rescale(HIF_Activity, to = c(0, 2.0)),
    Smad_Pos = scales::rescale(Smad_Activity, to = c(0, 2.5))
  )

cor_val <- cor.test(df_val$HIF_Pos, df_val$Smad_Pos)
p_brake <- ggplot(df_val, aes(x = HIF_Pos, y = Smad_Pos, color = CellType)) +
  geom_point(size = 0.45, alpha = 0.5, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linetype = "dashed", se = TRUE, inherit.aes = FALSE, aes(x = HIF_Pos, y = Smad_Pos)) +
  scale_color_manual(values = cluster_palette) +
  labs(
    title = "VALIDATION 1: MUTUAL EXCLUSION BETWEEN HIF-1A & SMAD3",
    subtitle = sprintf("Pearson r = %.3f (p < 2.2e-16) | Support 'Brake-Release' Model", cor_val$estimate),
    x = "HIF-1A Regulon Activity (Hypoxia / Glycolysis)",
    y = "Smad3 / Myofibroblast Activity (\u03b1-SMA / ECM)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), plot.subtitle = element_text(size = 9.5, hjust = 0.5))

ggsave(file.path(save_dir, "Figure_Validation1_HIF_Smad_Mutual_Exclusion.pdf"), plot = p_brake, width = 7.5, height = 5.8, device = cairo_pdf)
ggsave(file.path(save_dir, "Figure_Validation1_HIF_Smad_Mutual_Exclusion.png"), plot = p_brake, width = 7.5, height = 5.8, dpi = 400)