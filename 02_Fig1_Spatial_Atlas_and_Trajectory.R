```r
#!/usr/bin/env Rscript
# ==============================================================================
# 02_Fig1_Spatial_Atlas_and_Trajectory.R
# 功能：整合绘制 11 亚群空间高亮图、Marker DotPlot、8 亚群山峦图及 4 模块时序热图
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(scales)
  library(ggridges)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# 1. 设置路径并加载 01 生成的数据
base_dir <- "/Users/kenny/Desktop/凡/空转/人标本/q_d_2/analysis"
data_dir <- file.path(base_dir, "outs")
save_dir <- file.path(base_dir, "outs")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

seu_anno  <- readRDS(file.path(base_dir, "seu_anno_8um.rds"))
fib_clean <- readRDS(file.path(base_dir, "fib_clean_8um.rds"))
clean_df  <- readRDS(file.path(base_dir, "clean_df_coords.rds"))

target_levels <- levels(seu_anno$CellType)
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

# ==============================================================================
# 模块 A：全局空间分布图、UMAP 与 Marker Gene DotPlot
# ==============================================================================
p_spatial_all <- ggplot(clean_df, aes(x = x, y = -y, color = CellType)) +
  geom_point(size = 0.45, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = umap_colors) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(save_dir, "Fig1C_Global_Spatial_Distribution.png"), p_spatial_all, width = 9, height = 8, dpi = 400)

p_umap_all <- DimPlot(seu_anno, reduction = "umap", group.by = "CellType", cols = umap_colors, label = TRUE, repel = TRUE, pt.size = 0.25) +
  theme_classic()
ggsave(file.path(save_dir, "Fig1D_Global_UMAP.pdf"), p_umap_all, width = 8.5, height = 6.5)

# Marker Gene DotPlot
gene_groups <- list(
  "Tendon ECM"             = c("COL1A1", "COL1A2", "DCN", "SCX", "FMOD", "COL4A1", "LUM"),
  "Fibrotic / Contractile" = c("ACTA2", "TAGLN", "POSTN", "COL3A1", "FAP"),
  "ECM Remodeling"         = c("MMP2", "MMP14", "TIMP1", "LOX", "LOXL2"),
  "Iron Homeostasis"       = c("TFRC", "FTH1", "FTL", "GPX4", "HMOX1"),
  "Endothelial"            = c("PECAM1", "VWF", "CDH5"),
  "Chondrogenic"           = c("COMP", "CRTAC1", "MMP3"),
  "Myeloid / Resorption"   = c("CD68", "SPP1", "ACP5", "MMP9", "APOE")
)

gene2group <- do.call(rbind, lapply(names(gene_groups), function(grp) {
  data.frame(gene = gene_groups[[grp]], Group = grp, stringsAsFactors = FALSE)
})) %>% filter(gene %in% rownames(seu_anno))
gene2group$Group <- factor(gene2group$Group, levels = names(gene_groups))

ordered_levels <- c(
  "Active Fibroblasts", "Myofibroblasts", "Pro-fibrotic Fibroblasts", "ECM-Remodeling Fibroblasts",
  "Tendon Core Fibroblasts", "Quiescent Fibroblasts-1", "Quiescent Fibroblasts-2",
  "Fibro-Immune Transitional", "Endothelial Cells",
  "Ferritin-High Macrophages", "TFRC+ Macrophages"
)
ordered_levels <- intersect(ordered_levels, unique(as.character(seu_anno$CellType)))
seu_anno$CellType <- factor(seu_anno$CellType, levels = ordered_levels)

p_base <- DotPlot(seu_anno, features = gene2group$gene, group.by = "CellType")
plot_data <- p_base$data %>%
  rename(gene = features.plot) %>%
  mutate(gene = as.character(gene)) %>%
  inner_join(gene2group, by = "gene")

plot_data$id <- factor(plot_data$id, levels = rev(ordered_levels))
plot_data$gene <- factor(plot_data$gene, levels = gene2group$gene)

p_final_dot <- ggplot(plot_data, aes(x = gene, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, color = "#333333", stroke = 0.35) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_size(range = c(0.8, 5.5), breaks = c(0, 25, 50, 75, 100), limits = c(0, 100), name = "Percent Expressed") +
  scale_fill_gradient2(low = "#74A9CF", mid = "#F7F7F7", high = "#E26D8C", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Avg Expr (z-score)") +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black", face = "italic"),
    axis.text.y = element_text(size = 9.5, color = "black", face = "bold"),
    panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#999999", linewidth = 0.5),
    strip.background = element_rect(fill = "#D8D8D8", color = "#999999", linewidth = 0.5),
    strip.text = element_text(size = 9, face = "bold", color = "#222222", margin = margin(t = 4, b = 4)),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "right"
  )

out_dot <- file.path(save_dir, "Figure_Markers_DotPlot_Final_Publication")
ggsave(paste0(out_dot, ".pdf"), plot = p_final_dot, device = cairo_pdf, width = 16.5, height = 6.4)
ggsave(paste0(out_dot, ".png"), plot = p_final_dot, width = 16.5, height = 6.4, dpi = 400)

# ==============================================================================
# 模块 B：11 分群空间高亮对比图 (Figure_All_11_Subclusters_Highlighted)
# ==============================================================================
focus_types <- target_levels
facet_list <- lapply(focus_types, function(ct) {
  df <- clean_df
  df$Focus_Target <- ct
  df$Is_Focus <- (df$CellType == ct)
  return(df)
})
facet_df <- bind_rows(facet_list)
facet_df$Focus_Target <- factor(facet_df$Focus_Target, levels = focus_types)

p_focus_11 <- ggplot() +
  geom_point(data = subset(facet_df, !Is_Focus), aes(x = x, y = -y), color = "#E9ECEF", size = 0.4, alpha = 0.4, stroke = 0) +
  geom_point(data = subset(facet_df, Is_Focus), aes(x = x, y = -y, color = Focus_Target), size = 0.7, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = umap_colors) +
  facet_wrap(~Focus_Target, ncol = 4) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 9, face = "bold", color = "#111111", margin = margin(b = 4)),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(file.path(save_dir, "Figure_All_11_Subclusters_Highlighted.png"), plot = p_focus_11, width = 14, height = 10, dpi = 400)

# ==============================================================================
# 模块 C：8 个成纤维亚群时序山峦图 (Fig2_Human_Tendon_RidgePlot_Flawless)
# ==============================================================================
df_all8_raw <- data.frame(
  CellType = as.character(fib_clean$CellType),
  Pseudotime = fib_clean$global_pt
) %>%
  mutate(CellType = ifelse(CellType == "Pro-fibrotic Fibroblas", "Pro-fibrotic Fibroblasts", CellType)) %>%
  filter(!is.na(Pseudotime)) %>%
  filter(Pseudotime <= quantile(Pseudotime, 0.985))

ordered_stats <- df_all8_raw %>%
  group_by(CellType) %>%
  summarise(med_pt = median(Pseudotime)) %>%
  arrange(med_pt)

sorted_subsets <- ordered_stats$CellType
df_all8 <- df_all8_raw %>% mutate(CellType = factor(CellType, levels = rev(sorted_subsets)))

fresh_palette_dynamic <- setNames(
  colorRampPalette(c("#3A7D8C", "#55A88B", "#78BA89", "#F3A953", "#E8835C", "#DC6C73", "#BA557E", "#8B3E69"))(length(sorted_subsets)),
  sorted_subsets
)

x_min <- min(df_all8$Pseudotime) - 0.4
x_max <- max(df_all8$Pseudotime) + 0.4

df_phase <- data.frame(
  Phase = factor(c("Phase I: Homeostasis", "Phase II: Plasticity Window", "Phase III: Fibrogenesis"),
                 levels = c("Phase I: Homeostasis", "Phase II: Plasticity Window", "Phase III: Fibrogenesis")),
  xmin = c(x_min, 4.5, 11.5),
  xmax = c(4.5, 11.5, x_max),
  ymin = 0.08, ymax = 0.92,
  label_x = c((x_min + 4.5) / 2, (4.5 + 11.5) / 2, (11.5 + x_max) / 2)
)

phase_box_fills <- c(
  "Phase I: Homeostasis"        = "#E6F2FA",
  "Phase II: Plasticity Window" = "#FFF3E6",
  "Phase III: Fibrogenesis"     = "#FDECEE"
)

p_phase_track <- ggplot(df_phase) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Phase), color = "#333333", linewidth = 0.65) +
  geom_text(aes(x = label_x, y = 0.5, label = Phase), fontface = "bold", size = 3.3, color = "#222222") +
  scale_fill_manual(values = phase_box_fills) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(t = 4, r = 18, b = 2, l = 10))

p_ridge_main <- ggplot(df_all8, aes(x = Pseudotime, y = CellType, fill = CellType)) +
  geom_vline(xintercept = c(4.5, 11.5), linetype = "dashed", color = "#B0B0B0", linewidth = 0.45, alpha = 0.75) +
  geom_density_ridges(
    scale = 1.85, alpha = 0.82, color = "#2A2A2A", linewidth = 0.55, rel_min_height = 0.015,
    quantile_lines = TRUE, quantiles = 2, vline_color = "#111111", vline_linetype = "solid", vline_width = 0.6
  ) +
  scale_fill_manual(values = fresh_palette_dynamic) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 18, 3), limits = c(x_min, x_max)) +
  theme_classic() +
  labs(x = "Pseudotime (Progression to Fibrosis)", y = "") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#333333", linewidth = 0.55),
    axis.title.x = element_text(face = "bold", size = 10.5, color = "#222222", margin = margin(t = 6)),
    axis.text.y = element_text(face = "bold", size = 9.2, color = "#222222"),
    axis.text.x = element_text(face = "bold", size = 9.0, color = "#333333"),
    axis.line = element_line(linewidth = 0.55, color = "#333333"),
    plot.margin = margin(t = 2, r = 18, b = 8, l = 10)
  )

title_theme <- ggdraw() + draw_label("Fibroblast Lineage Transition Dynamics in Human Rotator Cuff Tendon", fontface = "bold", size = 12.5, x = 0.5, y = 0.65, color = "#111111")
p_final_flawless <- plot_grid(title_theme, p_phase_track, p_ridge_main, ncol = 1, rel_heights = c(0.08, 0.07, 0.85), align = "v", axis = "lr")

ggsave(file.path(save_dir, "Fig2_Human_Tendon_RidgePlot_Flawless.pdf"), plot = p_final_flawless, width = 7.8, height = 5.6)
ggsave(file.path(save_dir, "Fig2_Human_Tendon_RidgePlot_Flawless.png"), plot = p_final_flawless, width = 7.8, height = 5.6, dpi = 300)

# ==============================================================================
# 模块 D：四大模块拟时序动态基因热图 (Fig2D_Pseudotime_Heatmap_AllGenes_Labeled)
# ==============================================================================
valid_cells <- rownames(fib_clean@meta.data)[!is.na(fib_clean$global_pt)]
pt_vals <- fib_clean$global_pt[valid_cells]
pt_cutoff <- quantile(pt_vals, 0.985)
keep_cells <- valid_cells[pt_vals <= pt_cutoff]
pt_vec <- fib_clean$global_pt[keep_cells]

m1_homeostasis  <- c("Scx", "Mkx", "Tnmd", "Dcn", "Fmod", "Col1a1", "Col1a2", "Comp", "Bgn", "Lum")
m2_iron_stress  <- c("Fth1", "Ftl", "Ftl1", "Tfrc", "Slc39a14", "Slc40a1", "Hmox1", "Acsl4", "Ptgs2", "Sod2", "Prdx1")
m3_plasticity   <- c("Lox", "Loxl1", "Loxl2", "Bmp1", "Pcolce", "Postn", "Mmp2", "Timp1", "Timp2", "Cthrc1", "Sparc", "Fn1", "Ccn2", "Ctgf")
m4_fibrogenesis <- c("Acta2", "Tagln", "Col3a1", "Col5a1", "Col5a2", "Tgfb1", "Serpine1", "Myl9", "Tpm2", "Cnn1")

all_features <- rownames(fib_clean)
match_genes <- function(gene_list) {
  matched <- character()
  for (g in gene_list) {
    hit <- all_features[tolower(all_features) == tolower(g)]
    if (length(hit) > 0) matched <- c(matched, hit[1])
  }
  return(unique(matched))
}

m1_genes <- match_genes(m1_homeostasis)
m2_genes <- match_genes(m2_iron_stress)
m3_genes <- match_genes(m3_plasticity)
m4_genes <- match_genes(m4_fibrogenesis)

candidate_genes <- unique(c(m1_genes, m2_genes, m3_genes, m4_genes))
expr_raw <- FetchData(fib_clean, vars = candidate_genes, cells = keep_cells)

gene_nonzero <- colnames(expr_raw)[colSums(expr_raw > 0) >= 5]
m1_genes <- intersect(m1_genes, gene_nonzero)
m2_genes <- intersect(m2_genes, gene_nonzero)
m3_genes <- intersect(m3_genes, gene_nonzero)
m4_genes <- intersect(m4_genes, gene_nonzero)

final_genes <- c(m1_genes, m2_genes, m3_genes, m4_genes)
gene_modules <- factor(
  c(rep("Module 1: Physiological Homeostasis", length(m1_genes)),
    rep("Module 2: Iron Accumulation & Redox Stress", length(m2_genes)),
    rep("Module 3: Plasticity & Crosslinking (Target)", length(m3_genes)),
    rep("Module 4: Terminal Fibrogenesis & Scarring", length(m4_genes))),
  levels = c("Module 1: Physiological Homeostasis", "Module 2: Iron Accumulation & Redox Stress",
             "Module 3: Plasticity & Crosslinking (Target)", "Module 4: Terminal Fibrogenesis & Scarring")
)

n_points <- 100
pt_grid <- seq(min(pt_vec), max(pt_vec), length.out = n_points)
smoothed_matrix <- matrix(NA, nrow = length(final_genes), ncol = n_points, dimnames = list(final_genes, paste0("PT_", 1:n_points)))

for (i in seq_along(final_genes)) {
  y <- expr_raw[, final_genes[i]]
  fit <- loess(y ~ pt_vec, span = 0.5)
  pred <- predict(fit, newdata = data.frame(pt_vec = pt_grid))
  if (any(is.na(pred))) {
    pred <- zoo::na.locf(pred, na.rm = FALSE)
    pred <- zoo::na.locf(pred, fromLast = TRUE)
  }
  smoothed_matrix[i, ] <- pred
}

row_sds <- apply(smoothed_matrix, 1, sd)
valid_rows <- which(!is.na(row_sds) & row_sds > 1e-5)
smoothed_matrix <- smoothed_matrix[valid_rows, , drop = FALSE]
gene_modules <- gene_modules[valid_rows]

z_matrix <- t(scale(t(smoothed_matrix)))
z_matrix[is.na(z_matrix)] <- 0
z_matrix[z_matrix > 2] <- 2
z_matrix[z_matrix < -2] <- -2
rownames(z_matrix) <- toupper(rownames(z_matrix))
rownames(z_matrix)[rownames(z_matrix) %in% c("CTGF", "CCN2")] <- "CCN2/CTGF"

phase_annot <- cut(pt_grid, breaks = c(-Inf, 4.5, 11.5, Inf), labels = c("Phase I: Homeostasis", "Phase II: Plasticity Window", "Phase III: Fibrogenesis"))
top_anno <- HeatmapAnnotation(
  Phase = phase_annot,
  Pseudotime = pt_grid,
  col = list(
    Phase = c("Phase I: Homeostasis" = "#74C0FC", "Phase II: Plasticity Window" = "#FFA94D", "Phase III: Fibrogenesis" = "#FF6B6B"),
    Pseudotime = colorRamp2(seq(min(pt_grid), max(pt_grid), length.out = 3), c("#74C0FC", "#FFA94D", "#FF6B6B"))
  ),
  annotation_height = unit(c(5.5, 5.5), "mm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8.5, fontface = "bold", col = "#333333")
)

x1_ratio <- which.min(abs(pt_grid - 4.5)) / n_points
x2_ratio <- which.min(abs(pt_grid - 11.5)) / n_points

ht_full <- Heatmap(
  z_matrix, name = "Z-Score",
  col = colorRamp2(c(-2, 0, 2), c("#5C93C4", "#FFFFFF", "#D36387")),
  row_split = gene_modules, cluster_rows = TRUE, cluster_row_slices = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 7.8, fontface = "bold.italic", col = "#1A1A1A"),
  top_annotation = top_anno, row_title_gp = gpar(fontsize = 9.2, fontface = "bold", col = "#222222"),
  row_title_rot = 0, border = gpar(col = "#444444", lwd = 0.7)
)

draw_and_decorate_all <- function() {
  draw(ht_full, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = "Dynamic Transcriptional Cascades Coupling Iron Overload & Fibrogenesis",
       column_title_gp = gpar(fontsize = 12, fontface = "bold", col = "#111111"),
       padding = unit(c(5, 5, 5, 10), "mm"))
  for (slice_i in 1:4) {
    decorate_heatmap_body("Z-Score", slice = slice_i, {
      grid.lines(c(x1_ratio, x1_ratio), c(0, 1), gp = gpar(lty = "dashed", col = "#D36387", lwd = 1.1))
      grid.lines(c(x2_ratio, x2_ratio), c(0, 1), gp = gpar(lty = "dashed", col = "#D36387", lwd = 1.1))
    })
  }
}

pdf(file.path(save_dir, "Fig2D_Pseudotime_Heatmap_AllGenes_Labeled.pdf"), width = 9.5, height = 9.2)
draw_and_decorate_all()
dev.off()

png(file.path(save_dir, "Fig2D_Pseudotime_Heatmap_AllGenes_Labeled.png"), width = 9.5, height = 9.2, units = "in", res = 300)
draw_and_decorate_all()
dev.off()