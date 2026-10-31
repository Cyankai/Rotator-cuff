```r
#!/usr/bin/env Rscript
# ==============================================================================
# 03_Fig2_DEG_Enrichment_and_Pathways.R
# 功能：生成空间基因矩阵、GraphPad 风格火山图、GO GSEA 背靠背三图、KEGG 与 11 亚群通路热图
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggrepel)
  library(ggnewscale)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
})

# 1. 路径配置与数据加载
base_dir <- "/Users/kenny/Desktop/凡/空转/人标本/q_d_2/analysis"
data_dir <- file.path(base_dir, "outs")
save_dir <- file.path(base_dir, "outs")
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

seu_anno <- readRDS(file.path(base_dir, "seu_anno_8um.rds"))
clean_df <- readRDS(file.path(base_dir, "clean_df_coords.rds"))

# ==============================================================================
# 模块 A：Fig 2 Top 3行 × 7列 空间基因表达矩阵 (Healthy / Fibrosis / Iron)
# ==============================================================================
healthy_genes  <- c("COL1A1", "COL1A2", "DCN", "SCX", "COMP", "FMOD", "TNC")
fibrosis_genes <- c("COL3A1", "ACTA2", "POSTN", "MMP2", "LOX", "MMP9", "FN1")
iron_genes     <- c("GPX4", "FTL", "FTH1", "NCOA4", "TFRC", "SLC40A1", "HMOX1")

plot_spatial_gene_row <- function(genes, high_col) {
  plist <- lapply(genes, function(g) {
    if (!g %in% rownames(seu_anno)) return(ggplot() + ggtitle(g) + theme_void())
    df_g <- clean_df
    df_g$Expr <- FetchData(seu_anno, vars = g)[df_g$Barcode, 1]
    df_g <- df_g %>% arrange(Expr)
    ggplot(df_g, aes(x = x, y = -y, color = Expr)) +
      geom_point(size = 0.35, stroke = 0) +
      scale_color_gradient(low = "#EFEFEF", high = high_col) +
      coord_fixed(ratio = 1) +
      ggtitle(g) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 9),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width  = unit(0.12, "cm"),
        legend.title = element_blank()
      )
  })
  wrap_plots(plist, nrow = 1)
}

p_grid_3x7 <- plot_spatial_gene_row(healthy_genes, "#08519C") /
  plot_spatial_gene_row(fibrosis_genes, "#CB181D") /
  plot_spatial_gene_row(iron_genes, "#238B45")

ggsave(file.path(save_dir, "Figure_Spatial_3x7_Gene_Grid.png"), p_grid_3x7, width = 20, height = 8.5, dpi = 350)

# ==============================================================================
# 模块 B：计算或加载差异表达结果 (deg_region & deg_df) 并绘制 GraphPad 火山图
# ==============================================================================
deg_csv <- file.path(save_dir, "deg_results.csv")

if (file.exists(deg_csv)) {
  deg_df <- read.csv(deg_csv, stringsAsFactors = FALSE)
  if (!"avg_log2FC" %in% colnames(deg_df)) deg_df$avg_log2FC <- deg_df$log2FC
  rownames(deg_df) <- deg_df$gene
  deg_region <- deg_df
} else {
  message(">>> 正在计算 Fibrotic vs Healthy 区域差异表达基因...")
  Idents(seu_anno) <- "Region"
  deg_region <- FindMarkers(seu_anno, ident.1 = "Fibrotic", ident.2 = "Healthy", logfc.threshold = 0.05, min.pct = 0.02)
  deg_region$gene <- rownames(deg_region)
  deg_df <- deg_region %>%
    transmute(
      gene = gene,
      log2FC = avg_log2FC,
      avg_log2FC = avg_log2FC,
      padj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)
    )
  write.csv(deg_df, deg_csv, row.names = FALSE)
}

fc_cutoff <- 0.5
pval_log_cutoff <- 50

deg_df <- deg_df %>%
  mutate(
    padj = ifelse(padj == 0, 1e-300, padj),
    neg_log10_pval = -log10(padj),
    change = case_when(
      log2FC > fc_cutoff & neg_log10_pval > pval_log_cutoff ~ "Upregulated",
      log2FC < -fc_cutoff & neg_log10_pval > pval_log_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

target_labels <- c("COL1A1", "COL1A2", "POSTN", "COL3A1", "TFRC", "LOXL2", "HMOX1", "CHI3L2", "SOD2", "CLU", "CRTAC1", "CRTACS1", "MMP3", "DCN")
deg_df <- deg_df %>% mutate(is_label = gene %in% target_labels)

p_volcano <- ggplot(deg_df, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = change), size = 1.8, alpha = 0.75, stroke = 0) +
  scale_color_manual(values = c("Upregulated" = "#E26D8C", "Downregulated" = "#74A9CF", "Not Significant" = "#D3D3D3")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted", color = "#333333", linewidth = 0.7) +
  geom_text_repel(
    data = subset(deg_df, is_label), aes(label = gene),
    size = 3.5, fontface = "italic", color = "black",
    box.padding = 0.4, point.padding = 0.25, segment.color = "#666666", segment.linewidth = 0.4, show.legend = FALSE
  ) +
  labs(x = expression(bold(log[2]~Fold~Change)), y = expression(bold(-log[10]~Adjusted~italic(P)~value))) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1.1),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", face = "bold", size = 11),
    axis.title = element_text(color = "black", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(file.path(save_dir, "Figure_Volcano_GraphPad_Style.pdf"), plot = p_volcano, width = 7.0, height = 5.5, device = cairo_pdf)
ggsave(file.path(save_dir, "Figure_Volcano_GraphPad_Style.png"), plot = p_volcano, width = 7.0, height = 5.5, dpi = 400)

# ==============================================================================
# 模块 C：GO GSEA 背靠背柱状图 (BP, CC, MF 三图) 与 KEGG 散点图
# ==============================================================================
deg_clean <- deg_region %>%
  filter(!is.na(avg_log2FC) & !is.na(gene)) %>%
  arrange(desc(avg_log2FC))

gene_ids <- bitr(deg_clean$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
deg_mapped <- deg_clean %>% inner_join(gene_ids, by = c("gene" = "SYMBOL")) %>% distinct(ENTREZID, .keep_all = TRUE)

gsea_input <- deg_mapped$avg_log2FC
names(gsea_input) <- deg_mapped$ENTREZID
gsea_input <- sort(gsea_input, decreasing = TRUE)

plot_gsea_bimodal <- function(gsea_df, ont_name, title_label) {
  if (nrow(gsea_df) == 0) return(NULL)
  top_fibrosis <- gsea_df %>% filter(NES > 0) %>% slice_max(order_by = NES, n = 7) %>% mutate(Dominant = "Fibrosis", Score_Fibrosis = NES, Score_Healthy = 0)
  top_healthy  <- gsea_df %>% filter(NES < 0) %>% distinct(Description, .keep_all = TRUE) %>% slice_min(order_by = NES, n = 6) %>% mutate(Dominant = "Healthy", Score_Fibrosis = 0, Score_Healthy = abs(NES))
  
  plot_df <- bind_rows(top_fibrosis, top_healthy) %>%
    mutate(Term_clean = str_to_lower(Description), Term_wrap = str_wrap(Term_clean, width = 34)) %>%
    arrange(factor(Dominant, levels = c("Healthy", "Fibrosis")), ifelse(Dominant == "Healthy", -Score_Healthy, Score_Fibrosis)) %>%
    mutate(y_id = row_number())
  
  max_score <- max(c(plot_df$Score_Fibrosis, plot_df$Score_Healthy), na.rm = TRUE)
  axis_lim <- ceiling(max_score * 1.15)
  breaks_seq <- seq(0, axis_lim, by = 1)
  all_breaks <- c(-rev(breaks_seq[-1]), breaks_seq)
  all_labels <- as.character(abs(all_breaks))
  
  ggplot(plot_df) +
    geom_rect(data = subset(plot_df, Dominant == "Healthy"), aes(xmin = -Score_Healthy, xmax = 0, ymin = y_id - 0.38, ymax = y_id + 0.38, fill = Score_Healthy), color = "black", linewidth = 0.35) +
    scale_fill_gradient(low = "#BDD7E7", high = "#08519C", guide = "none") +
    new_scale_fill() +
    geom_rect(data = subset(plot_df, Dominant == "Fibrosis"), aes(xmin = 0, xmax = Score_Fibrosis, ymin = y_id - 0.38, ymax = y_id + 0.38, fill = Score_Fibrosis), color = "black", linewidth = 0.35) +
    scale_fill_gradient(low = "#FCAE91", high = "#67000D", guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.65, color = "black") +
    geom_text(data = subset(plot_df, Dominant == "Fibrosis"), aes(x = -0.06, y = y_id, label = Term_wrap), hjust = 1, vjust = 0.5, size = 3.2, color = "#111111", lineheight = 0.85) +
    geom_text(data = subset(plot_df, Dominant == "Healthy"), aes(x = 0.06, y = y_id, label = Term_wrap), hjust = 0, vjust = 0.5, size = 3.2, color = "#111111", lineheight = 0.85) +
    scale_x_continuous(breaks = all_breaks, labels = all_labels, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-axis_lim * 1.25, axis_lim * 1.25), ylim = c(0.3, nrow(plot_df) + 0.8), clip = "off") +
    labs(x = NULL, y = NULL, title = title_label, caption = "Left: Enriched in Healthy Tendon   |   Right: Enriched in Fibrotic Region") +
    theme_classic(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 15)),
      plot.caption = element_text(size = 10, face = "bold.italic", colour = "grey30", hjust = 0.5, margin = margin(t = 12)),
      axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.55, colour = "black"),
      axis.ticks.x = element_line(linewidth = 0.55, colour = "black"),
      axis.text.x = element_text(size = 10, colour = "black", face = "bold"),
      plot.margin = margin(15, 20, 20, 20)
    )
}

ontologies <- c("BP", "CC", "MF")
titles <- c("A     Biological Process Enrichment", "B     Cellular Component Enrichment", "C     Molecular Function Enrichment")
filenames <- c("Figure_Reference_Bimodal_GO_BP", "Figure_Reference_Bimodal_GO_CC", "Figure_Reference_Bimodal_GO_MF")

set.seed(1234)
for (i in seq_along(ontologies)) {
  ont <- ontologies[i]
  gsea_res <- gseGO(geneList = gsea_input, OrgDb = org.Hs.eg.db, ont = ont, minGSSize = 5, maxGSSize = 600, pvalueCutoff = 0.05, pAdjustMethod = "none", verbose = FALSE)
  gsea_df <- as.data.frame(gsea_res)
  write.csv(gsea_df, file = file.path(save_dir, paste0("GSEA_GO_", ont, "_All.csv")), row.names = FALSE)
  p_plot <- plot_gsea_bimodal(gsea_df, ont, titles[i])
  if (!is.null(p_plot)) {
    ggsave(file.path(save_dir, paste0(filenames[i], ".pdf")), plot = p_plot, width = 8.0, height = 6.4, device = "pdf")
    ggsave(file.path(save_dir, paste0(filenames[i], ".png")), plot = p_plot, width = 8.0, height = 6.4, dpi = 300)
  }
}

# ==============================================================================
# 模块 D：7 条核心通路 × 11 细胞亚群活性热图 (Figure_Pathway_Heatmap_11Clusters_Spaced)
# ==============================================================================
Idents(seu_anno) <- "CellType"
pathway_list <- list(
  "WNT / \beta-Catenin Signaling"     = c("CTNNB1", "WNT5A", "FZD1", "LEF1", "AXIN2", "TCF7L2"),
  "NF-kb Inflammation"            = c("NFKB1", "RELA", "IL6", "CCL2", "CXCL8", "TNFAIP3"),
  "Iron Overload & Ferroptosis"        = c("TFRC", "ACSL4", "SAT1", "NCOA4", "HMOX1", "PTGS2"),
  "Ferritin Iron Sequestration"        = c("FTL", "FTH1", "SLC40A1", "GPX4", "SOD2"),
  "TGF-beta Signaling"               = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1"),
  "Myofibroblast Activation"           = c("ACTA2", "TAGLN", "MYL9", "CTHRC1", "POSTN", "FAP"),
  "ECM Crosslinking & Remodeling"      = c("LOX", "LOXL2", "PLOD2", "COL3A1", "FN1", "MMP14")
)

seu_anno <- AddModuleScore(seu_anno, features = pathway_list, name = "PW_")
pw_cols <- paste0("PW_", seq_along(pathway_list))

mat_raw <- seu_anno@meta.data %>%
  group_by(CellType) %>%
  summarise(across(all_of(pw_cols), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("CellType") %>%
  as.matrix() %>%
  t()
rownames(mat_raw) <- names(pathway_list)

mat_zscore <- t(scale(t(mat_raw)))
target_celltypes <- colnames(mat_zscore)

celltype_palette <- c(
  "Active Fibroblasts"          = "#EE8A82",
  "ECM-Remodeling Fibroblasts"  = "#7FB2D9",
  "Endothelial Cells"           = "#F27E8A",
  "Ferritin-High Macrophages"   = "#A7CDA8",
  "Fibro-Immune Transitional"   = "#D99B82",
  "Myofibroblasts"              = "#D991C2",
  "Pro-fibrotic Fibroblasts"    = "#E8BA73",
  "Quiescent Fibroblasts-1"     = "#8287C6",
  "Quiescent Fibroblasts-2"     = "#57B8B2",
  "Tendon Core Fibroblasts"     = "#8EBE94",
  "TFRC+ Macrophages"           = "#F9E79F"
)

ha_col <- HeatmapAnnotation(
  Cell_Type = factor(target_celltypes, levels = target_celltypes),
  col = list(Cell_Type = celltype_palette[target_celltypes]),
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.32, "cm"),
  gp = gpar(col = "white", lwd = 1.2)
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#5C93C4", "#FFFFFF", "#D05C82"))

ht_opt <- Heatmap(
  mat_zscore,
  name = "Enrichment\nZ-score",
  col = col_fun,
  top_annotation = ha_col,
  cluster_columns = TRUE, cluster_rows = TRUE,
  clustering_distance_columns = "pearson", clustering_method_columns = "ward.D2",
  column_dend_height = unit(1.2, "cm"), row_dend_width = unit(1.2, "cm"),
  width = unit(ncol(mat_zscore) * 0.85, "cm"),
  height = unit(nrow(mat_zscore) * 0.55, "cm"),
  rect_gp = gpar(col = "#F1F3F5", lwd = 1.0),
  row_names_gp = gpar(fontsize = 9.5, fontface = "plain"),
  column_names_gp = gpar(fontsize = 8.5, fontface = "plain"),
  column_names_rot = 45, column_names_centered = FALSE
)

pdf(file.path(save_dir, "Figure_Pathway_Heatmap_11Clusters_Spaced.pdf"), width = 10.5, height = 6.2)
draw(ht_opt, padding = unit(c(8, 16, 22, 10), "mm"))
dev.off()

png(file.path(save_dir, "Figure_Pathway_Heatmap_11Clusters_Spaced.png"), width = 10.5, height = 6.2, units = "in", res = 400)
draw(ht_opt, padding = unit(c(8, 16, 22, 10), "mm"))
dev.off()