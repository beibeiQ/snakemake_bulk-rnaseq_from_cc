#!/usr/bin/env Rscript
# ============================================================
# DESeq2 差异基因分析脚本
# 输入: Salmon quant.sf 文件 + metadata
# 输出: PCA 图 + 差异基因表
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ── 参数获取 ──────────────────────────────────────────────
salmon_dir   <- snakemake@params[["salmon_dir"]]
metadata_file<- snakemake@input[["metadata"]]
gtf_file     <- snakemake@input[["gtf"]]
outdir       <- snakemake@params[["outdir"]]
contrasts    <- snakemake@params[["contrasts"]]
padj_cutoff  <- snakemake@params[["padj_cutoff"]]
lfc_cutoff   <- snakemake@params[["lfc_cutoff"]]
count_type   <- snakemake@params[["count_type"]]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── 1. 构建 tx2gene 映射（从 GTF 提取）──────────────────
cat("📖 解析 GTF 文件构建 transcript-to-gene 映射...\n")
gtf <- read.table(gtf_file, sep = "\t", comment.char = "#",
                  stringsAsFactors = FALSE, quote = "")
gtf_tx <- gtf[gtf$V3 == "transcript", ]

extract_attr <- function(attr_str, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  match <- str_extract(attr_str, pattern)
  if (is.na(match)) return(NA)
  str_replace(match, paste0(key, ' "'), "") %>% str_replace('"', "")
}

tx2gene <- data.frame(
  TXNAME = sapply(gtf_tx$V9, extract_attr, key = "transcript_id"),
  GENEID = sapply(gtf_tx$V9, extract_attr, key = "gene_id"),
  stringsAsFactors = FALSE
)
tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), ]
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), ]

cat(sprintf("✅ 成功提取 %d 个 transcript-to-gene 映射\n", nrow(tx2gene)))

# ── 2. 读取 metadata ──────────────────────────────────────
metadata <- read.table(metadata_file, header = TRUE, sep = "\t",
                       row.names = 1, stringsAsFactors = FALSE)
metadata$condition <- factor(metadata$condition)

cat(sprintf("📋 读取 %d 个样本的 metadata\n", nrow(metadata)))
print(metadata)

# ── 3. tximport 导入 Salmon 定量结果 ──────────────────────
samples <- rownames(metadata)
files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("❌ 以下 Salmon 输出文件缺失:\n", paste(missing, collapse = "\n"))
}

cat("📥 使用 tximport 导入 Salmon 定量数据...\n")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = count_type)

cat(sprintf("✅ 成功导入 %d 个基因的表达数据\n", nrow(txi$counts)))

# ── 4. 构建 DESeq2 对象 ───────────────────────────────────
dds <- DESeqDataSetFromTximport(txi, colData = metadata,
                                 design = ~ condition)

# 过滤低表达基因（至少 3 个样本中 count >= 10）
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat(sprintf("🔍 过滤后保留 %d 个基因\n", nrow(dds)))

# ── 5. 运行 DESeq2 ────────────────────────────────────────
cat("🧬 运行 DESeq2 差异分析...\n")
dds <- DESeq(dds)

# ── 6. PCA 图 ─────────────────────────────────────────────
cat("📊 生成 PCA 图...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  ggtitle("PCA Plot - Sample Clustering")

ggsave(file.path(outdir, "pca_plot.pdf"), p_pca, width = 8, height = 6)
cat("✅ PCA 图已保存至 results/deseq2/pca_plot.pdf\n")

# ── 7. 差异基因分析（支持多个比较组）──────────────────────
for (contrast_pair in contrasts) {
  treatment <- contrast_pair[1]
  control   <- contrast_pair[2]
  contrast_name <- paste0(treatment, "_vs_", control)

  cat(sprintf("\n🔬 分析比较组: %s vs %s\n", treatment, control))

  res <- results(dds, contrast = c("condition", treatment, control),
                 alpha = padj_cutoff)
  res <- res[order(res$padj), ]

  # 添加差异标签
  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange > lfc_cutoff & res$padj < padj_cutoff] <- "UP"
  res$diffexpressed[res$log2FoldChange < -lfc_cutoff & res$padj < padj_cutoff] <- "DOWN"

  # 保存完整结果
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj", "diffexpressed")]

  output_file <- file.path(outdir, paste0(contrast_name, "_DEG_results.csv"))
  write.csv(res_df, output_file, row.names = FALSE, quote = FALSE)

  # 统计
  n_up   <- sum(res$diffexpressed == "UP", na.rm = TRUE)
  n_down <- sum(res$diffexpressed == "DOWN", na.rm = TRUE)
  cat(sprintf("  ⬆️  上调基因: %d\n", n_up))
  cat(sprintf("  ⬇️  下调基因: %d\n", n_down))
  cat(sprintf("  💾 结果已保存至: %s\n", output_file))

  # Volcano plot
  res_plot <- res_df[!is.na(res_df$padj), ]
  res_plot$log10padj <- -log10(res_plot$padj)

  # 标记 top 基因
  res_plot <- res_plot %>%
    arrange(padj) %>%
    mutate(label = ifelse(row_number() <= 10 & diffexpressed != "NO", gene_id, ""))

  p_volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = log10padj,
                                     color = diffexpressed, label = label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(size = 3, max.overlaps = 15, color = "black") +
    scale_color_manual(values = c("UP" = "#E64B35", "DOWN" = "#4DBBD5", "NO" = "grey")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 14) +
    labs(title = paste0("Volcano Plot: ", contrast_name),
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)",
         color = "Differential Expression") +
    theme(legend.position = "top")

  volcano_file <- file.path(outdir, paste0(contrast_name, "_volcano.pdf"))
  ggsave(volcano_file, p_volcano, width = 10, height = 8)
  cat(sprintf("  📈 Volcano 图已保存至: %s\n", volcano_file))

  # Heatmap（top 50 DEGs）
  if (n_up + n_down >= 10) {
    top_genes <- res_df %>%
      filter(diffexpressed != "NO") %>%
      arrange(padj) %>%
      head(50) %>%
      pull(gene_id)

    mat <- assay(vsd)[top_genes, ]
    mat <- t(scale(t(mat)))  # Z-score 标准化

    annotation_col <- data.frame(
      Condition = metadata$condition,
      row.names = rownames(metadata)
    )

    heatmap_file <- file.path(outdir, paste0(contrast_name, "_heatmap.pdf"))
    pdf(heatmap_file, width = 10, height = 12)
    pheatmap(mat,
             annotation_col = annotation_col,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 8,
             main = paste0("Top 50 DEGs: ", contrast_name))
    dev.off()
    cat(sprintf("  🔥 Heatmap 已保存至: %s\n", heatmap_file))
  }
}

cat("\n✅ DESeq2 分析完成！\n")
