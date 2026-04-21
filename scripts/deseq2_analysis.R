#!~/.conda/envs/snakemake/bin/R
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
padj_cutoff          <- snakemake@params[["padj_cutoff"]]
lfc_cutoff           <- snakemake@params[["lfc_cutoff"]]
count_type           <- snakemake@params[["count_type"]]
use_tissue_covariate <- snakemake@params[["use_tissue_covariate"]]
split_by_tissue      <- snakemake@params[["split_by_tissue"]]

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

# 检查是否有 tissue 列
has_tissue <- "tissue" %in% colnames(metadata)
if (has_tissue) {
  metadata$tissue <- factor(metadata$tissue)
  cat(sprintf("📋 读取 %d 个样本的 metadata（包含 %d 个组织）\n",
              nrow(metadata), length(unique(metadata$tissue))))
} else {
  cat(sprintf("📋 读取 %d 个样本的 metadata\n", nrow(metadata)))
}
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
# 有 tissue 列且开启协变量时，用 ~ tissue + condition 控制组织间差异
if (has_tissue && isTRUE(use_tissue_covariate)) {
  cat("📐 Design formula: ~ tissue + condition\n")
  design_formula <- ~ tissue + condition
} else {
  cat("📐 Design formula: ~ condition\n")
  design_formula <- ~ condition
}

dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = design_formula)

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

intgroup <- if (has_tissue) c("condition", "tissue") else "condition"
pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

if (has_tissue) {
  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = condition, shape = tissue, label = name)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_shape_manual(values = seq(15, 15 + length(unique(pca_data$tissue)))) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top") +
    ggtitle("PCA Plot - Sample Clustering (color=condition, shape=tissue)")
} else {
  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top") +
    ggtitle("PCA Plot - Sample Clustering")
}

ggsave(file.path(outdir, "pca_plot.pdf"), p_pca, width = 8, height = 6)
cat("✅ PCA 图已保存至 results/deseq2/pca_plot.pdf\n")

# ── 7. 差异基因分析函数 ───────────────────────────────────
run_deseq2_contrast <- function(dds, vsd, metadata, contrast_pair,
                                 outdir, padj_cutoff, lfc_cutoff, has_tissue) {
  treatment     <- contrast_pair[1]
  control       <- contrast_pair[2]
  contrast_name <- paste0(treatment, "_vs_", control)

  cat(sprintf("\n🔬 分析比较组: %s vs %s\n", treatment, control))

  res <- results(dds, contrast = c("condition", treatment, control),
                 alpha = padj_cutoff)
  res <- res[order(res$padj), ]

  res$diffexpressed <- "NO"
  res$diffexpressed[res$log2FoldChange >  lfc_cutoff & res$padj < padj_cutoff] <- "UP"
  res$diffexpressed[res$log2FoldChange < -lfc_cutoff & res$padj < padj_cutoff] <- "DOWN"

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                       "stat", "pvalue", "padj", "diffexpressed")]

  output_file <- file.path(outdir, paste0(contrast_name, "_DEG_results.csv"))
  write.csv(res_df, output_file, row.names = FALSE, quote = FALSE)

  n_up   <- sum(res_df$diffexpressed == "UP",   na.rm = TRUE)
  n_down <- sum(res_df$diffexpressed == "DOWN",  na.rm = TRUE)
  cat(sprintf("  ⬆️  上调基因: %d\n", n_up))
  cat(sprintf("  ⬇️  下调基因: %d\n", n_down))
  cat(sprintf("  💾 结果已保存至: %s\n", output_file))

  # Volcano plot
  res_plot <- res_df[!is.na(res_df$padj), ] %>%
    mutate(log10padj = -log10(padj)) %>%
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
         x = "log2 Fold Change", y = "-log10(adjusted p-value)",
         color = "Differential Expression") +
    theme(legend.position = "top")

  ggsave(file.path(outdir, paste0(contrast_name, "_volcano.pdf")),
         p_volcano, width = 10, height = 8)

  # Heatmap（top 50 DEGs）
  if (n_up + n_down >= 10) {
    top_genes <- res_df %>%
      filter(diffexpressed != "NO") %>%
      arrange(padj) %>%
      head(50) %>%
      pull(gene_id)

    mat <- assay(vsd)[top_genes, rownames(metadata)]
    mat <- t(scale(t(mat)))

    # Heatmap 注释：有 tissue 列时同时显示 condition 和 tissue
    if (has_tissue) {
      annotation_col <- data.frame(
        Condition = metadata$condition,
        Tissue    = metadata$tissue,
        row.names = rownames(metadata)
      )
    } else {
      annotation_col <- data.frame(
        Condition = metadata$condition,
        row.names = rownames(metadata)
      )
    }

    heatmap_file <- file.path(outdir, paste0(contrast_name, "_heatmap.pdf"))
    pdf(heatmap_file, width = 10, height = 12)
    pheatmap(mat,
             annotation_col = annotation_col,
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 8,
             main = paste0("Top 50 DEGs: ", contrast_name))
    dev.off()
    cat(sprintf("  🔥 Heatmap 已保存至: %s\n", heatmap_file))
  }
}

# ── 8. 执行差异分析 ───────────────────────────────────────
# 模式一：split_by_tissue = true，每个组织单独跑一套 DESeq2
if (has_tissue && isTRUE(split_by_tissue)) {
  tissues <- levels(metadata$tissue)
  cat(sprintf("\n🔀 按组织分组分析，共 %d 个组织: %s\n",
              length(tissues), paste(tissues, collapse = ", ")))

  for (tis in tissues) {
    cat(sprintf("\n━━━ 组织: %s ━━━\n", tis))
    tis_samples  <- rownames(metadata)[metadata$tissue == tis]
    tis_meta     <- metadata[tis_samples, , drop = FALSE]
    tis_counts   <- txi$counts[, tis_samples]
    tis_txi      <- txi
    tis_txi$counts    <- txi$counts[, tis_samples]
    tis_txi$abundance <- txi$abundance[, tis_samples]
    tis_txi$length    <- txi$length[, tis_samples]

    dds_tis <- DESeqDataSetFromTximport(tis_txi, colData = tis_meta,
                                        design = ~ condition)
    keep_tis <- rowSums(counts(dds_tis) >= 10) >= max(2, floor(nrow(tis_meta) * 0.5))
    dds_tis  <- dds_tis[keep_tis, ]
    dds_tis  <- DESeq(dds_tis)
    vsd_tis  <- vst(dds_tis, blind = FALSE)

    tis_outdir <- file.path(outdir, tis)
    dir.create(tis_outdir, showWarnings = FALSE, recursive = TRUE)

    for (contrast_pair in contrasts) {
      run_deseq2_contrast(dds_tis, vsd_tis, tis_meta, contrast_pair,
                          tis_outdir, padj_cutoff, lfc_cutoff, has_tissue = FALSE)
    }
  }

# 模式二：所有样本一起跑，tissue 作为协变量（默认）
} else {
  for (contrast_pair in contrasts) {
    run_deseq2_contrast(dds, vsd, metadata, contrast_pair,
                        outdir, padj_cutoff, lfc_cutoff, has_tissue)
  }
}

cat("\n✅ DESeq2 分析完成！\n")
