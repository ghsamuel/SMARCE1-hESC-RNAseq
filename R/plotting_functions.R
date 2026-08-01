# Simple plotting helpers used by the Quarto report

make_volcano <- function(res_lfc, res_raw, title) {
  df <- as.data.frame(res_raw) %>% rownames_to_column("gene_id") %>%
    mutate(gene_symbol = sym[gene_id])
  top <- df %>% filter(!is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc_cutoff) %>%
    arrange(padj) %>% slice_head(n = 15)
  EnhancedVolcano(res_lfc, lab = sym[rownames(res_lfc)],
                   x = "log2FoldChange", y = "padj", title = title,
                   pCutoff = alpha, FCcutoff = lfc_cutoff,
                   pointSize = 2.0, labSize = 3.5,
                   selectLab = top$gene_symbol, drawConnectors = TRUE,
                   widthConnectors = 0.5, max.overlaps = 20)
}

make_ma <- function(res_lfc, title) {
  df <- as.data.frame(res_lfc) %>% rownames_to_column("gene_id") %>%
    mutate(sig = case_when(
      !is.na(padj) & padj < alpha & log2FoldChange >  lfc_cutoff ~ "Up",
      !is.na(padj) & padj < alpha & log2FoldChange < -lfc_cutoff ~ "Down",
      TRUE ~ "NS"
    ))
  ggplot(df, aes(baseMean, log2FoldChange, color = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_x_log10() +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_cartesian(ylim = c(-5, 5)) +
    theme_bw() +
    labs(title = title, x = "Mean of normalized counts", y = "log2 fold change", color = NULL)
}

plot_gene_by_clone <- function(gene_id_target, title_gene) {
  if (is.na(gene_id_target) || !(gene_id_target %in% rownames(norm_counts_all))) return(invisible(NULL))
  df <- data.frame(genotype = samples$genotype, count = norm_counts_all[gene_id_target, ])
  ggplot(df, aes(genotype, count, color = genotype)) +
    geom_jitter(width = 0.1, size = 2, show.legend = FALSE) +
    theme_bw(base_size = 8) +
    labs(y = NULL, x = NULL, title = title_gene) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

top_n_genes <- function(de_df, dir, n = 10) {
  d <- de_df %>% filter(!is.na(padj))
  if (dir == "up")   d <- d %>% filter(log2FoldChange > 0) %>% arrange(padj)
  if (dir == "down") d <- d %>% filter(log2FoldChange < 0) %>% arrange(padj)
  d %>% slice_head(n = n)
}

plot_top_panel <- function(de_df, tag, dir, n = 10) {
  genes <- top_n_genes(de_df, dir, n)
  plots <- map2(genes$gene_id, genes$gene_symbol, ~plot_gene_by_clone(.x, .y))
  plots <- plots[!map_lgl(plots, is.null)]
  wrap_plots(plots, ncol = 5) + plot_annotation(title = paste0(tag, " top ", n, " ", toupper(dir)))
}
