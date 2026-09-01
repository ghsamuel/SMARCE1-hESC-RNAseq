suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(tximport)
  library(here)
})

dir.create(here::here("objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results", "tables"), recursive = TRUE, showWarnings = FALSE)

# Import counts
samples <- read.csv(here::here("data", "metadata.csv"), header = TRUE, stringsAsFactors = FALSE)

files <- file.path(here::here("data", "quant_trimmed"), samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

tx2gene <- read.delim(here::here("data", "tx2gene.tsv"), header = FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID", "GENENAME")

txi <- tximport(files, type = "salmon",
                tx2gene = tx2gene[, c("TXNAME", "GENEID")],
                countsFromAbundance = "no")

sym <- setNames(tx2gene$GENENAME, tx2gene$GENEID)

# Check sample order before analysis
stopifnot(all(samples$sample_id %in% names(files)))
stopifnot(identical(samples$sample_id, names(files)))

# Filter and run clone-level DESeq2
samples$family <- case_when(
  samples$genotype == "WT" ~ "WT",
  samples$genotype %in% c("null1", "null2") ~ "KO",
  samples$genotype %in% c("del1", "del2")   ~ "DEL"
)

# Set factors explicitly so DESeq2 does not convert character columns silently
samples$genotype <- relevel(factor(samples$genotype), ref = "WT")
samples$family <- factor(samples$family, levels = c("WT", "KO", "DEL"))

# Confirm imported count columns and metadata are in the same order
stopifnot(identical(samples$sample_id, colnames(txi$counts)))

ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ genotype)

sumcounts <- rowSums(counts(ddsTxi))

keep <- sumcounts > 50
ddsTxi <- ddsTxi[keep, ]
cat("Genes retained after filtering:", nrow(ddsTxi), "\n")

vsd <- vst(ddsTxi, blind = FALSE)
pca_data <- plotPCA(vsd, returnData = TRUE, intgroup = "genotype", ntop = 500)
percentVar <- round(100 * attr(pca_data, "percentVar"))


dds <- DESeq(ddsTxi)
resultsNames(dds)

alpha <- 0.01
lfc_cutoff <- 1

res_KO56 <- results(dds, contrast = c("genotype", "null1", "WT"), alpha = alpha)
res_KO67 <- results(dds, contrast = c("genotype", "null2", "WT"), alpha = alpha)
res_SL25 <- results(dds, contrast = c("genotype", "del1", "WT"), alpha = alpha)
res_SL45 <- results(dds, contrast = c("genotype", "del2", "WT"), alpha = alpha)

res_KO56_df <- as.data.frame(res_KO56) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))
res_KO67_df <- as.data.frame(res_KO67) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))
res_SL25_df <- as.data.frame(res_SL25) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))
res_SL45_df <- as.data.frame(res_SL45) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))

res_KO56_LFC <- lfcShrink(dds, coef = "genotype_null1_vs_WT", type = "apeglm")
res_KO67_LFC <- lfcShrink(dds, coef = "genotype_null2_vs_WT", type = "apeglm")
res_SL25_LFC <- lfcShrink(dds, coef = "genotype_del1_vs_WT", type = "apeglm")
res_SL45_LFC <- lfcShrink(dds, coef = "genotype_del2_vs_WT", type = "apeglm")

de_summary <- data.frame(
  Comparison = c("KO56 vs WT", "KO67 vs WT", "SL25 vs WT", "SL45 vs WT"),
  Total_DE = c(sum(res_KO56$padj < alpha, na.rm = TRUE),
               sum(res_KO67$padj < alpha, na.rm = TRUE),
               sum(res_SL25$padj < alpha, na.rm = TRUE),
               sum(res_SL45$padj < alpha, na.rm = TRUE))
)


# Family-level model
dds_family <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ family)
dds_family <- dds_family[rowSums(counts(dds_family)) > 50, ]
dds_family <- DESeq(dds_family)

res_KO  <- results(dds_family, contrast = c("family", "KO",  "WT"), alpha = alpha)
res_DEL <- results(dds_family, contrast = c("family", "DEL", "WT"), alpha = alpha)

res_KO_df  <- as.data.frame(res_KO)  %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))
res_DEL_df <- as.data.frame(res_DEL) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))

cat("KO family DE genes (padj<0.01):", sum(res_KO_df$padj < alpha, na.rm = TRUE), "\n")
cat("DEL family DE genes (padj<0.01):", sum(res_DEL_df$padj < alpha, na.rm = TRUE), "\n")

# Direct KO versus DEL model
ko_del_samples <- samples$sample_id[samples$family %in% c("KO", "DEL")]
cd_kd <- samples[samples$sample_id %in% ko_del_samples, ]
rownames(cd_kd) <- cd_kd$sample_id
cd_kd$family <- factor(cd_kd$family, levels = c("DEL", "KO"))   # DEL = reference

txi_kd <- tximport(files[ko_del_samples], type = "salmon",
                    tx2gene = tx2gene[, c("TXNAME", "GENEID")], countsFromAbundance = "no")
dds_kd <- DESeqDataSetFromTximport(txi_kd, colData = cd_kd, design = ~ family)
dds_kd <- dds_kd[rowSums(counts(dds_kd)) > 50, ]
dds_kd <- DESeq(dds_kd)

res_KO_vs_DEL_lfc <- lfcShrink(dds_kd, coef = "family_KO_vs_DEL", type = "apeglm")
res_KO_vs_DEL_raw <- results(dds_kd, name = "family_KO_vs_DEL")
res_KO_vs_DEL_df <- as.data.frame(res_KO_vs_DEL_lfc) %>% rownames_to_column("gene_id") %>%
  mutate(gene_symbol = coalesce(na_if(unname(sym[gene_id]), ""), gene_id))

cat("KO vs DEL: ", sum(res_KO_vs_DEL_df$padj < alpha, na.rm = TRUE), "genes differ directly (padj<", alpha, ")\n")


# Save the main analysis objects
analysis_results <- list(
  samples = samples, files = files, tx2gene = tx2gene, sym = sym, txi = txi,
  ddsTxi = ddsTxi, dds = dds, dds_family = dds_family, dds_kd = dds_kd,
  vsd = vsd, pca_data = pca_data, percentVar = percentVar, sumcounts = sumcounts,
  alpha = alpha, lfc_cutoff = lfc_cutoff, de_summary = de_summary,
  res_KO56 = res_KO56, res_KO67 = res_KO67, res_SL25 = res_SL25, res_SL45 = res_SL45,
  res_KO56_df = res_KO56_df, res_KO67_df = res_KO67_df,
  res_SL25_df = res_SL25_df, res_SL45_df = res_SL45_df,
  res_KO56_LFC = res_KO56_LFC, res_KO67_LFC = res_KO67_LFC,
  res_SL25_LFC = res_SL25_LFC, res_SL45_LFC = res_SL45_LFC,
  res_KO = res_KO, res_DEL = res_DEL, res_KO_df = res_KO_df, res_DEL_df = res_DEL_df,
  res_KO_vs_DEL_lfc = res_KO_vs_DEL_lfc, res_KO_vs_DEL_raw = res_KO_vs_DEL_raw,
  res_KO_vs_DEL_df = res_KO_vs_DEL_df,
  settings = list(alpha = alpha, lfc_cutoff = lfc_cutoff,
                  filtering = "rowSums(counts) > 50",
                  clone_design = "~ genotype", family_design = "~ family",
                  generated = Sys.time())
)
saveRDS(analysis_results, here::here("objects", "analysis_results.rds"))
writeLines(capture.output(sessionInfo()), here::here("results", "sessionInfo_analysis.txt"))
cat("Saved objects/analysis_results.rds\n")
