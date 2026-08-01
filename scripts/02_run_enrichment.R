suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(rrvgo)
  library(enrichplot)
  library(here)
})

a <- readRDS(here::here("objects", "analysis_results.rds"))
list2env(a, envir = environment())

directional_go <- function(de_df, tag, lfc_cut = 1, padj_cut = 0.01) {
  de_df$ensembl_base <- sub("\\..*$", "", de_df$gene_id)
  universe <- unique(de_df$ensembl_base[!is.na(de_df$padj)])
  up   <- unique(de_df$ensembl_base[!is.na(de_df$padj) & de_df$padj < padj_cut & de_df$log2FoldChange >  lfc_cut])
  down <- unique(de_df$ensembl_base[!is.na(de_df$padj) & de_df$padj < padj_cut & de_df$log2FoldChange < -lfc_cut])
  cat(tag, ":", length(up), "UP,", length(down), "DOWN (padj<", padj_cut, ", |LFC|>", lfc_cut, ")\n")

  ego_up   <- enrichGO(up,   universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                        ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego_down <- enrichGO(down, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                        ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  list(up = ego_up, down = ego_down)
}

run_gsea <- function(de_df, ranking_column) {
  de_df$ensembl_base <- sub("\\..*$", "", de_df$gene_id)
  df <- de_df %>% filter(!is.na(.data[[ranking_column]]), !is.na(ensembl_base)) %>%
    distinct(ensembl_base, .keep_all = TRUE)
  ranks <- setNames(df[[ranking_column]], df$ensembl_base)
  ranks <- sort(ranks, decreasing = TRUE)
  gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
        pvalueCutoff = 0.05, minGSSize = 15, maxGSSize = 500, eps = 0, seed = TRUE)
}

# Run enrichment once
go_KO <- directional_go(res_KO_df, "KO")
go_DEL <- directional_go(res_DEL_df, "DEL")
gsea_KO <- run_gsea(res_KO_df, "stat")
gsea_DEL <- run_gsea(res_DEL_df, "stat")

merge_go_gsea <- function(ora_up, ora_down, gsea_obj, tag) {
  ora_all <- bind_rows(
    as.data.frame(ora_up)   %>% mutate(direction = "UP"),
    as.data.frame(ora_down) %>% mutate(direction = "DOWN")
  ) %>% dplyr::select(ID, Description, direction, ora_padj = p.adjust)

  gsea_df <- as.data.frame(gsea_obj) %>%
    dplyr::select(ID, Description, NES, gsea_padj = p.adjust) %>%
    mutate(direction = ifelse(NES > 0, "UP", "DOWN"))

  merged <- full_join(ora_all, gsea_df, by = c("ID", "Description", "direction")) %>%
    mutate(in_both = !is.na(ora_padj) & !is.na(gsea_padj),
           condition = tag) %>%
    arrange(desc(in_both), ora_padj, gsea_padj)
  merged
}
merged_KO  <- merge_go_gsea(go_KO$up,  go_KO$down,  gsea_KO,  "KO")
merged_DEL <- merge_go_gsea(go_DEL$up, go_DEL$down, gsea_DEL, "DEL")
merged_all <- bind_rows(merged_KO, merged_DEL)

cat("KO: ", sum(merged_KO$in_both, na.rm = TRUE), "terms significant in both ORA and GSEA\n")
cat("DEL:", sum(merged_DEL$in_both, na.rm = TRUE), "terms significant in both ORA and GSEA\n")

# Semantic reduction
hsGO <- godata('org.Hs.eg.db', ont = "BP")

gsea_KO_simplified <- clusterProfiler::simplify(gsea_KO, cutoff = 0.7, by = "p.adjust",
                                                 select_fun = min, semData = hsGO)




go_ids <- gsea_KO@result$ID
simMatrix <- calculateSimMatrix(go_ids, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")

scores <- setNames(-log10(gsea_KO@result$p.adjust), gsea_KO@result$ID)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")

cat("Original terms:", nrow(simMatrix), "| Reduced to", length(unique(reducedTerms$parentTerm)), "parent clusters\n")



gsea_KO_pt <- pairwise_termsim(gsea_KO)


enrichment_results <- list(
  go_KO = go_KO, go_DEL = go_DEL, gsea_KO = gsea_KO, gsea_DEL = gsea_DEL,
  merged_KO = merged_KO, merged_DEL = merged_DEL, merged_all = merged_all,
  hsGO = hsGO, gsea_KO_simplified = gsea_KO_simplified,
  reducedTerms = reducedTerms, gsea_KO_pt = gsea_KO_pt,
  settings = list(go_ontology = "ALL for ORA; BP for GSEA",
                  gsea_rank = "Wald statistic", generated = Sys.time())
)
saveRDS(enrichment_results, here::here("objects", "enrichment_results.rds"))
writeLines(capture.output(sessionInfo()), here::here("results", "sessionInfo_enrichment.txt"))
cat("Saved objects/enrichment_results.rds\n")
