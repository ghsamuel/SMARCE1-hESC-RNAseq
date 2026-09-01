# One-off: relabel gene-member columns in the saved enrichment objects from
# ENSEMBL to SYMBOL. Seconds, not a recompute. Safe to re-run (setReadable is
# a no-op on an already-readable object).
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(here)
})

e <- readRDS(here::here("objects", "enrichment_results.rds"))

mk <- function(x) {
  if (is.null(x)) return(x)
  tryCatch(setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL"),
           error = function(err) { message("skipped: ", conditionMessage(err)); x })
}

e$gsea_KO            <- mk(e$gsea_KO)
e$gsea_DEL           <- mk(e$gsea_DEL)
e$go_KO$up           <- mk(e$go_KO$up)
e$go_KO$down         <- mk(e$go_KO$down)
e$go_DEL$up          <- mk(e$go_DEL$up)
e$go_DEL$down        <- mk(e$go_DEL$down)
e$gsea_KO_simplified <- mk(e$gsea_KO_simplified)
e$gsea_KO_pt         <- mk(e$gsea_KO_pt)

cat("core_enrichment now reads:\n  ",
    substr(as.data.frame(e$gsea_KO)$core_enrichment[1], 1, 90), "\n")

saveRDS(e, here::here("objects", "enrichment_results.rds"))
cat("Saved objects/enrichment_results.rds\n")
