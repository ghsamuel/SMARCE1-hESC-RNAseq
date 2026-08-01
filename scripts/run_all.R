# Run this only when the data, model, thresholds, or enrichment analysis changes.
source("scripts/01_run_analysis.R")
source("scripts/02_run_enrichment.R")
cat("Analysis complete. You can now render hesc_KO_DEL.qmd.\n")
