suppressPackageStartupMessages(library(survival))

# Independent R rerun for the 1,200 retained validation datasets. The source
# CSVs contain synthetic records only. Results are written into the current
# audit folder so the archived validation files remain immutable.
args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
source_dir <- file.path(repo_root, "data", "large", "cox_validation_source")
output_dir <- file.path(repo_root, "verification", "cox_reverification")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- vector("list", 24L * 50L)
position <- 1L
for (combo in 0:23) {
  data <- read.csv(sprintf("%s/data_combo_%02d.csv", source_dir, combo))
  for (dataset_id in sort(unique(data$dataset_id))) {
    subset_data <- data[data$dataset_id == dataset_id, , drop = FALSE]
    fit <- coxph(
      Surv(time, event_c) ~ C,
      data = subset_data,
      weights = w,
      robust = TRUE,
      ties = "efron"
    )
    results[[position]] <- data.frame(
      combo = combo,
      dataset_id = dataset_id,
      beta_r = unname(coef(fit)[["C"]]),
      se_r = sqrt(fit$var[1, 1])
    )
    position <- position + 1L
  }
  message(sprintf("R coxph validation combo %02d/23 complete", combo))
}

result <- do.call(rbind, results)
stopifnot(nrow(result) == 1200L, all(is.finite(result$beta_r)), all(is.finite(result$se_r)))
write.csv(result, file.path(output_dir, "R_coxph_estimates_rerun.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
message("R coxph validation finished: 1,200 fits")
