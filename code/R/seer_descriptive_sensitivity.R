suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

# This script deliberately reuses the independently audited cohort-reconstruction
# function, but replaces its final cleanup/return expressions in memory so the
# de-identified patient-level table can be used for additional sensitivity models.
# No patient-level data are written to disk.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript seer_descriptive_sensitivity.R <crc.csv> <bladder.csv> <output_dir>")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source_file <- file.path(script_dir, "seer_cohort_reconstruction.R")
source(source_file, local = .GlobalEnv)
reconstruct_keep <- reconstruct
old_body <- as.list(body(reconstruct_keep))
stopifnot(length(old_body) == 115L)
new_return <- quote(invisible(list(
  data = d, raw_n = raw_n, t1_n = t1_n,
  flow = as.data.table(do.call(rbind, steps)),
  outcomes = outcomes, ratio = ratio, overall = overall,
  ess = ess, unknown = unk,
  max_weighted_smd = max(bvals, na.rm = TRUE)
)))
body(reconstruct_keep) <- as.call(c(old_body[1:112], list(new_return)))

out_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_ps <- function(d) {
  factor_cols <- names(d)[vapply(d, is.factor, logical(1))]
  for (v in factor_cols) set(d, j = v, value = droplevels(d[[v]]))
  fit <- glm(
    exposure01 ~ ns(age, df = 4) + sex + race + marital + sitegrp +
      sizegrp + gradegrp + histgrp + ns(year, df = 4),
    data = d, family = binomial()
  )
  d[, ps := pmin(pmax(predict(fit, type = "response"), 1e-6), 1 - 1e-6)]
  d[, ow := ifelse(exposure01 == 1, 1 - ps, ps)]
  invisible(d)
}

endpoint_model <- function(d, event_var, endpoint, cohort, analysis) {
  ev <- as.integer(d$surv_days <= 5 * 365.25 & d[[event_var]] == 1)
  fit <- coxph(
    Surv(pmin(d$surv_days / 365.25, 5), ev) ~ exposure01,
    data = d, weights = ow, robust = TRUE, ties = "efron"
  )
  b <- unname(coef(fit)[1]); se <- sqrt(vcov(fit)[1, 1])
  ph <- tryCatch(cox.zph(fit, transform = "km")$table[1, "p"], error = function(e) NA_real_)
  data.table(
    cohort = cohort, analysis = analysis, endpoint = endpoint,
    n = nrow(d), events = sum(ev), log_hr = b, robust_se = se,
    hr = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
    ph_p = ph
  )
}

joint_model <- function(d, cohort, analysis) {
  z <- d[, .(
    pid = .I, exposure01, ow,
    t5 = pmin(surv_days / 365.25, 5),
    css = as.integer(surv_days <= 5 * 365.25 & css_event == 1),
    oth = as.integer(surv_days <= 5 * 365.25 & other_event == 1)
  )]
  st <- rbind(
    z[, .(pid, exposure01, ow, t5, endpoint = "c", event = css)],
    z[, .(pid, exposure01, ow, t5, endpoint = "o", event = oth)]
  )
  st[, endpoint := factor(endpoint, levels = c("c", "o"))]
  st[, `:=`(xc = exposure01 * (endpoint == "c"), xo = exposure01 * (endpoint == "o"))]
  fit <- coxph(
    Surv(t5, event) ~ xc + xo + strata(endpoint) + cluster(pid),
    data = st, weights = ow, robust = TRUE, ties = "efron"
  )
  b <- coef(fit)[c("xc", "xo")]
  v <- vcov(fit)[c("xc", "xo"), c("xc", "xo")]
  contrast <- unname(b[1] - b[2])
  se <- sqrt(as.numeric(c(1, -1) %*% v %*% c(1, -1)))
  data.table(
    cohort = cohort, analysis = analysis, n = nrow(d),
    log_hr_cancer = unname(b[1]), log_hr_other = unname(b[2]),
    covariance = unname(v[1, 2]), log_ratio = contrast,
    robust_se = se, ratio = exp(contrast),
    lower = exp(contrast - 1.96 * se), upper = exp(contrast + 1.96 * se)
  )
}

piecewise_joint <- function(d, cohort, cut = 2) {
  base <- d[, .(
    pid = .I, exposure01, ow,
    t5 = pmin(surv_days / 365.25, 5),
    css = as.integer(surv_days <= 5 * 365.25 & css_event == 1),
    oth = as.integer(surv_days <= 5 * 365.25 & other_event == 1)
  )]
  stacked <- rbind(
    base[, .(pid, exposure01, ow, t5, endpoint = "cancer", event = css)],
    base[, .(pid, exposure01, ow, t5, endpoint = "other", event = oth)]
  )
  split <- as.data.table(survSplit(
    Surv(t5, event) ~ ., data = as.data.frame(stacked),
    cut = cut, episode = "period", start = "tstart"
  ))
  split[, period_label := factor(ifelse(period == 1, "0-2 years", "2-5 years"),
                                 levels = c("0-2 years", "2-5 years"))]
  split[, ep_period := interaction(endpoint, period_label, drop = TRUE)]
  split[, `:=`(
    c_early = exposure01 * (endpoint == "cancer" & period == 1),
    c_late = exposure01 * (endpoint == "cancer" & period == 2),
    o_early = exposure01 * (endpoint == "other" & period == 1),
    o_late = exposure01 * (endpoint == "other" & period == 2)
  )]
  fit <- coxph(
    Surv(tstart, t5, event) ~ c_early + c_late + o_early + o_late +
      strata(ep_period) + cluster(pid),
    data = split, weights = ow, robust = TRUE, ties = "efron"
  )
  b <- coef(fit)[c("c_early", "c_late", "o_early", "o_late")]
  v <- vcov(fit)[names(b), names(b)]
  rows <- lapply(list(c("c_early", "o_early", "0-2 years"),
                      c("c_late", "o_late", "2-5 years")), function(keys) {
    bc <- unname(b[keys[1]]); bo <- unname(b[keys[2]])
    idx <- match(keys[1:2], names(b))
    contrast <- bc - bo
    se_ratio <- sqrt(as.numeric(c(1, -1) %*% v[idx, idx] %*% c(1, -1)))
    sec <- sqrt(v[idx[1], idx[1]]); seo <- sqrt(v[idx[2], idx[2]])
    data.table(
      cohort = cohort, period = keys[3], cut_years = cut,
      cancer_events = split[endpoint == "cancer" & period_label == keys[3], sum(event)],
      other_events = split[endpoint == "other" & period_label == keys[3], sum(event)],
      persons_at_period_start = uniqueN(split[period_label == keys[3], pid]),
      log_hr_cancer = bc, se_cancer = sec, hr_cancer = exp(bc),
      cancer_lower = exp(bc - 1.96 * sec), cancer_upper = exp(bc + 1.96 * sec),
      log_hr_other = bo, se_other = seo, hr_other = exp(bo),
      other_lower = exp(bo - 1.96 * seo), other_upper = exp(bo + 1.96 * seo),
      covariance = v[idx[1], idx[2]], log_ratio = contrast, se_ratio = se_ratio,
      ratio = exp(contrast), ratio_lower = exp(contrast - 1.96 * se_ratio),
      ratio_upper = exp(contrast + 1.96 * se_ratio)
    )
  })
  rbindlist(rows)
}

baseline_long <- function(d, cohort) {
  rows <- list()
  for (vn in c("age", "year")) {
    for (g in levels(d$exposure)) {
      z <- d[exposure == g]
      rows[[length(rows) + 1L]] <- data.table(
        cohort = cohort, variable = vn, level = "mean (SD)", exposure = g,
        n = nrow(z), unweighted_value = sprintf("%.2f (%.2f)", mean(z[[vn]]), sd(z[[vn]])),
        weighted_value = sprintf("%.2f", weighted.mean(z[[vn]], z$ow)),
        crude_smd = smd(d[[vn]], d$exposure01, rep(1, nrow(d))),
        weighted_smd = smd(d[[vn]], d$exposure01, d$ow)
      )
    }
  }
  for (vn in c("sex", "race", "marital", "sitegrp", "sizegrp", "gradegrp", "histgrp")) {
    for (lv in levels(droplevels(d[[vn]]))) {
      x <- as.numeric(d[[vn]] == lv)
      for (g in levels(d$exposure)) {
        z <- d[exposure == g]
        zg <- as.numeric(z[[vn]] == lv)
        rows[[length(rows) + 1L]] <- data.table(
          cohort = cohort, variable = vn, level = lv, exposure = g,
          n = sum(zg), unweighted_value = sprintf("%d (%.1f%%)", sum(zg), 100 * mean(zg)),
          weighted_value = sprintf("%.1f%%", 100 * weighted.mean(zg, z$ow)),
          crude_smd = smd(x, d$exposure01, rep(1, nrow(d))),
          weighted_smd = smd(x, d$exposure01, d$ow)
        )
      }
    }
  }
  rbindlist(rows)
}

weighted_cif5 <- function(d, cohort) {
  d <- copy(d)
  d[, `:=`(
    t5 = pmin(surv_days / 365.25, 5),
    ec = as.integer(surv_days <= 5 * 365.25 & css_event == 1),
    eo = as.integer(surv_days <= 5 * 365.25 & other_event == 1)
  )]
  rbindlist(lapply(levels(d$exposure), function(g) {
    z <- d[exposure == g]
    times <- sort(unique(z[(ec == 1 | eo == 1) & t5 <= 5, t5]))
    surv <- 1; cifc <- 0; cifo <- 0
    for (tt in times) {
      risk <- z[t5 >= tt, sum(ow)]
      dc <- z[t5 == tt & ec == 1, sum(ow)]
      do <- z[t5 == tt & eo == 1, sum(ow)]
      if (is.finite(risk) && risk > 0) {
        cifc <- cifc + surv * dc / risk
        cifo <- cifo + surv * do / risk
        surv <- surv * (1 - (dc + do) / risk)
      }
    }
    data.table(cohort = cohort, exposure = g, cancer_cif_5y = cifc,
               other_cif_5y = cifo, overall_survival_5y = surv)
  }))
}

diagnostics <- function(obj, cohort) {
  d <- obj$data
  crude_weights <- rep(1, nrow(d))
  raw_smd <- c(smd(d$age, d$exposure01, crude_weights), smd(d$year, d$exposure01, crude_weights))
  weighted_smd <- c(smd(d$age, d$exposure01, d$ow), smd(d$year, d$exposure01, d$ow))
  for (vn in c("sex", "race", "marital", "sitegrp", "sizegrp", "gradegrp", "histgrp")) {
    for (lv in levels(droplevels(d[[vn]]))) {
      x <- as.numeric(d[[vn]] == lv)
      raw_smd <- c(raw_smd, smd(x, d$exposure01, crude_weights))
      weighted_smd <- c(weighted_smd, smd(x, d$exposure01, d$ow))
    }
  }
  follow <- survfit(Surv(surv_days / 365.25, 1 - os_event) ~ 1, data = d)
  median_follow <- unname(summary(follow)$table["median"])
  data.table(
    cohort = cohort, n = nrow(d), median_reverse_km_followup_years = median_follow,
    max_crude_smd = max(raw_smd, na.rm = TRUE),
    max_weighted_smd = max(weighted_smd, na.rm = TRUE),
    unknown_size_n = sum(is.na(d$size)), unknown_size_pct = mean(is.na(d$size)),
    unknown_grade_n = sum(d$gradegrp == "Unknown"), unknown_grade_pct = mean(d$gradegrp == "Unknown"),
    nx_n = sum(d$N == "NX"), nx_pct = mean(d$N == "NX")
  )
}

run_cohort <- function(path, cohort) {
  obj <- reconstruct_keep(path, cohort)
  d <- obj$data
  primary_models <- rbind(
    endpoint_model(d, "css_event", "cancer death", cohort, "primary unknown-category"),
    endpoint_model(d, "other_event", "other-cause death", cohort, "primary unknown-category")
  )
  primary_joint <- joint_model(d, cohort, "primary unknown-category")
  cif <- weighted_cif5(d, cohort)
  diag <- diagnostics(obj, cohort)
  ess <- copy(obj$ess); ess[, cohort := cohort]
  unknown <- copy(obj$unknown); unknown[, cohort := cohort]
  flow <- copy(obj$flow); setnames(flow, c("V1", "V2"), c("step", "n")); flow[, cohort := cohort]

  known <- copy(d[!is.na(size) & gradegrp != "Unknown"])
  known <- fit_ps(known)
  known_models <- rbind(
    endpoint_model(known, "css_event", "cancer death", cohort, "known size and grade"),
    endpoint_model(known, "other_event", "other-cause death", cohort, "known size and grade")
  )
  known_joint <- joint_model(known, cohort, "known size and grade")
  known_summary <- data.table(
    cohort = cohort, primary_n = nrow(d), known_size_grade_n = nrow(known),
    retained_fraction = nrow(known) / nrow(d)
  )
  piece <- piecewise_joint(d, cohort, 2)

  event_counts <- d[, .(
    n = .N,
    cancer_deaths_5y = sum(surv_days <= 5 * 365.25 & css_event == 1),
    other_deaths_5y = sum(surv_days <= 5 * 365.25 & other_event == 1),
    all_deaths_5y = sum(surv_days <= 5 * 365.25 & os_event == 1)
  ), by = exposure]
  event_counts[, cohort := cohort]
  surgery_codes <- d[, .(n = .N), by = .(exposure, surg_code)][order(exposure, surg_code)]
  surgery_codes[, `:=`(cohort = cohort, percent_within_group = 100 * n / sum(n)), by = exposure]
  baseline <- baseline_long(d, cohort)

  rm(d, known); gc()
  list(models = rbind(primary_models, known_models),
       joint = rbind(primary_joint, known_joint), cif = cif, diagnostics = diag,
       ess = ess, unknown = unknown, flow = flow, known_summary = known_summary,
       piecewise = piece, event_counts = event_counts,
       surgery_codes = surgery_codes, baseline = baseline)
}

crc <- run_cohort(normalizePath(args[[1]], mustWork = TRUE), "CRC")
bladder <- run_cohort(normalizePath(args[[2]], mustWork = TRUE), "Bladder")

fwrite(rbind(crc$models, bladder$models), file.path(out_dir, "seer_endpoint_models.csv"))
fwrite(rbind(crc$joint, bladder$joint), file.path(out_dir, "seer_joint_contrasts.csv"))
fwrite(rbind(crc$cif, bladder$cif), file.path(out_dir, "seer_weighted_cif_5y.csv"))
fwrite(rbind(crc$diagnostics, bladder$diagnostics), file.path(out_dir, "seer_diagnostics.csv"))
fwrite(rbind(crc$ess, bladder$ess), file.path(out_dir, "seer_overlap_ess.csv"))
fwrite(rbind(crc$unknown, bladder$unknown), file.path(out_dir, "seer_unknown_by_group.csv"))
fwrite(rbind(crc$flow, bladder$flow), file.path(out_dir, "seer_cohort_flow.csv"))
fwrite(rbind(crc$known_summary, bladder$known_summary), file.path(out_dir, "seer_known_size_grade_sensitivity_counts.csv"))
fwrite(rbind(crc$piecewise, bladder$piecewise), file.path(out_dir, "seer_piecewise_0_2_5_joint.csv"))
fwrite(rbind(crc$event_counts, bladder$event_counts), file.path(out_dir, "seer_event_counts_by_group.csv"))
fwrite(rbind(crc$surgery_codes, bladder$surgery_codes), file.path(out_dir, "seer_surgery_code_composition.csv"))
fwrite(rbind(crc$baseline, bladder$baseline), file.path(out_dir, "seer_baseline_characteristics_long.csv"))

joint <- rbind(crc$joint, bladder$joint)[analysis == "primary unknown-category"]
ld <- joint[cohort == "Bladder", log_ratio] - joint[cohort == "CRC", log_ratio]
se_between <- sqrt(sum(joint$robust_se^2))
between <- data.table(
  contrast = "Bladder ratio / CRC ratio", log_ratio = ld, robust_se = se_between,
  ratio = exp(ld), lower = exp(ld - 1.96 * se_between), upper = exp(ld + 1.96 * se_between),
  p_value = 2 * pnorm(-abs(ld / se_between)), interpretation = "exploratory descriptive contrast"
)
fwrite(between, file.path(out_dir, "seer_between_cohort_contrast.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
cat("Completed aggregate-only SEER sensitivity reanalysis in", out_dir, "\n")
