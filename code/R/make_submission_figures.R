suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
fig_dir <- file.path(root, "figures")
verification_dir <- file.path(root, "verification")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(verification_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#222222"
blue <- "#0072B2"
sky <- "#56B4E9"
orange <- "#D55E00"
amber <- "#E69F00"
purple <- "#6A3D9A"
green <- "#009E73"
grey <- "#6B7280"
light_grey <- "#E5E7EB"

theme_pub <- function(base_size = 8.5) {
  theme_minimal(base_family = font_family, base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6E8EB", linewidth = 0.3),
      axis.title = element_text(colour = ink, face = "plain"),
      axis.text = element_text(colour = ink),
      strip.text = element_text(face = "bold", colour = ink),
      strip.background = element_rect(fill = "#F3F4F6", colour = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = rel(1.02), hjust = 0),
      plot.subtitle = element_text(colour = "#4B5563"),
      plot.tag = element_text(face = "bold", size = rel(1.2), colour = ink),
      plot.margin = margin(5, 6, 5, 6)
    )
}

save_bundle <- function(plot, stem, width, height) {
  png_file <- file.path(fig_dir, paste0(stem, ".png"))
  tif_file <- file.path(fig_dir, paste0(stem, ".tiff"))
  pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(fig_dir, paste0(stem, ".svg"))

  agg_png(png_file, width = width, height = height, units = "in", res = 300,
          background = "white", scaling = 1)
  print(plot)
  dev.off()

  agg_tiff(tif_file, width = width, height = height, units = "in", res = 600,
           compression = "lzw", background = "white", scaling = 1)
  print(plot)
  dev.off()

  cairo_pdf(pdf_file, width = width, height = height, family = "Arial",
            onefile = FALSE)
  print(plot)
  dev.off()

  svg(svg_file, width = width, height = height, family = "Arial",
      onefile = FALSE, bg = "white")
  print(plot)
  dev.off()
}

panel_tags <- theme(plot.tag.position = c(0.01, 0.99),
                    plot.tag = element_text(face = "bold", size = 11,
                                            hjust = 0, vjust = 1))

# Figure 1: approved standalone vector workflow ---------------------------
figure1_env <- new.env(parent=globalenv())
sys.source(file.path(root, "code", "R", "make_figure1.R"), envir=figure1_env)
figure1_env$make_figure1(root)

# Figure 2: main-family bias and coverage ------------------------------------
primary <- fread(file.path(root, "analysis", "strict_cox_rerun",
                           "strict_primary_summary_72_rows.csv"))
d2 <- primary[family == "main_hi" & abs(eff - 0.85) < 1e-10]
d2[, se_sp := factor(sprintf("%.2f / %.2f", se, sp),
                     levels = c("1.00 / 1.00", "0.90 / 0.90", "0.75 / 0.85"))]
d2[, uhr_f := factor(sprintf("%.2f", uhr), levels = c("1.00", "0.50", "0.33"))]

bias_panel <- function(n_value) {
  z <- d2[n == n_value]
  ggplot(z, aes(se_sp, uhr_f, fill = median_HR_percent_error_c)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%+.2f%%", median_HR_percent_error_c),
                  colour = abs(median_HR_percent_error_c) > 8),
              size = 3.0, family = font_family, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white")) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
                         midpoint = 0, limits = c(-16, 8), oob = squish,
                         name = "Relative HR\nerror (%)") +
    labs(title = paste0("Relative error; n = ", comma(n_value)),
         x = "Treatment-code Se / Sp", y = "Health-reserve HR (uHR)") +
    coord_equal() + theme_pub() +
    theme(panel.grid = element_blank(), legend.key.height = unit(3.5, "mm"))
}

coverage_panel <- function(n_value) {
  z <- d2[n == n_value]
  ggplot(z, aes(se_sp, uhr_f, fill = coverage_c)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", coverage_c),
                  colour = coverage_c < 0.35),
              size = 3.0, family = font_family, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white")) +
    scale_fill_gradientn(colours = c("#7F0000", "#D7301F", "#FEE8C8", "#D8ECF3", "#0072B2"),
                         values = rescale(c(0, 0.25, 0.5, 0.8, 0.95)), limits = c(0, 0.96),
                         oob = squish, name = "95% coverage") +
    labs(title = paste0("Coverage; n = ", comma(n_value)),
         x = "Treatment-code Se / Sp", y = "Health-reserve HR (uHR)") +
    coord_equal() + theme_pub() +
    theme(panel.grid = element_blank(), legend.key.height = unit(3.5, "mm"))
}

legend_space <- theme(
  legend.position="bottom", legend.title=element_text(margin=margin(b=6)),
  legend.margin=margin(t=8,r=6,b=10,l=6), legend.box.spacing=unit(3,"mm"),
  plot.margin=margin(t=8,r=10,b=6,l=10))
bar_guide <- guides(fill=guide_colorbar(title.position="top", title.hjust=0.5,
                                      barwidth=unit(38,"mm"),barheight=unit(3.5,"mm")))
b1 <- bias_panel(5000) + bar_guide + legend_space
b2 <- bias_panel(50000) + bar_guide + legend_space
c1 <- coverage_panel(5000) + bar_guide + legend_space
c2 <- coverage_panel(50000) + bar_guide + legend_space
top2 <- (b1+b2+plot_layout(guides="collect")) & legend_space
bottom2 <- (c1+c2+plot_layout(guides="collect")) & legend_space
fig2 <- top2 / bottom2 + plot_annotation(tag_levels="a") & panel_tags
save_bundle(fig2,"Figure2_primary_bias_coverage",7.2,7.0)

# Figure 3: directional negative-control signal ------------------------------
sig <- fread(file.path(root, "source_data", "figure3_directional_signal_summary.csv"))
sig[, `:=`(mean = 100 * mean, min = 100 * min, max = 100 * max)]
sig[, uhr_f := factor(sprintf("%.2f", uhr), levels = c("1.00", "0.50", "0.33"),
                      labels = c("None\n1.00", "Moderate\n0.50", "Strong\n0.33"))]

signal_panel <- function(z, title, color_map = NULL) {
  if (is.null(color_map)) {
    z[, series := title]
    color_map <- setNames(blue, title)
  }
  ggplot(z, aes(uhr_f, mean, group = series, colour = series)) +
    geom_hline(yintercept = 2.5, linetype = 3, colour = grey, linewidth = 0.55) +
    geom_linerange(aes(ymin = min, ymax = max), linewidth = 0.85,
                   position = position_dodge(width = 0.20)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.20)) +
    geom_line(linewidth = 0.65, position = position_dodge(width = 0.20)) +
    scale_colour_manual(values = color_map, name = NULL) +
    scale_y_continuous(limits = c(0, 103), breaks = c(0, 25, 50, 75, 100),
                       labels = label_percent(scale = 1)) +
    labs(title = title, x = "Shared health-reserve strength (uHR)",
         y = "Protective-direction signal") +
    theme_pub() +
    theme(legend.position = if (length(color_map) > 1) "bottom" else "none",
          plot.title = element_text(margin = margin(l = 13)))
}

za <- copy(sig[family == "main_hi" & n == 5000]); za[, series := "Main"]
zb <- copy(sig[family == "main_hi" & n == 50000]); zb[, series := "Main"]
zc <- copy(sig[family %chin% c("low_crc", "low_bladder") & n == 50000])
zc[, series := fifelse(family == "low_crc", "CRC anchor", "Bladder anchor")]
p3a <- signal_panel(za, "Main family; n = 5,000", c(Main = blue))
p3b <- signal_panel(zb, "Main family; n = 50,000", c(Main = blue))
p3c <- signal_panel(zc, "Case-informed low-event families; n = 50,000",
                    c(`CRC anchor` = green, `Bladder anchor` = orange))
p3a <- p3a + labs(title="Main family\nn = 5,000")
p3b <- p3b + labs(title="Main family\nn = 50,000")
p3c <- p3c + labs(title="Case-informed low-event\nfamilies; n = 50,000")
fig3 <- p3a + plot_spacer() + p3b + plot_spacer() + p3c +
  plot_layout(ncol=5,widths=c(1,.12,1,.12,1),guides="collect",axes="collect_y",axis_titles="collect") +
  plot_annotation(tag_levels="a") & panel_tags &
  theme(legend.position="bottom",plot.margin=margin(t=8,r=8,b=8,l=8),plot.title=element_text(margin=margin(l=13,b=9)))
save_bundle(fig3,"Figure3_directional_signal",7.2,3.75)

# Figure 4: QBA ---------------------------------------------------------------
qba <- fread(file.path(root, "analysis", "strict_qba", "strict_qba_summary_36_rows.csv"))
qba[, prior_label := factor(prior,
  levels = c("oracle_cell_specific", "global_mechanism_prior", "oracle_mean_x1.25", "oracle_mean_x1.50"),
  labels = c("Cell oracle", "Global mechanism", "Oracle mean x1.25", "Oracle mean x1.50"))]
qba[, uhr_f := factor(sprintf("%.2f", uhr), levels = c("1.00", "0.50", "0.33"))]
qba[, se_sp := factor(sprintf("Se/Sp %.2f/%.2f", se, sp),
                      levels = c("Se/Sp 1.00/1.00", "Se/Sp 0.90/0.90", "Se/Sp 0.75/0.85"))]
prior_cols <- c("Cell oracle" = blue, "Global mechanism" = orange,
                "Oracle mean x1.25" = green, "Oracle mean x1.50" = purple)

p4a <- ggplot(qba, aes(uhr_f, coverage, colour = prior_label, group = prior_label)) +
  geom_hline(yintercept = 0.95, linetype = 3, colour = grey, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, coverage - 1.96 * mcse_coverage),
                    ymax = pmin(1, coverage + 1.96 * mcse_coverage)),
                width = 0.12, position = position_dodge(width = 0.32), linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.32), size = 1.8) +
  geom_line(position = position_dodge(width = 0.32), linewidth = 0.45) +
  facet_wrap(~se_sp, nrow = 1) +
  scale_colour_manual(values = prior_cols, name = "Bias distribution") +
  scale_y_continuous(limits = c(0.2, 1.01), breaks = c(0.25, 0.5, 0.75, 0.95),
                     labels = label_percent()) +
  labs(title = "Grid-level 95% coverage", x = "uHR", y = "Coverage") +
  theme_pub(8) + theme(legend.position = "bottom")

qsum <- qba[, .(bias = mean(bias_loghr), rmse = mean(rmse_loghr),
                coverage = mean(coverage), width = mean(width_ratio)), by = prior_label]
p4b <- ggplot(qsum, aes(bias, prior_label, colour = prior_label)) +
  geom_vline(xintercept = 0, linetype = 3, colour = grey) +
  geom_segment(aes(x = 0, xend = bias, yend = prior_label), linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sprintf("%+.3f", bias)), nudge_y = 0.20,
            size = 2.6, family = font_family, show.legend = FALSE) +
  scale_colour_manual(values = prior_cols, guide = "none") +
  scale_x_continuous(limits = c(-0.006, 0.031), breaks = c(0, 0.01, 0.02, 0.03)) +
  labs(title = "Pooled log-HR bias", x = "Bias", y = NULL) + theme_pub(8)

p4c <- ggplot(qsum, aes(width, prior_label, colour = prior_label)) +
  geom_vline(xintercept = 1, linetype = 3, colour = grey) +
  geom_segment(aes(x = 1, xend = width, yend = prior_label), linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sprintf("%.3f", width)), nudge_y = 0.20,
            size = 2.6, family = font_family, show.legend = FALSE) +
  scale_colour_manual(values = prior_cols, guide = "none") +
  scale_x_continuous(limits = c(0.98, 1.24), breaks = c(1, 1.1, 1.2)) +
  labs(title = "Pooled interval-width ratio", x = "Adjusted / unadjusted", y = NULL) +
  theme_pub(8)

fig4 <- p4a / (p4b + p4c) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "a") & panel_tags
save_bundle(fig4, "Figure4_QBA_performance", 7.2, 5.7)

# Figure 5: H2 and pRMST ------------------------------------------------------
h2 <- fread(file.path(root, "analysis", "strict_cox_rerun", "strict_h2_summary_7_rows.csv"))
h2[, scenario_label := factor(name, levels = rev(name), labels = rev(c(
  "No direct effect; confounding", "Early harm", "Sustained harm",
  "Early protection", "Sustained protection", "Early harm; low event", "Null; no confounding")))]
h2long <- melt(h2, id.vars = c("scenario_id", "scenario_label", "median_hr_o"),
               measure.vars = c("two_sided_detection", "protective_direction_signal", "harmful_direction_signal"),
               variable.name = "rule", value.name = "rate")
h2long[, rate := 100 * rate]
h2long[, rule_label := factor(rule,
  levels = c("two_sided_detection", "protective_direction_signal", "harmful_direction_signal"),
  labels = c("Two-sided", "Protective direction", "Harmful direction"))]
rule_cols <- c("Two-sided" = ink, "Protective direction" = blue, "Harmful direction" = orange)
p5a <- ggplot(h2long, aes(rate, scenario_label, colour = rule_label, shape = rule_label)) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  scale_colour_manual(values = rule_cols, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_x_continuous(limits = c(0, 104), breaks = c(0, 25, 50, 75, 100),
                     labels = label_percent(scale = 1)) +
  labs(title = "Negative-control boundary", x = "CI exclusion rate", y = NULL) +
  theme_pub(8) + theme(legend.position = "bottom", axis.text.y = element_text(size = 7.2))

rmst <- fread(file.path(root, "analysis", "rmst_standard_ato",
                        "rmst_standard_ato_rescored_12_scenarios.csv"))
core <- rmst[archived_tag == "core"]
core[, se_sp := factor(sprintf("%.2f / %.2f", se, sp),
                       levels = c("1.00 / 1.00", "0.90 / 0.90", "0.75 / 0.85"))]
core[, uhr_f := factor(sprintf("%.2f", uhr), levels = c("1.00", "0.50", "0.33"))]
p5b <- ggplot(core, aes(se_sp, uhr_f, fill = bias)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.3f", bias), colour = abs(bias) > 0.055),
            size = 2.8, family = font_family, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = ink, `TRUE` = "white")) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
                       midpoint = 0, limits = c(-0.09, 0.09), name = "Bias (years)") +
  labs(title = "Core 5-year pRMST bias", x = "Treatment-code Se / Sp", y = "uHR") +
  coord_equal() + theme_pub(8) + theme(panel.grid = element_blank())

rmst[, scenario_short := factor(scenario_id, levels = scenario_id,
  labels = c(paste0("C", 1:9), "EB", "DI", "CR"))]
covlong <- melt(rmst, id.vars = c("scenario_id", "scenario_short"),
                measure.vars = c("coverage", "bias_shifted_coverage_descriptive"),
                variable.name = "coverage_type", value.name = "coverage_value")
covlong[, coverage_label := factor(coverage_type,
  levels = c("coverage", "bias_shifted_coverage_descriptive"),
  labels = c("Raw", "Bias-recentered diagnostic"))]
p5c <- ggplot(covlong, aes(scenario_short, coverage_value, colour = coverage_label,
                           group = coverage_label)) +
  geom_hline(yintercept = 0.95, linetype = 3, colour = grey, linewidth = 0.5) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.8) +
  geom_vline(xintercept = 9.5, linetype = 2, colour = "#A0A4A8", linewidth = 0.45) +
  scale_colour_manual(values = c("Raw" = orange, "Bias-recentered diagnostic" = blue), name = NULL) +
  scale_y_continuous(limits = c(0.55, 1.005), breaks = c(0.6, 0.7, 0.8, 0.9, 0.95, 1),
                     labels = label_percent()) +
  labs(title = "Raw versus bias-recentered coverage", x = "Scenario", y = "Coverage") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.8),
        legend.position = "bottom")

fig5 <- p5a + (p5b / p5c) + plot_layout(widths = c(1.05, 1.35)) +
  plot_annotation(tag_levels = "a") & panel_tags
save_bundle(fig5, "Figure5_H2_pRMST_boundaries", 7.2, 6.0)

# Figure 6: SEER cohort flow and descriptive models --------------------------
seer_dir <- file.path(root, "analysis", "seer_descriptive_sensitivity")
flow <- fread(file.path(seer_dir, "seer_cohort_flow.csv"))
flow_levels <- c("T1", "first", "hist", "hist_confirm", "site", "M0", "N0NX",
                 "no_therapy", "adult_survival", "surgery")
flow_labels <- c("T1", "First\nprimary", "Eligible\nhistology", "Histology\nconfirmed",
                 "Eligible\nsite", "M0", "N0/NX", "No prior\ntherapy",
                 "Adult +\nsurvival", "Eligible\nsurgery")
flow[, step_f := factor(step, levels = flow_levels, labels = flow_labels)]
cohort_cols <- c("CRC" = green, "Bladder" = purple)
p6a <- ggplot(flow, aes(step_f, n, colour = cohort, group = cohort)) +
  geom_line(linewidth = 0.75) + geom_point(size = 2.0) +
  geom_text(aes(label = comma(n)), vjust = -0.75, size = 2.25,
            family = font_family, show.legend = FALSE) +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 125000), breaks = seq(0, 120000, 30000), labels = comma) +
  labs(title = "Post-T1 cohort selection", x = NULL, y = "Records retained") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.2),
                        legend.position = "bottom")

models <- fread(file.path(seer_dir, "seer_endpoint_models.csv"))
models[, endpoint_label := factor(endpoint,
  levels = c("cancer death", "other-cause death"),
  labels = c("Cancer death", "Other-cause death"))]
models[, analysis_label := factor(analysis,
  levels = c("primary unknown-category", "known size and grade"),
  labels = c("Primary", "Known size + grade"))]
models[, row_label := factor(paste(cohort, analysis_label, endpoint_label, sep = " | "),
                             levels = rev(paste(cohort, analysis_label, endpoint_label, sep = " | ")))]
p6b <- ggplot(models, aes(hr, row_label, colour = endpoint_label, shape = analysis_label)) +
  geom_vline(xintercept = 1, linetype = 3, colour = grey) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15,
                linewidth = 0.55, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("Cancer death" = blue, "Other-cause death" = orange), name = "Endpoint") +
  scale_shape_manual(values = c("Primary" = 16, "Known size + grade" = 17), name = "Analysis") +
  scale_x_log10(limits = c(0.50, 1.40), breaks = c(0.5, 0.7, 1, 1.4)) +
  labs(title = "Observed-code cause-specific HRs", x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub(7.8) + theme(axis.text.y = element_text(size = 6.6),
                         plot.title.position = "plot",
                         plot.title = element_text(hjust = 0, size = 7.7,
                                                   margin = margin(l = 18, b = 2)),
                         legend.position = "bottom")

joint <- fread(file.path(seer_dir, "seer_joint_contrasts.csv"))
j1 <- joint[, .(cohort, analysis_label = fifelse(analysis == "primary unknown-category", "Overall", "Known size + grade"),
                ratio, lower, upper)]
piece <- fread(file.path(seer_dir, "seer_piecewise_0_2_5_joint.csv"))
j2 <- piece[, .(cohort, analysis_label = period, ratio, lower = ratio_lower, upper = ratio_upper)]
rat <- rbind(j1, j2, fill = TRUE)
rat[, analysis_label := factor(analysis_label,
  levels = c("Overall", "known", "Known size + grade", "0-2 years", "2-5 years"))]
rat[, row_label := factor(paste(cohort, analysis_label, sep = " | "),
                          levels = rev(paste(cohort, analysis_label, sep = " | ")))]
p6c <- ggplot(rat, aes(ratio, row_label, colour = cohort)) +
  geom_vline(xintercept = 1, linetype = 3, colour = grey) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15,
                linewidth = 0.55, orientation = "y") +
  geom_point(size = 2.2) +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  scale_x_log10(limits = c(0.60, 1.85), breaks = c(0.6, 0.8, 1, 1.3, 1.8)) +
  labs(title = "Cancer-to-other HR ratios", x = "HR ratio (log scale)", y = NULL) +
  theme_pub(7.8) + theme(axis.text.y = element_text(size = 6.6),
                         plot.title.position = "plot",
                         plot.title = element_text(hjust = 0, size = 7.7,
                                                   margin = margin(l = 18, b = 2)),
                         legend.position = "bottom")

fig6 <- wrap_elements(full=p6a) / (p6b+p6c+plot_layout(widths=c(1.08,.92))) +
  plot_layout(heights=c(.9,1.35)) + plot_annotation(tag_levels="a") & panel_tags
save_bundle(fig6,"Figure6_SEER_descriptive",7.2,6.7)
# Export concise source-data checks used by the figures.
checks <- data.table(
  figure = paste0("Figure", 1:6),
  png = file.path(fig_dir, paste0(c("Figure1_DGM_observation_boundary", "Figure2_primary_bias_coverage",
    "Figure3_directional_signal", "Figure4_QBA_performance", "Figure5_H2_pRMST_boundaries",
    "Figure6_SEER_descriptive"), ".png")),
  tiff = file.path(fig_dir, paste0(c("Figure1_DGM_observation_boundary", "Figure2_primary_bias_coverage",
    "Figure3_directional_signal", "Figure4_QBA_performance", "Figure5_H2_pRMST_boundaries",
    "Figure6_SEER_descriptive"), ".tiff")),
  pdf = file.path(fig_dir, paste0(c("Figure1_DGM_observation_boundary", "Figure2_primary_bias_coverage",
    "Figure3_directional_signal", "Figure4_QBA_performance", "Figure5_H2_pRMST_boundaries",
    "Figure6_SEER_descriptive"), ".pdf")),
  svg = file.path(fig_dir, paste0(c("Figure1_DGM_observation_boundary", "Figure2_primary_bias_coverage",
    "Figure3_directional_signal", "Figure4_QBA_performance", "Figure5_H2_pRMST_boundaries",
    "Figure6_SEER_descriptive"), ".svg"))
)
checks[, c("png", "tiff", "pdf", "svg") := lapply(.SD, function(x) file.path("figures", basename(x))), .SDcols=c("png", "tiff", "pdf", "svg")]
fwrite(checks, file.path(verification_dir, "figure_export_manifest.csv"))
cat("Generated", nrow(checks), "figure bundles in", fig_dir, "\n")
