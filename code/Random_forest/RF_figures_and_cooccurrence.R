# ============================================================
# Figures, tables, and SNP co-occurrence for the random forest analysis.
#
# Run FINAL_RF_10_29_26.R first. This script reads its outputs from
# FINAL_RF_panC_stratified_70_30/ and writes every figure (PNG + SVG)
# and table to FINAL_RF_panC_stratified_70_30/figures/.
#
# Figures
#   Fig3. ROC curve on the test set, with AUC recomputed from the predictions
#   Fig4. Top 15 SNPs: heatmap of per-SNP classification metrics
#   Fig5. Top 15 SNPs: difference in mean cytotoxicity (present minus absent)
#   Fig6. Top 15 SNPs: co-occurrence (Pearson r, binary presence/absence)
#   Fig8. Top 15 SNPs: % of each panC group with the SNP (heatmap)
#
# Per-SNP statistics (Fig4, Fig5) use all 270 isolates, the same
# descriptive view as the paper figure. They are not out-of-sample results,
# and the SNPs were selected using these same isolates.
# ============================================================

library(ggplot2)

base_folder <- path.expand(paste0(
  "~/Library/CloudStorage/",
  "OneDrive-ThePennsylvaniaStateUniversity/",
  "Tyler's projects - Collaboration and External Projects/",
  "Cassidy_analysis/Random Forrest/RF_analysis"
))

model <- "SNP_status_no_panC"
rf_dir <- file.path(base_folder, "FINAL_RF_panC_stratified_70_30")
fig_dir <- file.path(rf_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

n_top <- 15
cytotoxicity_cutoff <- 0.70

# Every figure is saved as PNG and SVG
save_plot <- function(p, name, width, height, dpi = 600) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p,
         width = width, height = height, bg = "white")
}

theme_report <- theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))

# SNP names as in the model inputs, with the gene italic:
# hblB_SNP_T_270 -> <i>hblB</i>_SNP_T_270 (ggtext renders the <i> tag)
italic_gene_labels <- function(x) {
  sub("^(hbl[A-D]|nhe[A-C])", "<i>\\1</i>", x)
}

# Bold SNP names, with the gene italic: <b><i>hblB</i>_SNP_T_270</b>
bold_italic_labels <- function(x) {
  paste0("<b>", sub("^(hbl[A-D]|nhe[A-C])", "<i>\\1</i>", x), "</b>")
}


# ============================================================
# 1. LOAD MODEL OUTPUTS
# ============================================================

class_metrics <- read.csv(
  file.path(rf_dir, paste0("RF_", model, "_classification_metrics.csv")))
reg_metrics <- read.csv(
  file.path(rf_dir, paste0("RF_", model, "_test_performance.csv")))
model_summary <- read.csv(
  file.path(rf_dir, paste0("RF_", model, "_summary.csv")))
test_preds <- read.csv(
  file.path(rf_dir, paste0(model, "_classification_test_predictions.csv")))
importance <- read.csv(
  file.path(rf_dir, paste0(model, "_bootstrap_importance.csv")))

cm <- class_metrics[1, ]


# ============================================================
# 2. REPORT TABLE: ALL OUTPUTS
# ============================================================

report <- rbind(
  data.frame(Section = "Tuning and sample size",
             Metric = c("N_isolates", "N_training", "N_test", "N_predictors",
                        "mtry", "min.node.size", "num.trees",
                        "Best_tuning_CV_MSE"),
             Value = as.numeric(model_summary[1, c("N_isolates", "N_training",
                                                   "N_test", "N_predictors",
                                                   "mtry", "min.node.size",
                                                   "num.trees",
                                                   "Best_tuning_CV_MSE")])),
  data.frame(Section = "Regression (test set)",
             Metric = c("MSE", "RMSE", "MAE", "R2", "Bias"),
             Value = as.numeric(reg_metrics[1, c("MSE", "RMSE", "MAE", "R2", "Bias")])),
  data.frame(Section = "Classification (test set)",
             Metric = c("TP", "TN", "FP", "FN", "Accuracy", "Sensitivity",
                        "Specificity", "Precision", "F1", "Kappa", "AUC",
                        "N_high", "N_low"),
             Value = as.numeric(cm[1, c("TP", "TN", "FP", "FN", "Accuracy",
                                        "Sensitivity", "Specificity", "Precision",
                                        "F1", "Kappa", "AUC", "N_high", "N_low")]))
)

write.csv(report, file.path(fig_dir, "RF_report_all_outputs.csv"),
          row.names = FALSE)


# ============================================================
# 3. FIG 3: ROC CURVE
# ============================================================

# Trace the ROC from every observed probability as a threshold.
thresholds <- sort(unique(c(0, test_preds$Probability_high, 1)), decreasing = TRUE)

roc <- do.call(rbind, lapply(thresholds, function(t) {
  called_high <- test_preds$Probability_high >= t
  data.frame(
    FPR = sum(called_high & test_preds$Observed_class == "Low") /
      sum(test_preds$Observed_class == "Low"),
    TPR = sum(called_high & test_preds$Observed_class == "High") /
      sum(test_preds$Observed_class == "High")
  )
}))

# Trapezoidal AUC, checked against the value reported by the model
auc_check <- sum(diff(roc$FPR) * (head(roc$TPR, -1) + tail(roc$TPR, -1)) / 2)
cat(sprintf("AUC reported: %.3f | AUC from ROC curve: %.3f\n",
            cm$AUC, abs(auc_check)))

p_roc <- ggplot(roc, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_step(color = "#648ED1", linewidth = 0.8) +
  annotate("text", x = 0.62, y = 0.15,
           label = sprintf("AUC = %.2f", cm$AUC), size = 4.5, fontface = "bold") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "False positive rate (1 - specificity)",
       y = "True positive rate (sensitivity)") +
  theme_report

save_plot(p_roc, "Fig3_ROC_curve", width = 5, height = 5)


# ============================================================
# 4. TOP 15 SNPS: PER-SNP STATISTICS
# ============================================================

pres <- read.csv(
  file.path(base_folder, "snps_all_signif_pres_abs_102126_no_panC_sig.csv"),
  check.names = FALSE, na.strings = c("NA", ""), stringsAsFactors = FALSE)

safe_div <- function(a, b) if (b == 0) NA_real_ else a / b

# Predict High when the SNP is present (1), Low when absent (0).
snp_stats <- function(snp) {

  present <- as.numeric(pres[[snp]])
  cyto <- as.numeric(pres$cytotoxicity)

  ok <- !is.na(present)
  present <- present[ok]
  cyto <- cyto[ok]

  predicted_high <- present == 1
  actual_high <- cyto >= cytotoxicity_cutoff

  TP <- sum(predicted_high & actual_high)
  TN <- sum(!predicted_high & !actual_high)
  FP <- sum(predicted_high & !actual_high)
  FN <- sum(!predicted_high & actual_high)
  n <- length(cyto)

  data.frame(
    SNP = snp,
    N = n,
    TP = TP, TN = TN, FP = FP, FN = FN,
    Accuracy = safe_div(TP + TN, n),
    Precision = safe_div(TP, TP + FP),
    Sensitivity = safe_div(TP, TP + FN),
    Specificity = safe_div(TN, TN + FP),
    Mean_cytotoxicity_present = mean(cyto[present == 1]),
    Mean_cytotoxicity_absent = mean(cyto[present == 0]),
    Cytotoxicity_difference = mean(cyto[present == 1]) - mean(cyto[present == 0]),
    SD_present = sd(cyto[present == 1]),
    SD_absent = sd(cyto[present == 0]),
    N_present = sum(present == 1),
    N_absent = sum(present == 0),
    # Standard error of the difference in means (Welch)
    Difference_SE = sqrt(sd(cyto[present == 1])^2 / sum(present == 1) +
                         sd(cyto[present == 0])^2 / sum(present == 0))
  )
}

top <- head(importance[order(-importance$Mean), ], n_top)
top_stats <- do.call(rbind, lapply(top$Variable, snp_stats))

# Labels match the SNP names in the model inputs, e.g. nheA_up_SNP_t_312
top_stats$Label <- top_stats$SNP
top_stats$Importance_rank <- seq_len(nrow(top_stats))

write.csv(top_stats, file.path(fig_dir, "top15_SNP_statistics.csv"),
          row.names = FALSE)

# Plot rows from most to least important, top to bottom
row_levels <- rev(top_stats$Label)
top_stats$Label <- factor(top_stats$Label, levels = row_levels)


# ============================================================
# 5. FIG 4: HEATMAP OF PER-SNP CLASSIFICATION METRICS
# ============================================================

metric_cols <- c("Accuracy", "Precision", "Sensitivity", "Specificity")

heat <- data.frame(
  Label = factor(rep(top_stats$Label, times = length(metric_cols)),
                 levels = row_levels),
  Metric = factor(rep(metric_cols, each = nrow(top_stats)), levels = metric_cols),
  Value = unlist(lapply(metric_cols, function(m) top_stats[[m]]))
)

p_heat <- ggplot(heat, aes(x = Metric, y = Label, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(Value), "NA", sprintf("%.2f", Value))),
            size = 3.5) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      limits = c(0, 1), na.value = "grey90") +
  scale_y_discrete(labels = italic_gene_labels) +
  labs(x = NULL, y = NULL) +
  # theme_classic, not theme_minimal: theme_minimal stops ggtext from rendering the <i> tags
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(color = "black", angle = 60,
                                   hjust = 1, vjust = 1),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "none")

save_plot(p_heat, "Fig4_top15_SNP_heatmap",
          width = 5.5, height = 0.4 * n_top + 1.5)


# ============================================================
# 6. FIG 5: DIFFERENCE IN CYTOTOXICITY
# ============================================================

p_diff <- ggplot(top_stats, aes(x = Cytotoxicity_difference, y = Label,
                                fill = Cytotoxicity_difference)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.75) +
  geom_errorbar(aes(xmin = Cytotoxicity_difference - Difference_SE,
                    xmax = Cytotoxicity_difference + Difference_SE),
                width = 0.25, linewidth = 0.4, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", guide = "none") +
  scale_y_discrete(labels = italic_gene_labels) +
  labs(x = "Difference in cytotoxicity", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = element_text(face = "bold"))

save_plot(p_diff, "Fig5_top15_SNP_cytotoxicity_difference",
          width = 5.5, height = 0.4 * n_top + 1.5)


# ============================================================
# 7. FIG 6: SNP CO-OCCURRENCE (top 15 SNPs)
# ============================================================

x <- as.matrix(pres[, top_stats$SNP, drop = FALSE])
storage.mode(x) <- "numeric"

# 0 = no variant, 1 = variant present; NA (no call) excluded pairwise.
# Pairs need at least min_n isolates and variation in both SNPs.
make_corr <- function(x, stem, min_n = 20) {
  cols <- colnames(x)
  p <- ncol(x)
  R <- matrix(NA_real_, p, p, dimnames = list(cols, cols))
  N <- matrix(0L, p, p, dimnames = list(cols, cols))
  for (i in seq_len(p)) for (j in seq_len(p)) {
    ok <- is.finite(x[, i]) & is.finite(x[, j])
    N[i, j] <- sum(ok)
    if (N[i, j] >= min_n && length(unique(x[ok, i])) > 1 &&
        length(unique(x[ok, j])) > 1) {
      R[i, j] <- cor(x[ok, i], x[ok, j], method = "pearson")
    }
  }
  write.csv(R, file.path(fig_dir, paste0(stem, "_correlations.csv")))
  write.csv(N, file.path(fig_dir, paste0(stem, "_pairwise_N.csv")))

  long <- as.data.frame(as.table(R), stringsAsFactors = FALSE)
  names(long) <- c("SNP1", "SNP2", "r")
  long$SNP1 <- factor(long$SNP1, levels = cols)
  long$SNP2 <- factor(long$SNP2, levels = rev(cols))
  fig <- ggplot(long, aes(SNP1, SNP2, fill = r)) +
    geom_tile(color = "white", linewidth = .25) +
    coord_fixed() +
    scale_fill_gradient2(low = "#2867a2", mid = "white", high = "#b33c3c",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey85",
                         name = "Pearson r") +
    scale_x_discrete(labels = bold_italic_labels) +
    scale_y_discrete(labels = bold_italic_labels) +
    labs(x = NULL, y = NULL) +
    # theme_classic, not theme_minimal: theme_minimal stops ggtext from rendering the markup
    theme_classic(base_size = 11) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = ggtext::element_markdown(angle = 60, hjust = 1, vjust = 1,
                                                 size = 12, color = "black"),
          axis.text.y = ggtext::element_markdown(size = 12, color = "black"))
  save_plot(fig, stem, width = 11, height = 10, dpi = 320)

  pairs <- as.data.frame(as.table(R), stringsAsFactors = FALSE)
  names(pairs) <- c("SNP1", "SNP2", "r")
  pairs$n <- as.vector(N)
  pairs <- pairs[match(pairs$SNP1, cols) < match(pairs$SNP2, cols) & !is.na(pairs$r), ]
  pairs <- pairs[order(-abs(pairs$r)), ]
  write.csv(pairs, file.path(fig_dir, paste0(stem, "_ranked_pairs.csv")), row.names = FALSE)
  write.csv(subset(pairs, abs(r) >= 0.8),
            file.path(fig_dir, paste0(stem, "_strong_pairs.csv")), row.names = FALSE)
  invisible(list(r = R, n = N, pairs = pairs, plot = fig))
}

cooccur <- make_corr(x, "Fig6_top15_SNP_cooccurrence")


# ============================================================
# 8. FIG 8: % OF EACH PANC GROUP WITH THE SNP
# ============================================================

group_levels <- c("Group_I", "Group_II", "Group_III", "Group_IV", "Group_V",
                  "Group_VI", "Group_VII", "Group_VIII", "Group_clarus")
stopifnot(all(pres$panC_group %in% group_levels))

# Isolates with the variant present (1) and with a call (not NA), per group
counts <- do.call(rbind, lapply(top_stats$SNP, function(s) {
  v <- as.numeric(pres[[s]])
  do.call(rbind, lapply(group_levels, function(g) {
    in_group <- pres$panC_group == g
    data.frame(SNP = s, panC_group = g,
               Present = sum(v[in_group] == 1, na.rm = TRUE),
               Tested = sum(!is.na(v[in_group])))
  }))
}))
counts$Proportion <- ifelse(counts$Tested > 0, counts$Present / counts$Tested, NA_real_)
counts$Percent <- round(100 * counts$Proportion)
counts$SNP <- factor(counts$SNP, levels = row_levels)
counts$panC_group <- factor(counts$panC_group, levels = group_levels)
write.csv(counts, file.path(fig_dir, "top15_SNP_by_panC_group_counts.csv"), row.names = FALSE)

# Heatmap: same percent. NA (no isolates with a call) is labelled NA and drawn as 0.
fig_heat <- ggplot(counts, aes(panC_group, SNP, fill = Proportion)) +
  geom_tile(color = "white", linewidth = .4) +
  geom_text(aes(label = ifelse(is.na(Percent), "NA", Percent),
                color = !is.na(Proportion) & Proportion > .5), size = 3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2D6F85", limits = c(0, 1),
                      na.value = "#F7FBFF",
                      labels = function(x) paste0(100 * x, "%"),
                      name = "% of isolates with SNP") +
  scale_y_discrete(labels = italic_gene_labels) +
  labs(x = "panC group", y = NULL) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.title.x = element_text(face = "italic", color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"))
save_plot(fig_heat, "Fig8_top15_SNP_by_panC_group_heatmap", width = 8, height = 6, dpi = 320)


cat("\nWrote figures (PNG + SVG) and tables to:\n", fig_dir, "\n")
