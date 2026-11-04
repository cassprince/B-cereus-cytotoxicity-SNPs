# ============================================================
# Random forest regression and classification of SNP and gene data for cytotoxicity in the Bacillus cereus group.
#
# Discovery dataset: 270 isolates.
#
# Design
#   - Stratified 70/30 train/test split by panC group (seed 3).
#   - Tuning: 30 random hyperparameter sets, each scored by 10-fold
#   - panC-stratified cross-validation on the training isolates.
#   - Final models are fitted on all training isolates and evaluated once on the test isolates.
#   - Binary classification: High (cytotoxicity >= 0.70) vs Low (< 0.70), predicted High when the forest probability >= 0.50.
#   - Variable importance: permutation importance from 100 bootstrap refits of the training set, reported as 95% percentile intervals.
#
# Model
#   Categorical SNP status only (no panC group)
# ============================================================

# ============================================================
# 1. SETUP
# ============================================================

library(ranger)
library(ggplot2)

base_folder <- path.expand(paste0(
  "~/Library/CloudStorage/",
  "OneDrive-ThePennsylvaniaStateUniversity/",
  "Tyler's projects - Collaboration and External Projects/",
  "Cassidy_analysis/Random Forrest"
))

outdir <- file.path(
  base_folder,
  "FINAL_RF_panC_stratified_70_30"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_trees <- 1000
n_folds <- 10
n_tuning_evals <- 30
n_bootstraps <- 100
random_seed <- 3
train_ratio <- 0.70
cytotoxicity_cutoff <- 0.70
probability_cutoff <- 0.50

save_csv <- function(x, filename) {
  write.csv(
    x,
    file.path(outdir, filename),
    row.names = FALSE
  )
}

read_data <- function(filename) {
  read.csv(
    file.path(base_folder, filename),
    check.names = FALSE,
    na.strings = c("NA", ""),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 2. LOAD AND CHECK DATA
# ============================================================

pres <- read_data("RF_analysis/snps_all_signif_pres_abs_102126_no_panC_sig.csv")
summed <- read_data("RF_analysis/snps_all_signif_summed_102126_no_panC_sig.csv")
genes <- read_data("df_all_genes_090126.csv")

names(genes)[names(genes) == "acc"] <- "Isolate"

metadata_cols <- c("Isolate", "cytotoxicity", "panC_group")

snp_cols <- setdiff(names(pres), metadata_cols)

copy_genes <- c(
  "nheA", "nheB", "nheC",
  "hblA", "hblB", "hblC", "hblD"
)

presence_genes <- c(copy_genes, "cytK1", "cytK2")

stopifnot(
  all(metadata_cols %in% names(pres)),
  all(metadata_cols %in% names(summed)),
  all(metadata_cols %in% names(genes)),
  !anyDuplicated(pres$Isolate),
  !anyDuplicated(summed$Isolate),
  !anyDuplicated(genes$Isolate),
  setequal(pres$Isolate, summed$Isolate),
  setequal(pres$Isolate, genes$Isolate),
  length(snp_cols) == 214,
  all(snp_cols %in% names(summed)),
  all(presence_genes %in% names(genes))
)

# Align all three tables to the same isolate order.
summed <- summed[match(pres$Isolate, summed$Isolate), ]
genes <- genes[match(pres$Isolate, genes$Isolate), ]

stopifnot(
  isTRUE(all.equal(as.numeric(pres$cytotoxicity),
                   as.numeric(summed$cytotoxicity))),
  isTRUE(all.equal(as.numeric(pres$cytotoxicity),
                   as.numeric(genes$cytotoxicity))),
  identical(as.character(pres$panC_group),
            as.character(summed$panC_group)),
  identical(as.character(pres$panC_group),
            as.character(genes$panC_group)),
  !anyNA(pres$cytotoxicity),
  !anyNA(pres$panC_group)
)

for (gene in presence_genes) {

  x <- as.numeric(genes[[gene]])

  stopifnot(
    !anyNA(x),
    all(x >= 0),
    all(x == floor(x))
  )

  genes[[gene]] <- x
}

# Binary SNP status must agree with the SNP variant-copy counts.
for (snp in snp_cols) {

  a <- as.numeric(pres[[snp]])
  b <- as.numeric(summed[[snp]])

  valid <- !is.na(a)

  stopifnot(
    identical(is.na(a), is.na(b)),
    all(a[valid] %in% c(0, 1)),
    all(b[valid] >= 0),
    all(b[valid] == floor(b[valid])),
    all(a[valid] == as.numeric(b[valid] > 0))
  )
}

panC_levels <- sort(unique(as.character(pres$panC_group)))


# ============================================================
# 3. PREDICTOR MATRICES
# ============================================================

base <- data.frame(
  cytotoxicity = as.numeric(pres$cytotoxicity),
  panC_group = factor(pres$panC_group, levels = panC_levels),
  row.names = pres$Isolate
)

# SNP status encoding:
#   0 = gene or region absent
#   1 = reference allele
#   2 = variant allele
categorical_snps <- function(df, source, snps, keep = NULL) {

  if (is.null(keep)) {
    keep <- rep(TRUE, nrow(source))
  }

  for (snp in snps) {

    x <- as.numeric(source[[snp]][keep])

    x[is.na(x)] <- -1

    df[[snp]] <- factor(x + 1, levels = c(0, 1, 2))
  }

  df
}

# Categorical SNP status only (no panC group)
model_data <- list(
  SNP_status_no_panC = categorical_snps(
    base[, "cytotoxicity", drop = FALSE],
    pres,
    snp_cols
  )
)

for (nm in names(model_data)) {

  df <- model_data[[nm]]

  names(df) <- make.names(names(df), unique = TRUE)

  stopifnot(
    !anyNA(df),
    nrow(df) >= 10
  )

  model_data[[nm]] <- df
}


# ============================================================
# 4. STRATIFIED 70/30 TRAIN/TEST SPLIT
# ============================================================

set.seed(random_seed)

all_ids <- rownames(base)

train_ids <- unlist(lapply(
  split(all_ids, base$panC_group),
  function(ids) {
    sample(
      ids,
      size = max(1L, min(length(ids), round(length(ids) * train_ratio)))
    )
  }
), use.names = FALSE)

test_ids <- setdiff(all_ids, train_ids)

stopifnot(length(train_ids) > 0, length(test_ids) > 0)

save_csv(
  data.frame(
    Isolate = all_ids,
    panC_group = as.character(base$panC_group),
    Set = ifelse(all_ids %in% train_ids, "train", "test")
  ),
  "shared_split_assignment.csv"
)

panC_counts <- as.data.frame(table(base$panC_group))
names(panC_counts) <- c("panC_group", "N_discovery")
save_csv(panC_counts, "full_discovery_panC_counts.csv")
print(panC_counts)


# ============================================================
# 5. PANC-STRATIFIED CROSS-VALIDATION FOLDS
# ============================================================


stratified_folds <- function(ids, groups, k = 10, seed = 3) {

  set.seed(seed)

  fold_id <- rep(NA_integer_, length(ids))
  names(fold_id) <- ids

  for (g in sort(unique(as.character(groups)))) {

    group_ids <- ids[as.character(groups) == g]

    n <- length(group_ids)

    if (n == 1) {
      next
    }

    shuffled <- sample(group_ids, size = n, replace = FALSE)

    assignments <- rep(seq_len(k), length.out = n)

    # Rotate the fold allocation for each panC group.
    offset <- sample.int(k, 1) - 1

    assignments <- ((assignments + offset - 1) %% k) + 1

    fold_id[shuffled] <- assignments
  }

  fold_id
}

train_groups <- as.character(base[train_ids, "panC_group"])

global_fold_ids <- stratified_folds(
  ids = train_ids,
  groups = train_groups,
  k = n_folds,
  seed = random_seed
)

save_csv(
  data.frame(
    Isolate = names(global_fold_ids),
    panC_group = as.character(base[names(global_fold_ids), "panC_group"]),
    Fold = as.integer(global_fold_ids)
  ),
  "shared_panC_stratified_CV_folds.csv"
)


# ============================================================
# 6. RANDOM FOREST HELPERS
# ============================================================

# Predictors that vary within the supplied data.
informative_predictors <- function(df, exclude = "cytotoxicity") {

  predictors <- setdiff(names(df), exclude)

  predictors[vapply(
    df[predictors],
    function(x) length(unique(x)) > 1,
    logical(1)
  )]
}

# Fit a random forest.
fit_forest <- function(
    train_data,
    mtry,
    min_node,
    seed,
    outcome = "cytotoxicity",
    probability = FALSE
) {

  predictors <- informative_predictors(
    train_data,
    exclude = unique(c("cytotoxicity", outcome))
  )

  if (length(predictors) == 0) {
    stop("No variable predictors remain.")
  }

  ranger(
    dependent.variable.name = outcome,
    data = train_data[, c(outcome, predictors), drop = FALSE],
    probability = probability,
    mtry = min(mtry, length(predictors)),
    min.node.size = min_node,
    num.trees = n_trees,
    importance = "permutation",
    respect.unordered.factors = "partition",
    seed = seed
  )
}

# Predictions for new data.
predict_forest <- function(
    forest,
    new_data,
    positive = "High"
) {

  required <- forest$forest$independent.variable.names

  missing <- setdiff(required, names(new_data))

  if (length(missing) > 0) {
    stop(paste(
      "Prediction data missing:",
      paste(missing, collapse = ", ")
    ))
  }

  predictions <- predict(
    forest,
    data = new_data[, required, drop = FALSE]
  )$predictions

  if (is.matrix(predictions)) {
    return(as.numeric(predictions[, positive]))
  }

  as.numeric(predictions)
}

regression_metrics <- function(observed, predicted) {

  residual <- observed - predicted

  denominator <- sum((observed - mean(observed))^2)

  data.frame(
    MSE = mean(residual^2),
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    R2 = if (denominator > 0) {
      1 - sum(residual^2) / denominator
    } else {
      NA_real_
    },
    Bias = mean(predicted - observed)
  )
}

# Confusion-matrix metrics for a binary outcome.
# actual: 1 = High, 0 = Low. probability: predicted probability of High.
classification_metrics <- function(actual, probability, cutoff = 0.50) {

  actual <- as.integer(actual)
  predicted <- as.integer(probability >= cutoff)

  TP <- sum(predicted == 1L & actual == 1L)
  TN <- sum(predicted == 0L & actual == 0L)
  FP <- sum(predicted == 1L & actual == 0L)
  FN <- sum(predicted == 0L & actual == 1L)
  n <- length(actual)

  accuracy <- (TP + TN) / n
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

  f1 <- if (is.na(precision) || is.na(sensitivity) ||
            (precision + sensitivity) == 0) {
    NA_real_
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }

  # Cohen's kappa
  expected <- ((TP + FP) * (TP + FN) + (TN + FN) * (TN + FP)) / n^2
  kappa <- if (expected == 1) {
    NA_real_
  } else {
    (accuracy - expected) / (1 - expected)
  }

  # AUC as the Mann-Whitney probability (ties count one half)
  pos <- probability[actual == 1L]
  neg <- probability[actual == 0L]
  auc <- if (length(pos) > 0 && length(neg) > 0) {
    mean(outer(pos, neg, ">") + 0.5 * outer(pos, neg, "=="))
  } else {
    NA_real_
  }

  data.frame(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Kappa = kappa,
    AUC = auc,
    N_high = sum(actual == 1L),
    N_low = sum(actual == 0L)
  )
}


# ============================================================
# 7. HYPERPARAMETER TUNING (TRAINING SET ONLY)
# ============================================================

tune_forest <- function(label, train_data, fold_ids) {

  p <- length(informative_predictors(train_data))

  stopifnot(p > 0)

  # Random sample of mtry and minimum node size combinations.
  grid <- expand.grid(
    mtry = seq_len(p),
    min.node.size = 3:10
  )

  set.seed(random_seed)

  selected <- sample(
    seq_len(nrow(grid)),
    n_tuning_evals,
    replace = nrow(grid) < n_tuning_evals
  )

  candidates <- grid[selected, , drop = FALSE]

  tuning_results <- vector("list", nrow(candidates))

  model_ids <- rownames(train_data)

  available_folds <- sort(unique(fold_ids[!is.na(fold_ids)]))

  if (length(available_folds) < 2) {
    stop("Fewer than two usable CV folds.")
  }

  for (evaluation in seq_len(nrow(candidates))) {

    mtry_value <- candidates$mtry[evaluation]
    node_value <- candidates$min.node.size[evaluation]

    fold_errors <- rep(NA_real_, length(available_folds))

    for (j in seq_along(available_folds)) {

      fold <- available_folds[j]

      held_out_ids <- model_ids[
        !is.na(fold_ids[model_ids]) & fold_ids[model_ids] == fold
      ]

      fit_ids <- setdiff(model_ids, held_out_ids)

      if (length(held_out_ids) == 0 || length(fit_ids) < 5) {
        next
      }

      forest <- fit_forest(
        train_data = train_data[fit_ids, , drop = FALSE],
        mtry = mtry_value,
        min_node = node_value,
        seed = random_seed + evaluation * 100 + j
      )

      predicted <- predict_forest(
        forest,
        train_data[held_out_ids, , drop = FALSE]
      )

      observed <- train_data[held_out_ids, "cytotoxicity"]

      fold_errors[j] <- mean((observed - predicted)^2)
    }

    valid_errors <- fold_errors[is.finite(fold_errors)]

    cv_mse <- if (length(valid_errors) == length(available_folds)) {
      mean(valid_errors)
    } else {
      Inf
    }

    tuning_results[[evaluation]] <- data.frame(
      Evaluation = evaluation,
      mtry = mtry_value,
      min.node.size = node_value,
      CV_MSE = cv_mse,
      CV_RMSE = sqrt(cv_mse),
      N_folds = length(valid_errors)
    )

    cat(
      label,
      "| Evaluation", evaluation,
      "/", nrow(candidates),
      "| CV MSE =", round(cv_mse, 5),
      "\n"
    )
  }

  tuning_table <- do.call(rbind, tuning_results)

  save_csv(tuning_table, paste0(label, "_tuning.csv"))

  if (!any(is.finite(tuning_table$CV_MSE))) {
    stop("No successful tuning evaluation for ", label)
  }

  best <- tuning_table[which.min(tuning_table$CV_MSE), , drop = FALSE]

  list(
    mtry = best$mtry,
    min.node.size = best$min.node.size,
    tuning = tuning_table
  )
}


# ============================================================
# 8. BOOTSTRAP VARIABLE IMPORTANCE
# ============================================================

bootstrap_importance <- function(label, train_data, best_mtry, best_node) {

  predictors <- informative_predictors(train_data)

  imp_matrix <- matrix(
    NA_real_,
    nrow = n_bootstraps,
    ncol = length(predictors),
    dimnames = list(NULL, predictors)
  )

  set.seed(random_seed)

  bootstrap_ids <- lapply(
    seq_len(n_bootstraps),
    function(i) {
      sample.int(nrow(train_data), size = nrow(train_data), replace = TRUE)
    }
  )

  for (i in seq_len(n_bootstraps)) {

    boot <- train_data[bootstrap_ids[[i]], , drop = FALSE]

    # Predictors can become constant in a resample. Skip the refit if none remain.
    if (length(informative_predictors(boot)) == 0) {
      next
    }

    forest <- fit_forest(
      train_data = boot,
      mtry = best_mtry,
      min_node = best_node,
      seed = random_seed + i
    )

    vi <- forest$variable.importance

    imp_matrix[i, names(vi)] <- as.numeric(vi)

    if (i %% 10 == 0) {
      cat(label, "| Bootstrap", i, "/", n_bootstraps, "\n")
    }
  }

  quantile_or_na <- function(x, prob) {

    if (all(is.na(x))) {
      return(NA_real_)
    }

    as.numeric(quantile(x, probs = prob, na.rm = TRUE))
  }

  importance <- data.frame(
    Variable = colnames(imp_matrix),
    Mean = colMeans(imp_matrix, na.rm = TRUE),
    Lower = apply(imp_matrix, 2, quantile_or_na, prob = 0.025),
    Upper = apply(imp_matrix, 2, quantile_or_na, prob = 0.975),
    N_refits = colSums(!is.na(imp_matrix))
  )

  importance <- importance[order(-importance$Mean), ]

  save_csv(importance, paste0(label, "_bootstrap_importance.csv"))

  write.csv(
    imp_matrix,
    file.path(outdir, paste0(label, "_bootstrap_refits.csv")),
    row.names = FALSE
  )

  importance
}


# ============================================================
# 9. IMPORTANCE PLOTS
# ============================================================

italic_gene_labels <- function(x) {
  gsub(
    "^(hbl[A-D]?|nhe[A-C]?|cytK[12]?)(?=_|$)",
    "*\\1*",
    x,
    perl = TRUE
  )
}

# Dot plot of the top 30 predictors with 95% bootstrap intervals.
plot_importance <- function(label, importance) {

  plot_data <- head(importance, 15)

  plot_data <- plot_data[
    is.finite(plot_data$Mean) &
      is.finite(plot_data$Lower) &
      is.finite(plot_data$Upper),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    warning("No importance estimates for ", label)
    return(invisible(NULL))
  }

  plot_data$Variable <- factor(
    plot_data$Variable,
    levels = rev(plot_data$Variable)
  )

  p <- ggplot(plot_data, aes(x = Mean, y = Variable)) +

    geom_segment(
      aes(x = Lower, xend = Upper, yend = Variable),
      color = "#648ED1",
      linewidth = 0.6
    ) +

    geom_point(
      color = "#648ED1",
      size = 2.5
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.4
    ) +

    scale_y_discrete(labels = italic_gene_labels) +

    labs(
      x = "Permutation variable importance",
      y = "Predictor"
    ) +

    theme_classic(base_size = 12) +

    theme(
      axis.text.x = element_text(color = "black"),
      axis.text.y = ggtext::element_markdown(color = "black"),
      axis.title = element_text(face = "bold")
    )

  height <- max(5, nrow(plot_data) * 0.28 + 1.5)

  ggsave(
    file.path(outdir, paste0(label, "_importance.png")),
    plot = p,
    width = 9,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(outdir, paste0(label, "_importance.svg")),
    plot = p,
    width = 9,
    height = height,
    bg = "white"
  )

  invisible(p)
}


# ============================================================
# 10. TUNE, FIT AND EVALUATE ONE MODEL
# ============================================================

train_rf <- function(label, data) {

  cat("\n====================================\n",
      "Running:", label,
      "\n====================================\n")

  isolate_ids <- rownames(data)
  stopifnot(length(isolate_ids) >= 10)

  # Tuning and the final fit use training isolates only.
  # Test isolates are used once, for the final evaluation.
  train_data <- data[intersect(isolate_ids, train_ids), , drop = FALSE]
  test_data <- data[intersect(isolate_ids, test_ids), , drop = FALSE]
  stopifnot(nrow(train_data) >= 10, nrow(test_data) >= 1)

  # Drop predictors that are constant in the training set. Folds that
  # make a predictor constant are handled inside the CV fits.
  predictors <- informative_predictors(train_data)
  stopifnot(length(predictors) > 0)
  train_data <- train_data[, c("cytotoxicity", predictors), drop = FALSE]
  test_data <- test_data[, c("cytotoxicity", predictors), drop = FALSE]

  # Shared fold IDs restricted to this model's training isolates
  fold_ids <- global_fold_ids[rownames(train_data)]
  save_csv(
    data.frame(
      Isolate = rownames(train_data),
      panC_group = as.character(base[rownames(train_data), "panC_group"]),
      Fold = as.integer(fold_ids)
    ),
    paste0(label, "_panC_stratified_CV_folds.csv")
  )

  # Tune on the training folds.
  best <- tune_forest(
    label = label,
    train_data = train_data,
    fold_ids = fold_ids
  )

  cat("Best mtry:", best$mtry,
      "\nBest minimum node size:", best$min.node.size, "\n")

  # Final regression model, fitted on all training isolates
  final <- fit_forest(
    train_data = train_data,
    mtry = best$mtry,
    min_node = best$min.node.size,
    seed = random_seed
  )

  saveRDS(final, file.path(outdir, paste0(label, "_final.rds")))

  vi <- final$variable.importance
  final_importance <- data.frame(
    Variable = names(vi),
    Importance = as.numeric(vi)
  )
  final_importance <- final_importance[
    order(-final_importance$Importance), , drop = FALSE
  ]
  save_csv(final_importance, paste0(label, "_final_importance.csv"))

  importance <- bootstrap_importance(
    label = label,
    train_data = train_data,
    best_mtry = best$mtry,
    best_node = best$min.node.size
  )
  plot_importance(label, importance)

  # Regression performance on the held-out test isolates
  test_pred <- predict_forest(final, test_data)

  test_regression <- cbind(
    data.frame(
      Model = label,
      N_train = nrow(train_data),
      N_test = nrow(test_data)
    ),
    regression_metrics(test_data$cytotoxicity, test_pred)
  )

  save_csv(
    data.frame(
      Isolate = rownames(test_data),
      Observed = test_data$cytotoxicity,
      Predicted = test_pred,
      Residual = test_data$cytotoxicity - test_pred
    ),
    paste0(label, "_test_predictions.csv")
  )

  # Binary classification: High = cytotoxicity >= 0.70, Low = below 0.70.
  # Reuses the regression-selected mtry and minimum node size.
  class_train <- train_data
  class_train$cytotoxicity_class <- factor(
    ifelse(class_train$cytotoxicity >= cytotoxicity_cutoff, "High", "Low"),
    levels = c("Low", "High")
  )

  class_test <- test_data
  class_test$cytotoxicity_class <- factor(
    ifelse(class_test$cytotoxicity >= cytotoxicity_cutoff, "High", "Low"),
    levels = c("Low", "High")
  )

  final_class <- fit_forest(
    train_data = class_train,
    mtry = best$mtry,
    min_node = best$min.node.size,
    seed = random_seed + 3000,
    outcome = "cytotoxicity_class",
    probability = TRUE
  )

  probability_high <- predict_forest(final_class, class_test, positive = "High")

  test_classification <- cbind(
    data.frame(
      Model = label,
      N_train = nrow(train_data),
      N_test = nrow(test_data),
      Outcome_cutoff = cytotoxicity_cutoff,
      Probability_cutoff = probability_cutoff
    ),
    classification_metrics(
      as.integer(class_test$cytotoxicity_class == "High"),
      probability_high,
      cutoff = probability_cutoff
    )
  )

  saveRDS(
    final_class,
    file.path(outdir, paste0(label, "_classification_final.rds"))
  )

  save_csv(
    data.frame(
      Isolate = rownames(test_data),
      Observed_cytotoxicity = test_data$cytotoxicity,
      Observed_class = as.character(class_test$cytotoxicity_class),
      Probability_high = probability_high,
      Predicted_class = ifelse(probability_high >= probability_cutoff, "High", "Low")
    ),
    paste0(label, "_classification_test_predictions.csv")
  )

  # Metadata needed to reuse the model, including factor levels
  metadata_out <- list(
    model_name = label,
    predictors = predictors,
    factor_levels = lapply(
      train_data[vapply(train_data, is.factor, logical(1))],
      levels
    ),
    mtry = best$mtry,
    min.node.size = best$min.node.size,
    num.trees = n_trees,
    n_folds = n_folds,
    n_tuning_evals = n_tuning_evals,
    n_bootstraps = n_bootstraps,
    stratification = "panC_group",
    training_ids = rownames(train_data),
    test_ids = rownames(test_data),
    n_train = nrow(train_data),
    n_test = nrow(test_data),
    snp_names = snp_cols,
    copy_genes = copy_genes,
    presence_genes = presence_genes,
    training_panC_levels = panC_levels

  )
  saveRDS(metadata_out, file.path(outdir, paste0(label, "_metadata.rds")))

  summary_row <- data.frame(
    Model = label,
    N_isolates = nrow(data),
    N_training = nrow(train_data),
    N_test = nrow(test_data),
    N_predictors = length(predictors),
    mtry = best$mtry,
    min.node.size = best$min.node.size,
    num.trees = n_trees,
    Best_tuning_CV_MSE = min(best$tuning$CV_MSE, na.rm = TRUE)
  )

  list(
    importance = importance,
    summary_row = summary_row,
    test_regression = test_regression,
    test_classification = test_classification
  )
}


# ============================================================
# 11. RUN ALL MODELS
# ============================================================

results <- lapply(
  names(model_data),
  function(nm) {
    train_rf(nm, model_data[[nm]])
  }
)

names(results) <- names(model_data)


# ============================================================
# 12. EXPORT RESULTS
# ============================================================

model_summary <- do.call(
  rbind,
  lapply(results, function(x) x$summary_row)
)

all_importance <- do.call(
  rbind,
  lapply(
    names(results),
    function(nm) {
      data.frame(
        Model = nm,
        results[[nm]]$importance
      )
    }
  )
)

test_regression_all <- do.call(
  rbind,
  lapply(results, function(x) x$test_regression)
)

test_classification_all <- do.call(
  rbind,
  lapply(results, function(x) x$test_classification)
)

save_csv(model_summary, "RF_SNP_status_no_panC_summary.csv")
save_csv(all_importance, "RF_SNP_status_no_panC_variable_importance.csv")
save_csv(test_regression_all, "RF_SNP_status_no_panC_test_performance.csv")
save_csv(test_classification_all, "RF_SNP_status_no_panC_classification_metrics.csv")

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "session_info.txt")
)

print(model_summary)

cat(
  "\nAnalysis complete.\n",
  "Results saved to:\n",
  outdir,
  "\n"
)
