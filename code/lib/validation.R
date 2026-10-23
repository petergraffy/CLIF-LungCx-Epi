# All phenotype learning for Model C is repeated inside each training fold.
fold_phenotypes <- function(sequence, train, test, k_max) {
  alphabet <- sort(unique(as.vector(sequence)))
  sq <- TraMineR::seqdef(as.data.frame(sequence), alphabet = alphabet)
  training <- sq[train, ]
  costs <- TraMineR::seqsubm(training, method = "TRATE")
  dm <- as.matrix(TraMineR::seqdist(training, method = "OM", sm = costs, indel = 1))
  hc <- stats::hclust(stats::as.dist(dm), method = "ward.D2")
  ks <- seq.int(2, min(k_max, length(train) - 1L))
  widths <- vapply(ks, function(k) mean(cluster::silhouette(stats::cutree(hc, k), stats::as.dist(dm))[, 3]), numeric(1))
  k <- ks[which.max(widths)]; label <- stats::cutree(hc, k)
  medoid <- vapply(seq_len(k), function(j) {
    ix <- which(label == j); train[ix[which.min(rowSums(dm[ix, ix, drop = FALSE]))]]
  }, integer(1))
  # refseq computes distances to training medoids only, avoiding an N x N test matrix.
  distances <- as.matrix(TraMineR::seqdist(sq[c(test, medoid), ], method = "OM", sm = costs, indel = 1,
    refseq = list(seq_along(test), length(test) + seq_along(medoid))))
  predicted <- max.col(-distances[seq_along(test), , drop = FALSE], ties.method = "first")
  list(train = factor(label, levels = seq_len(k)), test = factor(predicted, levels = seq_len(k)), k = k)
}

prediction_metrics <- function(y, p) {
  p <- pmin(1 - 1e-6, pmax(1e-6, p)); lp <- stats::qlogis(p)
  positive <- sum(y == 1); negative <- sum(y == 0)
  auc <- if (positive && negative) (sum(rank(p)[y == 1]) - positive * (positive + 1) / 2) / (positive * negative) else NA_real_
  calibration <- if (positive && negative && stats::sd(lp) > 0) suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())) else NULL
  citl <- if (positive && negative) suppressWarnings(stats::glm(y ~ 1 + offset(lp), family = stats::binomial())) else NULL
  c(auc = auc, brier = mean((y - p)^2), calibration_intercept = if (is.null(citl)) NA_real_ else unname(stats::coef(citl)[1]),
    calibration_slope = if (is.null(calibration)) NA_real_ else unname(stats::coef(calibration)[2]))
}

validate_nested <- function(d, sequence, cfg) {
  baseline <- c("age_years", "sex_category", "race_category", "admission_year", "sofa_total", "code_at_admission")
  # Non-reference duration categories avoid the compositional sum-to-window alias.
  summaries <- c("support_at_landmark", paste0("hours_", setdiff(SUPPORT_STATES, "ROOM_AIR")))
  required <- c(baseline, summaries, "patient_id", "event_type", "time_from_landmark")
  keep <- stats::complete.cases(d[required])
  x <- d[keep, ]; sequence <- sequence[keep, , drop = FALSE]
  audit <- tibble::tibble(candidate_n = nrow(d), complete_case_n = nrow(x), excluded_n = sum(!keep),
    endpoint = "in_hospital_death_by_post_landmark_horizon", horizon_days = cfg$analysis$prediction_days,
    phenotype_learning = "training_fold_only", status = "pending")
  if (nrow(x) < 30 || dplyr::n_distinct(x$patient_id) < cfg$analysis$validation_folds * 2) {
    audit$status <- "insufficient_complete_cases"; return(list(audit = audit))
  }
  x$y <- as.integer(x$event_type == "death" & x$time_from_landmark <= cfg$analysis$prediction_days)
  # All encounters have a known terminal hospital disposition/time. Hospice and
  # other live discharge are competing outcomes coded 0 for in-hospital death.
  set.seed(cfg$analysis$seed)
  ids <- sample(unique(x$patient_id)); folds <- rep(seq_len(cfg$analysis$validation_folds), length.out = length(ids))
  x$fold <- folds[match(x$patient_id, ids)]
  for (v in c("sex_category", "race_category", "code_at_admission", "support_at_landmark")) x[[v]] <- factor(x[[v]])
  predictions <- matrix(NA_real_, nrow(x), 3, dimnames = list(NULL, c("A", "B", "C")))
  fold_audit <- list()
  for (fold in seq_len(cfg$analysis$validation_folds)) {
    tr <- which(x$fold != fold); te <- which(x$fold == fold)
    result <- tryCatch({
      if (length(unique(x$y[tr])) < 2) stop("Training fold has only one outcome class")
      phen <- fold_phenotypes(sequence, tr, te, cfg$analysis$k_max)
      train <- x[tr, ]; test <- x[te, ]; train$phenotype <- phen$train; test$phenotype <- phen$test
      lists <- list(A = baseline, B = c(baseline, summaries), C = c(baseline, summaries, "phenotype"))
      pp <- matrix(NA_real_, length(te), 3)
      for (j in seq_along(lists)) {
        vars <- lists[[j]][vapply(train[lists[[j]]], function(z) length(unique(z)) > 1, logical(1))]
        for (v in vars) if (is.factor(train[[v]]) && any(!as.character(test[[v]]) %in% as.character(train[[v]]))) stop("Unseen validation level in ", v)
        fit <- stats::glm(stats::reformulate(vars, "y"), data = train, family = stats::binomial())
        if (!fit$converged || any(abs(stats::coef(fit)[is.finite(stats::coef(fit))]) > 30)) stop("Nonconvergent or separated prediction model")
        pp[, j] <- suppressWarnings(stats::predict(fit, newdata = test, type = "response"))
        if (any(!is.finite(pp[, j]))) stop("Nonfinite validation predictions")
      }
      predictions[te, ] <- pp
      "ok"
    }, error = function(e) conditionMessage(e))
    fold_audit[[fold]] <- tibble::tibble(fold, training_n = length(tr), validation_n = length(te), status = result)
  }
  if (anyNA(predictions)) {
    audit$status <- "incomplete_validation_no_performance_claim"
    return(list(audit = audit, folds = dplyr::bind_rows(fold_audit)))
  }
  audit$status <- "ok"; audit$events <- sum(x$y); audit$non_events <- sum(!x$y)
  point <- vapply(seq_len(3), function(j) prediction_metrics(x$y, predictions[, j]), numeric(4))
  metrics <- tibble::tibble(model = rep(c("A", "B", "C"), each = 4), metric = rep(rownames(point), 3), estimate = as.vector(point))
  # Patient bootstrap on paired out-of-fold predictions gives uncertainty for
  # performance and differences, conditional on this fitted CV procedure.
  boot <- list(); set.seed(cfg$analysis$seed + 1L)
  for (b in seq_len(cfg$analysis$validation_bootstrap_reps)) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    ix <- unlist(lapply(sampled, function(id) which(x$patient_id == id)), use.names = FALSE)
    vals <- vapply(seq_len(3), function(j) prediction_metrics(x$y[ix], predictions[ix, j]), numeric(4))
    boot[[b]] <- vals
  }
  if (length(boot)) {
    arr <- simplify2array(boot)
    metrics$std_error <- as.vector(apply(arr, c(1, 2), stats::sd, na.rm = TRUE))
    metrics$conf_low <- as.vector(apply(arr, c(1, 2), stats::quantile, probs = .025, na.rm = TRUE))
    metrics$conf_high <- as.vector(apply(arr, c(1, 2), stats::quantile, probs = .975, na.rm = TRUE))
  }
  differences <- dplyr::bind_rows(lapply(list(c(2, 1), c(3, 2)), function(pair) {
    delta <- if (length(boot)) vapply(boot, function(z) z[, pair[1]] - z[, pair[2]], numeric(4)) else matrix(NA_real_, 4, 1)
    tibble::tibble(comparison = paste0(c("A", "B", "C")[pair[1]], "_minus_", c("A", "B", "C")[pair[2]]),
      metric = rownames(point), estimate = point[, pair[1]] - point[, pair[2]],
      std_error = apply(delta, 1, stats::sd, na.rm = TRUE),
      conf_low = apply(delta, 1, stats::quantile, probs = .025, na.rm = TRUE),
      conf_high = apply(delta, 1, stats::quantile, probs = .975, na.rm = TRUE))
  }))
  calibration <- dplyr::bind_rows(lapply(seq_len(3), function(j) {
    bins <- pmin(10L, floor(predictions[, j] * 10) + 1L)
    tibble::tibble(bin = bins, observed = x$y, predicted = predictions[, j]) |>
      dplyr::group_by(bin) |> dplyr::summarise(model = c("A", "B", "C")[j], n = dplyr::n(),
        events = sum(observed), non_events = sum(1 - observed), mean_predicted = mean(predicted),
        observed_risk = mean(observed), .groups = "drop")
  }))
  decision <- dplyr::bind_rows(lapply(seq_len(3), function(j) dplyr::bind_rows(lapply(seq(.05, .8, .05), function(t) {
    positive <- predictions[, j] >= t
    tibble::tibble(model = c("A", "B", "C")[j], threshold = t, n = nrow(x),
      tp = sum(positive & x$y == 1), fp = sum(positive & x$y == 0),
      net_benefit = mean(positive & x$y == 1) - mean(positive & x$y == 0) * t / (1 - t),
      treat_all = mean(x$y) - mean(1 - x$y) * t / (1 - t), treat_none = 0)
  }))))
  list(audit = audit, folds = dplyr::bind_rows(fold_audit), metrics = metrics,
       differences = differences, calibration = calibration, decision = decision)
}
