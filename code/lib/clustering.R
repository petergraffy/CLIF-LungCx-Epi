sequence_fit <- function(panel, field = "joint", k_max = 8L, fixed_k = NULL) {
  wide <- panel |> dplyr::select(hospitalization_id, h, state = dplyr::all_of(field)) |>
    tidyr::pivot_wider(names_from = h, values_from = state, names_prefix = "h") |>
    dplyr::arrange(hospitalization_id)
  cols <- paste0("h", sort(unique(panel$h)))
  if (nrow(wide) < 4) stop("Too few sequences for clustering")
  sq <- TraMineR::seqdef(as.data.frame(wide[cols]), alphabet = sort(unique(unlist(wide[cols]))),
                        right = NA, left = NA, gaps = NA)
  dm <- as.matrix(TraMineR::seqdist(sq, method = "OM", indel = 1, sm = "TRATE", with.missing = FALSE))
  if (!any(dm > 0)) stop("All sequences are identical; no supported phenotype partition")
  hc <- stats::hclust(stats::as.dist(dm), method = "ward.D2")
  candidates <- seq.int(2L, min(as.integer(k_max), nrow(wide) - 1L))
  sil <- tibble::tibble(k = candidates, avg_silhouette = vapply(candidates, function(k) {
    mean(cluster::silhouette(stats::cutree(hc, k), stats::as.dist(dm))[, "sil_width"])
  }, numeric(1)))
  k <- fixed_k %||% sil$k[which.max(sil$avg_silhouette)]
  labels <- as.integer(stats::cutree(hc, k))
  sil$selected <- sil$k == k
  # Cluster IDs deliberately carry no outcome or inferred severity information.
  list(ids = wide$hospitalization_id, labels = labels, distance = dm, k = k,
       silhouette = sil, field = field, sequence = as.matrix(wide[cols]))
}

adjusted_rand <- function(a, b) {
  choose2 <- function(x) x * (x - 1) / 2
  tab <- table(a, b); n <- sum(tab)
  observed <- sum(choose2(tab)); row <- sum(choose2(rowSums(tab))); col <- sum(choose2(colSums(tab)))
  expected <- row * col / choose2(n); denom <- (row + col) / 2 - expected
  if (denom == 0) return(ifelse(identical(as.integer(factor(a)), as.integer(factor(b))), 1, NA_real_))
  (observed - expected) / denom
}

rectangular_assignment <- function(cost) {
  # Padding considers every prototype, including when n_local < n_consensus.
  n <- nrow(cost); m <- ncol(cost)
  finite <- cost[is.finite(cost)]; penalty <- if (length(finite)) max(1, finite) * 1e6 else 1e6
  cost[!is.finite(cost)] <- penalty
  if (n > m) cost <- cbind(cost, matrix(penalty, n, n - m))
  assignment <- as.integer(clue::solve_LSAP(cost))
  assignment[assignment > m | cost[cbind(seq_len(n), assignment)] >= penalty] <- NA_integer_
  assignment
}

bootstrap_stability <- function(fit, patient_id, reps, seed) {
  set.seed(seed)
  patients <- unique(patient_id)
  out <- list()
  for (b in seq_len(reps)) {
    selected <- sample(patients, length(patients), replace = TRUE)
    ix <- unlist(lapply(selected, function(p) which(patient_id == p)), use.names = FALSE)
    if (length(unique(ix)) <= fit$k) next
    sub <- fit$distance[ix, ix, drop = FALSE]
    cl <- stats::cutree(stats::hclust(stats::as.dist(sub), method = "ward.D2"), fit$k)
    medoid <- vapply(seq_len(fit$k), function(k) {
      j <- which(cl == k); ix[j[which.min(rowSums(sub[j, j, drop = FALSE]))]]
    }, integer(1))
    predicted <- max.col(-fit$distance[, medoid, drop = FALSE], ties.method = "first")
    jaccard <- outer(seq_len(fit$k), seq_len(fit$k), Vectorize(function(a, b) {
      sum(fit$labels == a & predicted == b) / sum(fit$labels == a | predicted == b)
    }))
    assignment <- rectangular_assignment(1 - jaccard)
    out[[length(out) + 1L]] <- tibble::tibble(replicate = b, local_cluster = seq_len(fit$k),
      jaccard = jaccard[cbind(seq_len(fit$k), assignment)],
      adjusted_rand = adjusted_rand(fit$labels, predicted))
  }
  dplyr::bind_rows(out)
}

cluster_signature <- function(panel, assignment, field) {
  alphabet <- if (field == "support") c(SUPPORT_STATES, TERMINAL_STATES) else
    c(as.vector(outer(SUPPORT_STATES, PHYS_STATES, paste, sep = "|")), TERMINAL_STATES)
  p <- panel |> dplyr::inner_join(assignment, by = "hospitalization_id")
  sizes <- assignment |> dplyr::count(local_cluster, name = "cluster_n")
  p |> dplyr::count(local_cluster, h, state = .data[[field]], name = "n") |>
    tidyr::complete(local_cluster, h, state = alphabet, fill = list(n = 0L)) |>
    dplyr::left_join(sizes, by = "local_cluster") |>
    dplyr::mutate(probability = n / cluster_n)
}

signature_matrix <- function(signature) {
  signature <- signature |> dplyr::mutate(feature = paste(h, state, sep = ":"))
  x <- signature |> dplyr::select(site, local_cluster, feature, probability) |>
    tidyr::pivot_wider(names_from = feature, values_from = probability) |>
    dplyr::arrange(site, local_cluster)
  list(keys = x[c("site", "local_cluster")], values = as.matrix(x[setdiff(names(x), c("site", "local_cluster"))]))
}

profile_distance <- function(x, y, min_overlap = .5) {
  out <- matrix(Inf, nrow(x), nrow(y))
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(y))) {
    ok <- is.finite(x[i, ]) & is.finite(y[j, ])
    if (mean(ok) >= min_overlap) out[i, j] <- sqrt(mean((x[i, ok] - y[j, ok])^2))
  }
  out
}

fit_consensus <- function(signature, k = NULL, min_overlap = .5, max_iter = 100L) {
  sm <- signature_matrix(signature); x <- sm$values
  sizes <- signature |> dplyr::distinct(site, local_cluster, cluster_n)
  keys <- sm$keys |> dplyr::left_join(sizes, by = c("site", "local_cluster"))
  if (anyNA(keys$cluster_n)) stop("Suppressed cluster sizes cannot be used for consensus weights")
  counts <- table(keys$site)
  k <- k %||% as.integer(stats::median(counts))
  if (k < 2 || k > nrow(x)) stop("Invalid consensus k")
  # Deterministic seed site, selected only by representation coverage and size.
  coverage <- rowMeans(is.finite(x))
  candidates <- names(counts)[counts == k]
  if (!length(candidates)) stop("No site has the proposed consensus k; set consensus_k explicitly")
  site_coverage <- vapply(candidates, function(s) mean(coverage[keys$site == s]), numeric(1))
  seed_site <- candidates[order(-site_coverage, candidates)][1]
  prototypes <- x[keys$site == seed_site, , drop = FALSE]
  old <- rep(NA_integer_, nrow(x)); converged <- FALSE
  for (it in seq_len(max_iter)) {
    costs <- profile_distance(x, prototypes, min_overlap)
    labels <- rep(NA_integer_, nrow(x))
    for (s in unique(keys$site)) {
      ix <- which(keys$site == s)
      labels[ix] <- rectangular_assignment(costs[ix, , drop = FALSE])
    }
    if (all(is.na(labels))) stop("No profiles have sufficient unsuppressed overlap for matching")
    if (identical(labels, old)) { converged <- TRUE; break }
    old <- labels
    for (j in seq_len(k)) {
      ix <- which(labels == j)
      if (!length(ix)) next
      for (c in seq_len(ncol(x))) {
        ok <- ix[is.finite(x[ix, c])]
        prototypes[j, c] <- if (length(ok)) stats::weighted.mean(x[ok, c], keys$cluster_n[ok]) else NA_real_
      }
    }
  }
  if (!converged) stop("Consensus failed to converge; inspect signatures before continuing")
  # Deterministic reference and labels depend only on respiratory inputs.
  support <- sub("^[0-9]+:", "", colnames(prototypes))
  support <- sub("\\|.*$", "", support)
  rank <- match(support, SUPPORT_STATES) - 1L
  rank[support %in% c("OTHER", "UNKNOWN", TERMINAL_STATES)] <- NA_real_
  burden <- apply(prototypes, 1, function(p) sum(p * rank, na.rm = TRUE) / sum(p[!is.na(rank)], na.rm = TRUE))
  ordering <- order(burden, seq_len(k)); new_id <- match(seq_len(k), ordering)
  distance <- rep(NA_real_, nrow(x)); matched <- which(!is.na(labels))
  distance[matched] <- costs[cbind(matched, labels[matched])]
  keys$phenotype <- ifelse(is.na(labels), NA_character_, paste0("P", new_id[labels]))
  keys$matching_distance <- distance; keys$matched <- !is.na(labels)
  prototypes <- prototypes[ordering, , drop = FALSE]
  rownames(prototypes) <- paste0("P", seq_len(k))
  list(mapping = keys, prototypes = prototypes, seed_site = seed_site, iterations = it,
       reference = "P1", k = k, min_overlap = min_overlap)
}

consensus_loso <- function(signature, fit) {
  out <- list()
  for (s in unique(signature$site)) {
    z <- tryCatch(fit_consensus(signature[signature$site != s, ], fit$k, fit$min_overlap), error = identity)
    if (inherits(z, "error")) {
      out[[s]] <- tibble::tibble(held_out_site = s, status = conditionMessage(z), agreement = NA_real_)
      next
    }
    align <- rectangular_assignment(profile_distance(z$prototypes, fit$prototypes, fit$min_overlap))
    comparison <- dplyr::inner_join(z$mapping, fit$mapping, by = c("site", "local_cluster"), suffix = c("_loo", "_full"))
    pred <- paste0("P", align[as.integer(sub("P", "", comparison$phenotype_loo))])
    ok <- !is.na(comparison$phenotype_loo) & !is.na(comparison$phenotype_full)
    held <- signature_matrix(signature[signature$site == s, ])
    held_x <- held$values[, colnames(z$prototypes), drop = FALSE]
    held_assignment <- rectangular_assignment(profile_distance(held_x, z$prototypes, fit$min_overlap))
    held_pred <- ifelse(is.na(held_assignment), NA_character_, paste0("P", align[held_assignment]))
    truth <- fit$mapping |> dplyr::filter(site == s) |> dplyr::arrange(local_cluster)
    out[[s]] <- tibble::tibble(held_out_site = s, status = "ok",
      agreement = if (any(ok)) mean(pred[ok] == comparison$phenotype_full[ok]) else NA_real_,
      held_out_agreement = mean(held_pred == truth$phenotype, na.rm = TRUE))
  }
  dplyr::bind_rows(out)
}
