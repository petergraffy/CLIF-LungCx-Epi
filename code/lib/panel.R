SUPPORT_STATES <- c("ROOM_AIR", "LOW_O2", "HFNC", "CPAP", "NIPPV", "IMV", "OTHER", "UNKNOWN")
PHYS_STATES <- as.vector(outer(c("H0", "H1", "HU"), c("C0", "C1", "CU"), paste, sep = "_"))
TERMINAL_STATES <- c("DEATH", "HOSPICE", "LIVE_DISCHARGE")

support_state <- function(x) {
  x <- toupper(trimws(gsub("[_-]", " ", as.character(x))))
  dplyr::case_when(
    x %in% c("ROOM AIR", "RA") ~ "ROOM_AIR",
    x %in% c("NASAL CANNULA", "FACE MASK", "TRACH COLLAR", "LOW FLOW NC") ~ "LOW_O2",
    x %in% c("HIGH FLOW NC", "HIGHFLOW NC", "HIGH FLOW NASAL CANNULA", "HFNC") ~ "HFNC",
    x == "CPAP" ~ "CPAP", x %in% c("NIPPV", "NIV", "BIPAP") ~ "NIPPV",
    x == "IMV" ~ "IMV", x == "OTHER" ~ "OTHER", TRUE ~ "UNKNOWN")
}
standard_fio2 <- function(x) {
  x <- num(x)
  percent <- which(!is.na(x) & x > 1 & x <= 100)
  x[percent] <- x[percent] / 100
  x[which(!is.na(x) & (x < .21 | x > 1))] <- NA_real_
  x
}

# Nearest contemporaneous measurement; ties prefer the earlier observation.
nearest_value <- function(target, source_time, values, tolerance) {
  ans <- rep(NA_real_, length(target))
  if (!length(source_time)) return(ans)
  ok <- !is.na(source_time) & !is.na(values)
  source_time <- as.numeric(source_time[ok]); values <- values[ok]
  if (!length(source_time)) return(ans)
  ord <- order(source_time); source_time <- source_time[ord]; values <- values[ord]
  j <- findInterval(as.numeric(target), source_time)
  left <- pmax(1L, j); right <- pmin(length(source_time), j + 1L)
  choose <- ifelse(abs(as.numeric(target) - source_time[left]) <= abs(as.numeric(target) - source_time[right]), left, right)
  good <- !is.na(target) & abs(as.numeric(target) - source_time[choose]) <= tolerance * 3600
  ans[good] <- values[choose[good]]
  ans
}

physiology_observations <- function(rs, vit, labs, cfg, estimate_low_flow = FALSE) {
  rs <- with_columns(rs, c("fio2_set", "lpm_set", "flow_rate"))
  fio <- standard_fio2(rs$fio2_set)
  state <- support_state(rs$device_category)
  fio[is.na(fio) & state == "ROOM_AIR"] <- .21
  if (estimate_low_flow) {
    flow <- num(rs$lpm_set)
    # Legacy synthetic/source extracts may use flow_rate; CLIF oxygen flow is lpm_set.
    flow[is.na(flow)] <- num(rs$flow_rate)[is.na(flow)]
    use <- is.na(fio) & state == "LOW_O2" & !is.na(flow) & flow >= 1 & flow <= 6
    # Deliberately sensitivity-only: delivered FiO2 varies with patient factors.
    fio[use] <- pmin(.45, .21 + .04 * flow[use])
  }
  z <- data.frame(time = as.POSIXct(character(), tz = "UTC"), type = character(), positive = logical(), pf = numeric(), ventilated = logical(), available_time = as.POSIXct(character(), tz = "UTC"))
  add <- function(time, type, positive, pf = rep(NA_real_, length(time)), ventilated = rep(FALSE, length(time)), available_time = time) {
    data.frame(time, type = rep(type, length(time)), positive, pf, ventilated, available_time)
  }
  po <- labs[tolower(labs$lab_category) == "po2_arterial", ]
  po$value <- num(po$lab_value_numeric)
  po <- po[!is.na(po$value) & po$value >= 20 & po$value <= 800, ]
  fi <- nearest_value(po$lab_result_dttm, rs$recorded_dttm, fio, cfg$analysis$pairing_hours)
  fi_time <- utc(nearest_value(po$lab_result_dttm, rs$recorded_dttm,
    ifelse(is.na(fio), NA_real_, as.numeric(rs$recorded_dttm)), cfg$analysis$pairing_hours))
  ventilation <- nearest_value(po$lab_result_dttm, rs$recorded_dttm, as.numeric(state %in% c("IMV", "NIPPV", "CPAP")), cfg$analysis$pairing_hours) == 1
  valid <- !is.na(fi)
  z <- rbind(z, add(po$lab_result_dttm[valid], "hypox", po$value[valid] / fi[valid] <= 300,
                     po$value[valid] / fi[valid], ventilation[valid],
                     pmax(po$lab_result_dttm[valid], fi_time[valid])))
  sp <- vit[tolower(vit$vital_category) == "spo2", ]
  sp$value <- num(sp$vital_value)
  sp <- sp[!is.na(sp$value) & sp$value >= 50 & sp$value <= 100, ]
  fi <- nearest_value(sp$recorded_dttm, rs$recorded_dttm, fio, cfg$analysis$pairing_hours)
  fi_time <- utc(nearest_value(sp$recorded_dttm, rs$recorded_dttm,
    ifelse(is.na(fio), NA_real_, as.numeric(rs$recorded_dttm)), cfg$analysis$pairing_hours))
  valid <- !is.na(fi) & (fi == .21 | (estimate_low_flow & sp$value <= 97))
  positive <- ifelse(fi == .21, sp$value < 90, sp$value / fi <= cfg$analysis$sf_threshold)
  z <- rbind(z, add(sp$recorded_dttm[valid], "hypox", positive[valid],
    available_time = pmax(sp$recorded_dttm[valid], fi_time[valid])))
  co <- labs[tolower(labs$lab_category) == "pco2_arterial", ]
  ph <- labs[tolower(labs$lab_category) == "ph_arterial", ]
  co$value <- num(co$lab_value_numeric); ph$value <- num(ph$lab_value_numeric)
  co <- co[!is.na(co$value) & co$value >= 10 & co$value <= 200, ]
  ph <- ph[!is.na(ph$value) & ph$value >= 6.5 & ph$value <= 8, ]
  pp <- nearest_value(co$lab_result_dttm, ph$lab_result_dttm, ph$value, cfg$analysis$hyper_pairing_hours)
  ph_time <- utc(nearest_value(co$lab_result_dttm, ph$lab_result_dttm,
    as.numeric(ph$lab_result_dttm), cfg$analysis$hyper_pairing_hours))
  valid <- !is.na(pp)
  z <- rbind(z, add(co$lab_result_dttm[valid], "hyper", co$value[valid] >= 45 & pp[valid] < 7.35,
    available_time = pmax(co$lab_result_dttm[valid], ph_time[valid])))
  z[order(z$time), ]
}

build_panel <- function(d, tables, window, cfg, absorbing = FALSE, estimate_low_flow = FALSE, carry = NULL) {
  carry <- carry %||% cfg$analysis$physiology_carry_hours
  rs_by <- split(tables$respiratory_support, tables$respiratory_support$hospitalization_id)
  vit_by <- split(tables$vitals, tables$vitals$hospitalization_id)
  lab_by <- split(tables$labs, tables$labs$hospitalization_id)
  out <- vector("list", nrow(d))
  for (i in seq_len(nrow(d))) {
    id <- d$hospitalization_id[i]; t0 <- d$t0[i]; end <- t0 + window * 3600
    if (d$valid_outcome[i]) end <- min(end, d$event_time[i])
    select_window <- function(x, time_col, empty, pre_hours = 0) {
      if (is.null(x)) x <- empty[0, ]
      x[which(!is.na(x[[time_col]]) & x[[time_col]] >= max(d$admission_dttm[i], t0 - pre_hours * 3600) & x[[time_col]] < end), ]
    }
    r <- select_window(rs_by[[id]], "recorded_dttm", tables$respiratory_support, cfg$analysis$support_carry_hours)
    r <- r[order(r$recorded_dttm), ]
    v <- select_window(vit_by[[id]], "recorded_dttm", tables$vitals)
    l <- select_window(lab_by[[id]], "lab_result_dttm", tables$labs)
    obs <- physiology_observations(r, v, l, cfg, estimate_low_flow)
    states <- support_state(r$device_category)
    ans <- tibble::tibble(hospitalization_id = id, h = seq_len(window) - 1L,
      support = "UNKNOWN", hypox = "HU", hyper = "CU", support_observed = FALSE,
      support_carried = FALSE, terminal = FALSE)
    rh <- floor(hours_between(r$recorded_dttm, t0))
    for (j in seq_len(window)) {
      h <- j - 1L; start <- t0 + h * 3600
      if (absorbing && d$valid_outcome[i] && floor(hours_between(d$event_time[i], t0)) <= h) {
        ans$support[j] <- switch(d$event_type[i], death = "DEATH", hospice = "HOSPICE", live_discharge = "LIVE_DISCHARGE")
        ans$terminal[j] <- TRUE
        next
      }
      in_hour <- which(rh == h & states != "UNKNOWN")
      if (length(in_hour)) {
        # Highest support documented within [h,h+1); ties use latest record.
        choice <- in_hour[order(match(states[in_hour], SUPPORT_STATES), r$recorded_dttm[in_hour], decreasing = TRUE)][1]
        ans$support[j] <- states[choice]; ans$support_observed[j] <- TRUE
      } else {
        before <- which(r$recorded_dttm <= start & states != "UNKNOWN" &
                          hours_between(start, r$recorded_dttm) <= cfg$analysis$support_carry_hours)
        if (length(before)) {
          ans$support[j] <- states[tail(before, 1)]; ans$support_carried[j] <- TRUE
        }
      }
      # Carry recent evidence of either sign. Unknown is never a normal result.
      for (type in c("hypox", "hyper")) {
        z <- obs[obs$type == type & obs$time < start + 3600 & obs$time >= start - carry * 3600, ]
        if (nrow(z)) ans[[type]][j] <- paste0(ifelse(type == "hypox", "H", "C"), as.integer(any(z$positive)))
      }
    }
    ans$physiology <- paste(ans$hypox, ans$hyper, sep = "_")
    ans$joint <- ifelse(ans$terminal, ans$support, paste(ans$support, ans$physiology, sep = "|"))
    out[[i]] <- ans
  }
  dplyr::bind_rows(out)
}

support_summaries <- function(panel) {
  base <- panel |> dplyr::group_by(hospitalization_id) |>
    dplyr::summarise(support_coverage = mean(support != "UNKNOWN" & !terminal),
      observed_hours = sum(support_observed), carried_hours = sum(support_carried),
      physiology_unknown_hours = sum(hypox == "HU" | hyper == "CU"),
      support_at_landmark = support[which.max(h)],
      any_imv = as.integer(any(support == "IMV")),
      high_support = any(support %in% c("HFNC", "CPAP", "NIPPV", "IMV")),
      measured_arf = any(hypox == "H1" | hyper == "C1"), .groups = "drop")
  counts <- panel |> dplyr::count(hospitalization_id, support, name = "hours") |>
    tidyr::complete(hospitalization_id, support = SUPPORT_STATES, fill = list(hours = 0L)) |>
    tidyr::pivot_wider(names_from = support, values_from = hours, names_prefix = "hours_")
  dplyr::left_join(base, counts, by = "hospitalization_id")
}
