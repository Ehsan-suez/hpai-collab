# ============================================================
# HPAI Challenge — Phase 1
# Q2: Predict temporal + spatial evolution over next 4 weeks
#
# Uses YOUR exact SEIDR model structure:
#  - Spatial (exp kernel) + broiler movements + HRZ + env (CLC water/wetland proximity)
#  - Activity gating for transmission only (as in your code)
#  - Preventive culls from prev_culls.csv (as in your code)
#
# OUTPUTS (written to: DATA_DIR/Q2_outputs):
#  1) Raw simulated trajectories (daily) across n_sims:
#       - Q2_raw_simulated_trajectories_daily.csv
#         columns: sim, date, newE,newI,newD,newR, S,E,I,D,R
#  2) Daily forecast summaries (median + 50/90% intervals):
#       - Q2_summary_daily_quantiles_overall.csv
#  3) Spatial risk at end of 4 weeks (county-level):
#       - Q2_spatial_risk_by_county_end_of_4weeks.csv
#  4) Spatial risk at end of 4 weeks (county x production):
#       - Q2_spatial_risk_by_county_and_production_end_of_4weeks.csv
#  5) 4-week totals by production (new infections / new detections; with q05/q50/q95):
#       - Q2_summary_4week_totals_by_production.csv
#  6) Figure:
#       - Q2_forecast_newD_overall.png
#
# Notes:
#  - "Management remains as it is today" implemented by keeping:
#       activity windows, movement data, preventive culls schedule unchanged.
#  - Two-stage approach:
#       (A) Assimilation to last observed confirmation date (force_observed=TRUE)
#       (B) Forecast 28 days forward from that state (force_observed=FALSE)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(sf)
  library(RANN)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# 0) Helpers
# -----------------------------
stop_if_missing <- function(df, cols, df_name="data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) stop(df_name, " missing columns: ", paste(missing, collapse=", "))
}

pick_existing_file <- function(candidates) {
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f) || length(f) == 0) stop("None of these files exist: ", paste(candidates, collapse=", "))
  f
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
}

norm_str <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x
}

qsum_dt <- function(x) {
  qs <- quantile(x, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
  data.table(q05=qs[[1]], q25=qs[[2]], q50=qs[[3]], q75=qs[[4]], q95=qs[[5]])
}

# -----------------------------
# 1) Read input files
# -----------------------------
DATA_DIR <- "~/Desktop/hpai_challenge/Phase 1 documents"

csv_files <- list(
  pop   = "population.csv",
  cases = "cases.csv",
  act   = "activity.csv",
  prev  = "prev_culls.csv",
  mov   = "movement.csv"
)

geo_files <- list(
  hrz      = "hrz_32626.geojson",
  clc      = "clc_32626.geojson",
  counties = "counties_32626.geojson"
)

# Read CSV
pop   <- fread(file.path(DATA_DIR, csv_files$pop))
cases <- fread(file.path(DATA_DIR, csv_files$cases))
act   <- fread(file.path(DATA_DIR, csv_files$act))
prev  <- fread(file.path(DATA_DIR, csv_files$prev))
mov   <- fread(file.path(DATA_DIR, csv_files$mov))

# Read GeoJSON
hrz_file      <- st_read(file.path(DATA_DIR, geo_files$hrz), quiet = TRUE)
clc_file      <- st_read(file.path(DATA_DIR, geo_files$clc), quiet = TRUE)
counties_file <- st_read(file.path(DATA_DIR, geo_files$counties), quiet = TRUE)

# Output dir
OUT_Q2 <- file.path(DATA_DIR, "Q2_outputs")
dir.create(OUT_Q2, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) Parse/standardize columns
# -----------------------------
# population.csv
stop_if_missing(pop, c("farm_id","x","y","production"), "population.csv")
pop[, farm_id := as.character(farm_id)]
pop[, production := norm_str(production)]

# Optional capacity field (size scaling)
cap_col <- NULL
for (cc in c("capacity","cap","c_i","Capacity","CAPACITY")) {
  if (cc %in% names(pop)) { cap_col <- cc; break }
}

# cases.csv
stop_if_missing(cases, c("farm_id","date_confirmed"), "cases.csv")
cases[, farm_id := as.character(farm_id)]
date_cols_cases <- intersect(c("date_suspicious","date_confirmed","cull_start","cull_end"), names(cases))
for (cc in date_cols_cases) cases[[cc]] <- as_date_safe(cases[[cc]])

# activity.csv
stop_if_missing(act, c("farm_id","date_start","date_end"), "activity.csv")
act[, farm_id := as.character(farm_id)]
act[, date_start := as_date_safe(date_start)]
act[, date_end   := as_date_safe(date_end)]

# prev_culls.csv
stop_if_missing(prev, c("farm_id","cull_start"), "prev_culls.csv")
prev[, farm_id := as.character(farm_id)]
prev[, cull_start := as_date_safe(cull_start)]
if ("cull_end" %in% names(prev)) prev[, cull_end := as_date_safe(cull_end)]

# movements.csv  (BUGFIX: df_name should be string, not mov_file)
stop_if_missing(mov, c("date","source_farm","dest_farm","volume"), "movement.csv")
mov[, date := as_date_safe(date)]
mov[, source_farm := as.character(source_farm)]
mov[, dest_farm   := as.character(dest_farm)]
mov[, volume := as.numeric(volume)]
mov[is.na(volume), volume := 0]

# -----------------------------
# 3) Build farm index, coords, SF points
# -----------------------------
farm_ids <- pop$farm_id
N <- length(farm_ids)
id_to_idx <- setNames(seq_len(N), farm_ids)

coords <- as.matrix(pop[, .(x,y)])
pop_sf <- st_as_sf(pop, coords = c("x","y"), crs = 32626, remove = FALSE)

# Size scaling w_size (optional)
if (!is.null(cap_col)) {
  cap <- suppressWarnings(as.numeric(pop[[cap_col]]))
  ok <- is.finite(cap) & cap > 0
  if (!any(ok)) {
    w_size <- rep(1.0, N)
  } else {
    cap[!ok] <- median(cap[ok], na.rm = TRUE)
    w_size <- cap / mean(cap, na.rm = TRUE)
  }
} else {
  w_size <- rep(1.0, N)
}

# -----------------------------
# 4) HRZ membership H_i
# -----------------------------
hrz_sf <- hrz_file
if (!is.na(st_crs(hrz_sf)$epsg) && st_crs(hrz_sf)$epsg != 32626) hrz_sf <- st_transform(hrz_sf, 32626)
pop$H_i <- as.integer(lengths(st_intersects(pop_sf, hrz_sf)) > 0)

# -----------------------------
# 5) Environmental covariate Z_i from CLC wetlands/water proximity
# -----------------------------
clc <- clc_file
if (!is.na(st_crs(clc)$epsg) && st_crs(clc)$epsg != 32626) clc <- st_transform(clc, 32626)

stop_if_missing(st_drop_geometry(clc), c("CODE_12"), "clc_32626.geojson")
clc$CODE_12 <- suppressWarnings(as.integer(as.character(clc$CODE_12)))

wet_codes <- c(411,412,421,422,423,511,512,521,522,523)
clc_wet <- clc[clc$CODE_12 %in% wet_codes, ]
if (nrow(clc_wet) == 0) stop("CLC wetland/water filter returned 0 features; check CODE_12 values/types.")

nearest_idx <- st_nearest_feature(pop_sf, clc_wet)
nearest_geom <- clc_wet[nearest_idx, ]
pop$d_wat_m <- as.numeric(st_distance(pop_sf, nearest_geom, by_element = TRUE))

kappa_m <- 10000
pop$Z_i <- exp(-pop$d_wat_m / kappa_m)

# -----------------------------
# 6) County membership (reporting only)
# -----------------------------
counties <- counties_file
if (!is.na(st_crs(counties)$epsg) && st_crs(counties)$epsg != 32626) counties <- st_transform(counties, 32626)

pop_join <- st_join(pop_sf, counties, join = st_within, left = TRUE)
county_field_candidates <- c("NAME","name","COUNTY","county","County","COUNTY_NAM","COUNTYNAME")
county_field <- county_field_candidates[county_field_candidates %in% names(pop_join)][1]
if (is.na(county_field)) county_field <- setdiff(names(counties), attr(counties, "sf_column"))[1]
pop$county_name <- as.character(pop_join[[county_field]])

# -----------------------------
# 7) Activity indicator A_i(t)
# -----------------------------
act_list <- split(act[, .(farm_id, date_start, date_end)], by = "farm_id", keep.by = FALSE)

is_active_on <- function(fid, t) {
  w <- act_list[[as.character(fid)]]
  if (is.null(w)) return(FALSE)
  any(t >= w$date_start & t <= w$date_end)
}

active_vec_on <- function(t, farm_ids) {
  as.integer(vapply(farm_ids, function(fid) is_active_on(fid, t), logical(1)))
}

# -----------------------------
# 8) Preprocess movements + preventive culls for fast daily access
# -----------------------------
mov <- mov[source_farm %in% farm_ids & dest_farm %in% farm_ids]
mov[, src_idx := id_to_idx[source_farm]]
mov[, dst_idx := id_to_idx[dest_farm]]
mov_by_date <- split(mov, by = "date", keep.by = FALSE)

prev <- prev[farm_id %in% farm_ids]
prev[, idx := id_to_idx[farm_id]]
prev_by_date <- split(prev, by = "cull_start", keep.by = FALSE)

# -----------------------------
# 9) Spatial neighbor lists (radius search) with aligned idx/dist
# -----------------------------
spatial_cutoff_m <- 50000
nn <- nn2(data = coords, query = coords, searchtype = "radius", radius = spatial_cutoff_m)

neighbor_idx  <- vector("list", N)
neighbor_dist <- vector("list", N)
for (i in seq_len(N)) {
  idx <- nn$nn.idx[i, ]
  dst <- nn$nn.dists[i, ]
  keep <- idx > 0 & idx != i
  neighbor_idx[[i]]  <- idx[keep]
  neighbor_dist[[i]] <- dst[keep]
}

K_exp <- function(d, alpha) exp(-d / alpha)
g_vol <- function(v) log1p(pmax(v, 0))

# -----------------------------
# 10) SEIDR Simulator (YOUR model + minimal additions for Q2)
#     Add init_state + farm-level flags for spatial risk
# -----------------------------
# ============================================================
# Precompute neighbor edge list (run ONCE, outside simulate)
# ============================================================
N <- length(farm_ids)

edge_i <- rep.int(seq_len(N), lengths(neighbor_idx))          # target i repeated
edge_j <- unlist(neighbor_idx,  use.names = FALSE)            # source j
edge_d <- unlist(neighbor_dist, use.names = FALSE)            # distance d_ij
edge_w <- w_size[edge_j]                                      # w_size of source j

stopifnot(length(edge_i) == length(edge_j),
          length(edge_j) == length(edge_d),
          length(edge_d) == length(edge_w))




simulate_seidr <- function(start_date, end_date,
                           pop, farm_ids, id_to_idx,
                           mov_by_date, prev_by_date,
                           # neighbor_idx, neighbor_dist,  # not needed anymore inside
                           w_size,
                           params,
                           # NEW: edge list inputs
                           edge_i, edge_j, edge_d, edge_w,
                           seed_farms = character(0),
                           init_state = NULL,
                           force_observed = FALSE,
                           cases = NULL,
                           rng_seed = NULL) {
  
  if (!is.null(rng_seed)) set.seed(rng_seed)
  
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day")
  TT <- length(dates)
  N <- length(farm_ids)
  
  # State codes: 1=S, 2=E, 3=I, 4=D, 5=R
  S <- 1L; E <- 2L; I <- 3L; D <- 4L; R <- 5L
  
  if (!is.null(init_state)) {
    X <- as.integer(init_state)
    if (length(X) != N) stop("init_state must have length N")
  } else {
    X <- rep.int(S, N)
    if (length(seed_farms) > 0) {
      seed_idx <- id_to_idx[as.character(seed_farms)]
      seed_idx <- seed_idx[!is.na(seed_idx)]
      X[seed_idx] <- E
    }
  }
  
  baseline_nonS  <- (X != S)
  newly_infected <- rep(FALSE, N)
  newly_detected <- rep(FALSE, N)
  
  forced_confirm <- rep(as.Date(NA), N)
  forced_cull    <- rep(as.Date(NA), N)
  if (force_observed && !is.null(cases)) {
    cc <- cases[farm_id %in% farm_ids]
    cc[, idx := id_to_idx[farm_id]]
    cc <- cc[order(date_confirmed)]
    cc <- cc[!duplicated(idx)]
    forced_confirm[cc$idx] <- cc$date_confirmed
    if ("cull_start" %in% names(cc)) forced_cull[cc$idx] <- cc$cull_start
  }
  
  p_EI <- 1 / params$L_latent
  p_ID <- 1 / params$L_inf_to_det
  p_DR <- 1 / params$L_det_to_rem
  
  L_mv <- params$L_mv
  if (!is.finite(L_mv) || L_mv < 1) stop("params$L_mv must be >= 1")
  empty_mv <- data.table(src_idx=integer(), dst_idx=integer(), volume=numeric())
  mv_queue <- replicate(L_mv, empty_mv, simplify = FALSE)
  
  gate_progression_by_activity <- isTRUE(params$gate_progression_by_activity)
  gate_detection_by_activity   <- isTRUE(params$gate_detection_by_activity)
  gate_removal_by_activity     <- isTRUE(params$gate_removal_by_activity)
  
  # small speed: compute once
  is_b2 <- (pop$production == "broiler_2")
  
  out <- data.table(
    date = dates,
    newE = integer(TT), newI = integer(TT), newD = integer(TT), newR = integer(TT),
    S = integer(TT), E = integer(TT), I = integer(TT), D = integer(TT), R = integer(TT)
  )
  
  # reuse vectors (less allocation)
  inf_w  <- numeric(N)
  lambda <- numeric(N)
  
  for (tt in seq_along(dates)) {
    t <- dates[tt]
    
    A <- active_vec_on(t, farm_ids)
    
    pv <- prev_by_date[[as.character(t)]]
    if (!is.null(pv) && nrow(pv) > 0) X[pv$idx] <- R
    
    mv_today <- mov_by_date[[as.character(t)]]
    if (is.null(mv_today)) mv_today <- empty_mv
    mv_today <- mv_today[, .(src_idx, dst_idx, volume)]
    if (L_mv > 1) mv_queue <- c(list(mv_today), mv_queue[1:(L_mv-1)]) else mv_queue <- list(mv_today)
    mv_recent <- rbindlist(mv_queue, use.names = TRUE, fill = TRUE)
    
    # inf_w (reuse)
    inf_w[] <- 0.0
    inf_w[X == I] <- 1.0
    inf_w[X == D] <- params$eta
    inf_w[A == 0L] <- 0.0
    
    # baseline hazard (reuse)
    lambda[] <- params$beta_hrz * pop$H_i + params$beta_env * pop$Z_i
    
    # -----------------------------
    # Neighbor transmission (FAST)
    # -----------------------------
    iw_edge <- inf_w[edge_j]
    if (any(iw_edge > 0)) {
      k_edge <- K_exp(edge_d, params$alpha)
      edge_contrib <- iw_edge * edge_w * k_edge
      
      rs <- rowsum(edge_contrib, edge_i, reorder = FALSE)
      nb_sum <- numeric(N)
      nb_sum[as.integer(rownames(rs))] <- rs[, 1]
      
      lambda <- lambda + params$beta_sp * nb_sum
    }
    
    # -----------------------------
    # Movement transmission (unchanged logic, but is_b2 precomputed)
    # -----------------------------
    sus_idx <- which(X == S & A == 1L)
    
    if (nrow(mv_recent) > 0 && length(sus_idx) > 0) {
      mv2 <- mv_recent[dst_idx %in% sus_idx]
      if (nrow(mv2) > 0) {
        mv2[, src_inf := inf_w[src_idx]]
        mv2 <- mv2[src_inf > 0]
        if (nrow(mv2) > 0) {
          mv2 <- mv2[is_b2[dst_idx]]
          if (nrow(mv2) > 0) {
            mv2[, contrib := g_vol(volume) * src_inf]
            mv_sum <- mv2[, .(mv_contrib = sum(contrib)), by = dst_idx]
            lambda[mv_sum$dst_idx] <- lambda[mv_sum$dst_idx] + params$beta_mv * mv_sum$mv_contrib
          }
        }
      }
    }
    
    # S->E
    p_inf <- 1 - exp(-pmax(lambda, 0))
    newE_vec <- rbinom(N, 1, p_inf)
    newE_vec[!(X == S & A == 1L)] <- 0L
    if (any(newE_vec == 1L)) {
      idxE <- which(newE_vec == 1L)
      newly_infected[idxE] <- newly_infected[idxE] | (!baseline_nonS[idxE])
      X[idxE] <- E
    }
    
    # E->I
    newI_vec <- rbinom(N, 1, p_EI)
    if (gate_progression_by_activity) newI_vec[!(X == E & A == 1L)] <- 0L else newI_vec[X != E] <- 0L
    X[newI_vec == 1L] <- I
    
    # I->D
    newD_vec <- integer(N)
    if (force_observed) {
      forceD <- which(!is.na(forced_confirm) & forced_confirm == t)
      forceD <- forceD[X[forceD] != R]
      if (length(forceD) > 0) {
        before <- X[forceD]
        X[forceD] <- D
        newD_vec[forceD] <- as.integer(before != D)
      }
      cand <- setdiff(which(X == I), forceD)
    } else {
      cand <- which(X == I)
    }
    
    if (gate_detection_by_activity) cand <- cand[A[cand] == 1L]
    if (length(cand) > 0) {
      det <- rbinom(length(cand), 1, p_ID)
      toD <- cand[det == 1L]
      if (length(toD) > 0) {
        X[toD] <- D
        newD_vec[toD] <- 1L
        newly_detected[toD] <- TRUE
      }
    }
    
    # D->R
    newR_vec <- integer(N)
    if (force_observed) {
      forceR <- which(!is.na(forced_cull) & forced_cull == t)
      forceR <- forceR[X[forceR] != R]
      if (length(forceR) > 0) {
        X[forceR] <- R
        newR_vec[forceR] <- 1L
      }
      candR <- setdiff(which(X == D), forceR)
    } else {
      candR <- which(X == D)
    }
    
    if (gate_removal_by_activity) candR <- candR[A[candR] == 1L]
    if (length(candR) > 0) {
      rem <- rbinom(length(candR), 1, p_DR)
      toR <- candR[rem == 1L]
      X[toR] <- R
      newR_vec[toR] <- 1L
    }
    
    out[tt, `:=`(
      newE = sum(newE_vec),
      newI = sum(newI_vec),
      newD = sum(newD_vec),
      newR = sum(newR_vec),
      S = sum(X == 1L),
      E = sum(X == 2L),
      I = sum(X == 3L),
      D = sum(X == 4L),
      R = sum(X == 5L)
    )]
  }
  
  list(
    out_daily = out,
    final_state = X,
    newly_infected = newly_infected,
    newly_detected = newly_detected
  )
}

# -----------------------------
# 11) Baseline parameters
# -----------------------------
params <- list(
  beta_sp  = 0.08,
  alpha    = 12000,
  beta_mv  = 0.20,
  beta_hrz = 0.002,
  beta_env = 0.001,
  
  eta      = 0.4,
  L_latent = 2,
  L_inf_to_det = 4,
  L_det_to_rem = 3,
  L_mv     = 7,
  
  gate_progression_by_activity = FALSE,
  gate_detection_by_activity   = FALSE,
  gate_removal_by_activity     = FALSE
)

# -----------------------------
# 12) Seed farm (your rule)
# -----------------------------
seed_candidates <- pop[H_i == 1 & production == "broiler_1", farm_id]
seed <- if (length(seed_candidates) > 0) seed_candidates[1] else pop$farm_id[1]

# -----------------------------
# 13) Define "today" + forecast horizon
# -----------------------------
conf_dates <- cases$date_confirmed
conf_dates <- conf_dates[!is.na(conf_dates)]
if (length(conf_dates) == 0) stop("cases.csv has no non-NA date_confirmed; cannot define simulation window.")

start_sim <- min(conf_dates)
today     <- max(conf_dates)   # "as it is today" = last observed confirmation date in dataset
end_fc    <- today + 28

# -----------------------------
# 14) OPTIONAL: Grid-search calibration (SSE vs observed confirmations)
#     Turn on/off here
# -----------------------------
do_grid_search <- F

if (do_grid_search) {
  
  # ----------------------------------------
  # 1) Build observed series + last-21-day window
  # ----------------------------------------
  obs_full <- cases[!is.na(date_confirmed), .N, by = date_confirmed][order(date_confirmed)]
  setnames(obs_full, "date_confirmed", "date")
  obs_full[, date := as.Date(date)]
  
  end_date_cal   <- max(obs_full$date)
  start_date_cal <- end_date_cal - 21L
  
  # keep only last 21 days (inclusive)
  obs <- obs_full[date >= start_date_cal & date <= end_date_cal]
  
  cat("Calibration window:", as.character(start_date_cal), "to", as.character(end_date_cal),
      " (", nrow(obs), " observed days)\n", sep = "")
  
  # fixed daily grid for fast scoring
  all_dates <- seq(start_date_cal, end_date_cal, by = "day")
  obs_vec <- integer(length(all_dates))
  obs_vec[match(obs$date, all_dates)] <- obs$N
  
  # ----------------------------------------
  # 2) Fast scoring (NO merge)
  # ----------------------------------------
  score_run_fast <- function(sim_out) {
    df <- sim_out$out_daily[, .(date = as.Date(date), pred = newD)]
    
    pred_vec <- numeric(length(all_dates))
    ii <- match(df$date, all_dates)
    ok <- !is.na(ii)
    pred_vec[ii[ok]] <- df$pred[ok]
    
    sum((obs_vec - pred_vec)^2)
  }
  
  # ----------------------------------------
  # 3) Grid definition
  # ----------------------------------------
  grid <- CJ(
    beta_sp  = c(0.04, 0.08, 0.12),
    alpha    = c(8000, 12000, 20000),
    beta_mv  = c(0.10, 0.20, 0.30),
    beta_hrz = c(0.001, 0.002, 0.004),
    beta_env = c(0.0005, 0.0010, 0.0020)
  )
  grid[, score := NA_real_]
  
  n_sims_per_point <- 10
  
  # ----------------------------------------
  # 4) Grid search loop (short window => faster)
  # ----------------------------------------
  for (k in seq_len(nrow(grid))) {
    
    params_k <- params
    params_k$beta_sp  <- grid$beta_sp[k]
    params_k$alpha    <- grid$alpha[k]
    params_k$beta_mv  <- grid$beta_mv[k]
    params_k$beta_hrz <- grid$beta_hrz[k]
    params_k$beta_env <- grid$beta_env[k]
    
    scores <- numeric(n_sims_per_point)
    
    for (r in seq_len(n_sims_per_point)) {
      
      sim_k <- simulate_seidr(
        start_date = start_date_cal,
        end_date   = end_date_cal,
        pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
        mov_by_date = mov_by_date, prev_by_date = prev_by_date,
        w_size = w_size,
        params = params_k,
        edge_i = edge_i, edge_j = edge_j, edge_d = edge_d, edge_w = edge_w,
        seed_farms = seed,
        init_state = NULL,
        force_observed = FALSE,
        rng_seed = 1000 + 100*k + r
      )
      
      scores[r] <- score_run_fast(sim_k)
    }
    
    grid$score[k] <- mean(scores)
    
    if (k %% 10 == 0) {
      cat(sprintf("Grid %d / %d (%.1f%%)\n", k, nrow(grid), 100*k/nrow(grid)))
    }
  }
  
  # ----------------------------------------
  # 5) Pick best + set params
  # ----------------------------------------
  best <- grid[order(score)][1:10]
  print(best)
  
  if (nrow(best) >= 1) {
    params$beta_sp  <- best$beta_sp[1]
    params$alpha    <- best$alpha[1]
    params$beta_mv  <- best$beta_mv[1]
    params$beta_hrz <- best$beta_hrz[1]
    params$beta_env <- best$beta_env[1]
    cat("Using calibrated params from best grid row (last-21-day fit).\n")
  }
}


# -----------------------------
# 15) Q2 Ensemble Forecast (two-stage per replicate)
# -----------------------------
run_one_q2 <- function(sim_id, rng_seed) {
  
  # A) Assimilation: start_sim -> today, force observed
  assim <- simulate_seidr(
    start_date = start_sim,
    end_date   = today,
    pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
    mov_by_date = mov_by_date, prev_by_date = prev_by_date,
    w_size = w_size,
    params = params,
    edge_i = edge_i, edge_j = edge_j, edge_d = edge_d, edge_w = edge_w,  # <-- add
    seed_farms = seed,
    init_state = NULL,
    force_observed = TRUE,
    cases = cases,
    rng_seed = rng_seed
  )
  
  init_state <- assim$final_state
  
  # B) Forecast: today+1 -> today+28, no forcing
  fc <- simulate_seidr(
    start_date = today + 1,
    end_date   = end_fc,
    pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
    mov_by_date = mov_by_date, prev_by_date = prev_by_date,
    w_size = w_size,
    params = params,
    edge_i = edge_i, edge_j = edge_j, edge_d = edge_d, edge_w = edge_w,  # <-- add
    seed_farms = character(0),
    init_state = init_state,
    force_observed = FALSE,
    cases = NULL,
    rng_seed = rng_seed + 9999
  )
  
  fc$out_daily[, sim := sim_id]
  
  farmflags <- data.table(
    sim = sim_id,
    farm_id = farm_ids,
    county_name = pop$county_name,
    production = pop$production,
    newly_infected = fc$newly_infected,
    newly_detected = fc$newly_detected
  )
  
  list(traj = fc$out_daily, farmflags = farmflags)
}
# Choose number of stochastic trajectories
n_sims <- 10  # raise to 500–1000 if you have time
set.seed(123)

res <- vector("list", n_sims)
for (s in seq_len(n_sims)) {
  res[[s]] <- run_one_q2(sim_id = s, rng_seed = 10000 + s)
  if (s %% 25 == 0) cat("Completed", s, "of", n_sims, "\n")
}

traj_all <- rbindlist(lapply(res, `[[`, "traj"))
farmflags_all <- rbindlist(lapply(res, `[[`, "farmflags"))

# -----------------------------
# 16) Save RAW trajectories (requested)
# -----------------------------
fwrite(traj_all, file.path(OUT_Q2, "Q2_raw_simulated_trajectories_daily.csv"))
fwrite(farmflags_all, file.path(OUT_Q2, "Q2_raw_farmflags_new_events.csv"))

# -----------------------------
# 17) Temporal summary: daily quantiles (overall)
# -----------------------------
sum_daily <- traj_all[, c(qsum_dt(newE), qsum_dt(newD), qsum_dt(newR)), by = date]

# sum_daily should be: date + 15 quantile columns
stopifnot(ncol(sum_daily) == 16)

new_names <- c(
  "date",
  "newE_q05","newE_q25","newE_q50","newE_q75","newE_q95",
  "newD_q05","newD_q25","newD_q50","newD_q75","newD_q95",
  "newR_q05","newR_q25","newR_q50","newR_q75","newR_q95"
)

data.table::setnames(sum_daily, new = new_names)

stopifnot(anyDuplicated(names(sum_daily)) == 0)

fwrite(sum_daily, file.path(OUT_Q2, "Q2_summary_daily_quantiles_overall.csv"))

# -----------------------------
# 18) By production: 4-week totals with intervals
#      (count farms that become newly infected/detected over 28 days)
# -----------------------------
prod_sim <- farmflags_all[
  , .(
    n_new_infected = sum(newly_infected, na.rm = TRUE),
    n_new_detected = sum(newly_detected, na.rm = TRUE)
  ),
  by = .(sim, production)
]

prod_summary <- prod_sim[
  , .(
    new_infected_q05 = as.numeric(quantile(n_new_infected, 0.05, na.rm = TRUE)),
    new_infected_q50 = as.numeric(quantile(n_new_infected, 0.50, na.rm = TRUE)),
    new_infected_q95 = as.numeric(quantile(n_new_infected, 0.95, na.rm = TRUE)),
    new_detected_q05 = as.numeric(quantile(n_new_detected, 0.05, na.rm = TRUE)),
    new_detected_q50 = as.numeric(quantile(n_new_detected, 0.50, na.rm = TRUE)),
    new_detected_q95 = as.numeric(quantile(n_new_detected, 0.95, na.rm = TRUE))
  ),
  by = production
][order(-new_detected_q50)]

fwrite(prod_summary, file.path(OUT_Q2, "Q2_summary_4week_totals_by_production.csv"))

# -----------------------------
# 19) Spatial risk (county): probability of ≥1 new event over 4 weeks
# -----------------------------
# Clean county_name to avoid NA groups
farmflags_all[is.na(county_name) | county_name == "", county_name := "Unknown"]

county_sim <- farmflags_all[
  , .(
    any_new_infected = any(newly_infected, na.rm = TRUE),
    any_new_detected = any(newly_detected, na.rm = TRUE)
  ),
  by = .(sim, county_name)
]

county_risk <- county_sim[
  , .(
    p_any_new_infected = mean(any_new_infected, na.rm = TRUE),
    p_any_new_detected = mean(any_new_detected, na.rm = TRUE)
  ),
  by = county_name
][order(-p_any_new_detected)]

fwrite(county_risk, file.path(OUT_Q2, "Q2_spatial_risk_by_county_end_of_4weeks.csv"))

# County x production risk
county_prod_sim <- farmflags_all[
  , .(
    any_new_infected = any(newly_infected, na.rm = TRUE),
    any_new_detected = any(newly_detected, na.rm = TRUE)
  ),
  by = .(sim, county_name, production)
]

county_prod_risk <- county_prod_sim[
  , .(
    p_any_new_infected = mean(any_new_infected, na.rm = TRUE),
    p_any_new_detected = mean(any_new_detected, na.rm = TRUE)
  ),
  by = .(county_name, production)
][order(-p_any_new_detected)]

fwrite(county_prod_risk, file.path(OUT_Q2, "Q2_spatial_risk_by_county_and_production_end_of_4weeks.csv"))

# -----------------------------
# 20) Quick forecast figure: newD (median + 90% interval)
# -----------------------------
p_fc_newD <- ggplot(as.data.frame(sum_daily), aes(x = date)) +
  geom_ribbon(aes(ymin = newD_q05, ymax = newD_q95), alpha = 0.25) +
  geom_line(aes(y = newD_q50), linewidth = 1) +
  labs(
    title = "Forecast (4 weeks): daily detections (newD) — median and 90% interval",
    x = "Date", y = "Daily detections"
  ) +
  theme_minimal()

ggsave(
  file.path(OUT_Q2, "Q2_forecast_newD_overall.png"),
  p_fc_newD, width = 8, height = 4, dpi = 300, bg = "white"
)

# -----------------------------
# 21) Print what was created
# -----------------------------
cat("\n================ Q2 DONE ================\n")
cat("Outputs in: ", OUT_Q2, "\n\n")
cat("Raw trajectories: Q2_raw_simulated_trajectories_daily.csv\n")
cat("Farm flags:       Q2_raw_farmflags_new_events.csv\n")
cat("Daily quantiles:  Q2_summary_daily_quantiles_overall.csv\n")
cat("Prod totals:      Q2_summary_4week_totals_by_production.csv\n")
cat("County risk:      Q2_spatial_risk_by_county_end_of_4weeks.csv\n")
cat("County+prod risk: Q2_spatial_risk_by_county_and_production_end_of_4weeks.csv\n")
cat("Figure:           Q2_forecast_newD_overall.png\n")
cat("========================================\n")

