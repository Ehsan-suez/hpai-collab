# ============================================================
# Farm-level stochastic SEIDR model (daily)
# Spatial (exponential kernel) + broiler movements + HRZ + environment (CLC wetlands/water)
# Outputs: daily counts of states + new infections/detections/removals
# Includes: optional grid-search calibration (SSE vs observed confirmations)
#

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(sf)
  library(RANN)
  library(data.table)
  library(ggplot2)
  library(sf)
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

# -----------------------------
# 1) Read input files
# -----------------------------
pop_file   <- "population.csv"
cases_file <- "cases.csv"
act_file   <- "activity.csv"
prev_file  <- "prev_culls.csv"

mov_file <- pick_existing_file(c("movement.csv", "movements.csv"))

hrz_file      <- "hrz_32626.geojson"
clc_file      <- "clc_32626.geojson"
counties_file <- "counties_32626.geojson"

pop   <- fread(pop_file)
cases <- fread(cases_file)
act   <- fread(act_file)
prev  <- fread(prev_file)
mov   <- fread(mov_file)

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

# movements.csv
stop_if_missing(mov, c("date","source_farm","dest_farm","volume"), mov_file)
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
hrz_sf <- st_read(hrz_file, quiet = TRUE)
if (st_crs(hrz_sf)$epsg != 32626) hrz_sf <- st_transform(hrz_sf, 32626)
pop$H_i <- as.integer(lengths(st_intersects(pop_sf, hrz_sf)) > 0)

# -----------------------------
# 5) Environmental covariate Z_i from CLC wetlands/water proximity
# -----------------------------
clc <- st_read(clc_file, quiet = TRUE)
if (st_crs(clc)$epsg != 32626) clc <- st_transform(clc, 32626)

stop_if_missing(st_drop_geometry(clc), c("CODE_12"), "clc_32626.geojson")
clc$CODE_12 <- suppressWarnings(as.integer(as.character(clc$CODE_12)))

# CLC Level-3 wetlands + water (include 511)
wet_codes <- c(411,412,421,422,423,511,512,521,522,523)
clc_wet <- clc[clc$CODE_12 %in% wet_codes, ]
if (nrow(clc_wet) == 0) stop("CLC wetland/water filter returned 0 features; check CODE_12 values/types.")

nearest_idx <- st_nearest_feature(pop_sf, clc_wet)
nearest_geom <- clc_wet[nearest_idx, ]
pop$d_wat_m <- as.numeric(st_distance(pop_sf, nearest_geom, by_element = TRUE))

kappa_m <- 10000  # 10 km baseline (sensitivity parameter)
pop$Z_i <- exp(-pop$d_wat_m / kappa_m)

# -----------------------------
# 6) County membership (reporting only)
# -----------------------------
counties <- st_read(counties_file, quiet = TRUE)
if (st_crs(counties)$epsg != 32626) counties <- st_transform(counties, 32626)

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
spatial_cutoff_m <- 50000  # 50 km
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
# 10) SEIDR Simulator
# -----------------------------
simulate_seidr <- function(start_date, end_date,
                           pop, farm_ids, id_to_idx,
                           mov_by_date, prev_by_date,
                           neighbor_idx, neighbor_dist,
                           w_size,
                           params,
                           seed_farms = character(0),
                           force_observed = FALSE,
                           cases = NULL,
                           rng_seed = NULL) {
  
  if (!is.null(rng_seed)) set.seed(rng_seed)
  
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day")
  TT <- length(dates)
  N <- length(farm_ids)
  
  # State codes: 1=S, 2=E, 3=I, 4=D, 5=R
  S <- 1L; E <- 2L; I <- 3L; D <- 4L; R <- 5L
  X <- rep(S, N)
  
  # Seed infections
  if (length(seed_farms) > 0) {
    seed_idx <- id_to_idx[as.character(seed_farms)]
    seed_idx <- seed_idx[!is.na(seed_idx)]
    X[seed_idx] <- E
  }
  
  # Forced confirmation/cull dates (optional conditioning)
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
  
  # Transition probabilities (geometric via mean durations)
  p_EI <- 1 / params$L_latent
  p_ID <- 1 / params$L_inf_to_det
  p_DR <- 1 / params$L_det_to_rem
  
  # Movement window queue
  L_mv <- params$L_mv
  if (!is.finite(L_mv) || L_mv < 1) stop("params$L_mv must be >= 1")
  empty_mv <- data.table(src_idx=integer(), dst_idx=integer(), volume=numeric())
  mv_queue <- replicate(L_mv, empty_mv, simplify = FALSE)
  
  # optional gating toggles (default: transmission only)
  gate_progression_by_activity <- isTRUE(params$gate_progression_by_activity)
  gate_detection_by_activity   <- isTRUE(params$gate_detection_by_activity)
  gate_removal_by_activity     <- isTRUE(params$gate_removal_by_activity)
  
  out <- data.table(
    date = dates,
    newE = integer(TT), newI = integer(TT), newD = integer(TT), newR = integer(TT),
    S = integer(TT), E = integer(TT), I = integer(TT), D = integer(TT), R = integer(TT)
  )
  
  for (tt in seq_along(dates)) {
    t <- dates[tt]
    
    # activity today
    A <- active_vec_on(t, farm_ids)
    
    # preventive removals (recorded)
    pv <- prev_by_date[[as.character(t)]]
    if (!is.null(pv) && nrow(pv) > 0) X[pv$idx] <- R
    
    # update movement queue
    mv_today <- mov_by_date[[as.character(t)]]
    if (is.null(mv_today)) mv_today <- empty_mv
    mv_today <- mv_today[, .(src_idx, dst_idx, volume)]
    if (L_mv > 1) mv_queue <- c(list(mv_today), mv_queue[1:(L_mv-1)]) else mv_queue <- list(mv_today)
    mv_recent <- rbindlist(mv_queue, use.names = TRUE, fill = TRUE)
    
    # infectiousness modifier (transmission only)
    inf_w <- numeric(N)
    inf_w[X == I] <- 1.0
    inf_w[X == D] <- params$eta
    inf_w[A == 0L] <- 0.0  # inactive farms do not transmit
    
    # ---- FOI ----
    lambda <- numeric(N)
    
    # background (only matters for susceptible & active because S->E is gated)
    lambda <- lambda + params$beta_hrz * pop$H_i + params$beta_env * pop$Z_i
    
    # spatial: only compute for susceptible & active
    sus_idx <- which(X == S & A == 1L)
    for (i in sus_idx) {
      nbs <- neighbor_idx[[i]]
      if (length(nbs) == 0) next
      iw <- inf_w[nbs]
      if (all(iw == 0)) next
      dists <- neighbor_dist[[i]]
      lambda[i] <- lambda[i] + params$beta_sp * sum(iw * w_size[nbs] * K_exp(dists, params$alpha))
    }
    
    # movement: only broiler_2 recipients that are susceptible & active
    if (nrow(mv_recent) > 0 && length(sus_idx) > 0) {
      mv2 <- mv_recent[dst_idx %in% sus_idx]
      if (nrow(mv2) > 0) {
        mv2[, src_inf := inf_w[src_idx]]
        mv2 <- mv2[src_inf > 0]
        if (nrow(mv2) > 0) {
          is_b2 <- (pop$production == "broiler_2")
          mv2 <- mv2[is_b2[dst_idx]]
          if (nrow(mv2) > 0) {
            mv2[, contrib := g_vol(volume) * src_inf]
            mv_sum <- mv2[, .(mv_contrib = sum(contrib)), by = dst_idx]
            lambda[mv_sum$dst_idx] <- lambda[mv_sum$dst_idx] + params$beta_mv * mv_sum$mv_contrib
          }
        }
      }
    }
    
    # S->E (TRANSMISSION gated by activity)
    p_inf <- 1 - exp(-pmax(lambda, 0))
    newE_vec <- rbinom(N, 1, p_inf)
    newE_vec[!(X == S & A == 1L)] <- 0L
    X[newE_vec == 1L] <- E
    
    # E->I (default: NOT gated by activity)
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
      X[toD] <- D
      newD_vec[toD] <- 1L
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
  
  list(out_daily = out, final_state = X)
}

# -----------------------------
# 11) Baseline parameters
# -----------------------------
params <- list(
  beta_sp  = 0.08,     # estimate
  alpha    = 12000,    # estimate (meters)
  beta_mv  = 0.20,     # estimate
  beta_hrz = 0.002,    # estimate
  beta_env = 0.001,    # estimate
  
  eta      = 0.4,      # fixed/sensitivity
  L_latent = 2,        # fixed/sensitivity
  L_inf_to_det = 4,    # fixed/sensitivity
  L_det_to_rem = 3,    # fixed/sensitivity
  L_mv     = 7,        # fixed (>=1)
  
  # gating toggles (DEFAULT FALSE = transmission-only gating)
  gate_progression_by_activity = FALSE,
  gate_detection_by_activity   = FALSE,
  gate_removal_by_activity     = FALSE
)

# -----------------------------
# 12) Pick a seed farm (HRZ broiler_1 if possible)
# -----------------------------
seed_candidates <- pop[H_i == 1 & production == "broiler_1", farm_id]
seed <- if (length(seed_candidates) > 0) seed_candidates[1] else pop$farm_id[1]

# -----------------------------
# 13) Run one simulation example
# -----------------------------
conf_dates <- cases$date_confirmed
conf_dates <- conf_dates[!is.na(conf_dates)]
if (length(conf_dates) == 0) stop("cases.csv has no non-NA date_confirmed; cannot define simulation window.")

start_sim <- min(conf_dates)
last_obs <- max(conf_dates)
end_sim <- last_obs + 28

sim <- simulate_seidr(
  start_date = start_sim,
  end_date   = end_sim,
  pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
  mov_by_date = mov_by_date, prev_by_date = prev_by_date,
  neighbor_idx = neighbor_idx, neighbor_dist = neighbor_dist,
  w_size = w_size,
  params = params,
  seed_farms = seed,
  force_observed = TRUE,   # set TRUE to condition on confirmations/culls
  cases = cases,
  rng_seed = 123
)

print(sim$out_daily)

# -----------------------------
# 14) OPTIONAL: Grid-search calibration (SSE vs observed confirmations)
#     Average over multiple sims per grid point to reduce MC noise
# -----------------------------
do_grid_search <- FALSE

if (do_grid_search) {
  obs <- cases[!is.na(date_confirmed), .N, by = date_confirmed][order(date_confirmed)]
  setnames(obs, "date_confirmed", "date")
  
  score_run <- function(sim_out, obs) {
    df <- sim_out$out_daily[, .(date, pred = newD)]
    m <- merge(obs, df, by = "date", all.x = TRUE)
    m[is.na(pred), pred := 0]
    sum((m$N - m$pred)^2)
  }
  
  grid <- CJ(
    beta_sp  = c(0.04, 0.08, 0.12),
    alpha    = c(8000, 12000, 20000),
    beta_mv  = c(0.10, 0.20, 0.30),
    beta_hrz = c(0.001, 0.002, 0.004),
    beta_env = c(0.0005, 0.0010, 0.0020)
  )
  grid[, score := NA_real_]
  
  n_sims_per_point <- 10
  
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
        start_date = min(obs$date),
        end_date   = max(obs$date),
        pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
        mov_by_date = mov_by_date, prev_by_date = prev_by_date,
        neighbor_idx = neighbor_idx, neighbor_dist = neighbor_dist,
        w_size = w_size,
        params = params_k,
        seed_farms = seed,
        force_observed = FALSE,
        rng_seed = 1000 + 100*k + r
      )
      scores[r] <- score_run(sim_k, obs)
    }
    
    grid$score[k] <- mean(scores)
  }
  
  best <- grid[order(score)][1:10]
  print(best)
}

# -----------------------------
# 15) OPTIONAL: Save outputs
# -----------------------------
# fwrite(sim$out_daily, "sim_out_daily.csv")
# fwrite(pop[, .(farm_id, county_name, H_i, d_wat_m, Z_i)], "farm_covariates.csv")


# -----------------------------
# 0) Safety checks
# -----------------------------
stopifnot(is.list(sim), "out_daily" %in% names(sim), "final_state" %in% names(sim))
stopifnot(is.data.table(sim$out_daily))
stopifnot(is.data.table(cases))
stopifnot(inherits(pop_sf, "sf"))
stopifnot(inherits(counties, "sf"))

# Ensure date is Date
sim$out_daily[, date := as.Date(date)]
cases[, date_confirmed := as.Date(date_confirmed)]

# -----------------------------
# 1) States over time (E, I, D)
# -----------------------------
df <- copy(sim$out_daily)

df_long <- melt(df,
                id.vars = "date",
                measure.vars = c("E","I","D"),
                variable.name = "state",
                value.name = "count")

p_states <- ggplot(df_long, aes(x = date, y = count, group = state, color = state)) +
  geom_line(linewidth = 1) +
  labs(x = "Date", y = "Number of farms", title = "Farm states over time (E/I/D)") +
  theme_minimal()

print(p_states)

# -----------------------------
# 2) Daily events (newE, newD, newR)
# -----------------------------
events_long <- melt(df,
                    id.vars = "date",
                    measure.vars = c("newE","newD","newR"),
                    variable.name = "event",
                    value.name = "count")

# --- ADD: observed daily confirmations (from cases) ---
obs_daily <- cases[!is.na(date_confirmed), .N, by = .(date = as.Date(date_confirmed))][order(date)]
setnames(obs_daily, "N", "obs_confirmed")

# Merge observed confirmations into events_long as another "event" series
obs_long <- data.table(date = obs_daily$date,
                       event = "obs_confirmed",
                       count = obs_daily$obs_confirmed)

events_long2 <- rbindlist(list(events_long, obs_long), use.names = TRUE, fill = TRUE)

p_events <- ggplot(events_long2, aes(x = date, y = count, group = event, color = event)) +
  geom_line(linewidth = 1) +
  labs(x = "Date", y = "Count", title = "Daily events + observed confirmations") +
  theme_minimal()

print(p_events)

# -----------------------------
# 3) Observed confirmations vs simulated detections (newD)
# -----------------------------
obs <- cases[!is.na(date_confirmed), .N, by = date_confirmed][order(date_confirmed)]
setnames(obs, "date_confirmed", "date")

pred <- sim$out_daily[, .(date, pred = newD)]
obs[,  date := as.Date(date)]
pred[, date := as.Date(date)]

cmp <- merge(obs, pred, by = "date", all = TRUE)
cmp[is.na(N),    N := 0]
cmp[is.na(pred), pred := 0]

p_cmp <- ggplot(cmp, aes(x = date)) +
  geom_line(aes(y = N,    color = "Observed confirmations"), linewidth = 1) +
  geom_line(aes(y = pred, color = "Simulated detections (newD)"), linewidth = 1) +
  labs(x = "Date", y = "Daily count", title = "Observed vs simulated") +
  theme_minimal() +
  scale_color_discrete(name = "")

print(p_cmp)

# -----------------------------
# 4) Base R quick plot of daily events
# -----------------------------
plot(sim$out_daily$date, sim$out_daily$newE, type = "l",
     xlab = "Date", ylab = "Count", main = "Daily events (base R)")
lines(sim$out_daily$date, sim$out_daily$newD)
lines(sim$out_daily$date, sim$out_daily$newR)
legend("topright", legend = c("newE","newD","newR"), lty = 1, bty = "n")

# -----------------------------
# 5) Map: observed confirmed outbreak farms (points)
# -----------------------------
obs_farms <- unique(as.character(cases[!is.na(date_confirmed), farm_id]))
obs_sf <- pop_sf[pop_sf$farm_id %in% obs_farms, ]

p_obs_map <- ggplot() +
  geom_sf(data = counties, fill = NA) +
  geom_sf(data = obs_sf, size = 2) +
  labs(title = "Observed confirmed outbreak farms") +
  theme_minimal()

print(p_obs_map)

# -----------------------------
# 6) Map: simulated end-of-period state (points)
# -----------------------------
state_lab <- c("S","E","I","D","R")
pop_sf$state_end <- factor(state_lab[sim$final_state], levels = state_lab)

p_endstate_map <- ggplot() +
  geom_sf(data = counties, fill = NA) +
  geom_sf(data = pop_sf, aes(color = state_end), size = 1.5, alpha = 0.8) +
  labs(title = "Simulated farm state at end of simulation", color = "State") +
  theme_minimal()

print(p_endstate_map)

# -----------------------------
# 7) County choropleth: simulated removed farms (R) by county
#    Requires: end_tab has columns county, state_end
# -----------------------------
if (!exists("end_tab")) {
  # If you don't have end_tab yet, create it from pop_sf
  end_tab <- as.data.table(st_drop_geometry(pop_sf))
}

stopifnot(all(c("county","state_end") %in% names(end_tab)))
stopifnot("county" %in% names(counties))

county_counts <- end_tab[, .(n_R = sum(state_end == "R", na.rm = TRUE)), by = county]

counties_map <- counties %>%
  left_join(county_counts, by = "county") %>%
  mutate(n_R = ifelse(is.na(n_R), 0, n_R))

p_county_R <- ggplot(counties_map) +
  geom_sf(aes(fill = n_R), color = "grey40", linewidth = 0.2) +
  labs(title = "Simulated removed farms by county (end date)", fill = "Count (R)") +
  theme_minimal()

print(p_county_R)
