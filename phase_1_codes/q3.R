# ============================================================
# Q3: Relative contribution galliform vs palmiped
# - Uses FOI decomposition + knockouts
# - Uses SAME spatial mechanics as Q2 (edge list), but splits by source group
# ============================================================

# Classify farms by production/species labels
is_galliform <- function(prod) grepl("broiler|layer|chicken|turkey|hen|pullet|breeder|gall", prod, ignore.case = TRUE)
is_palmiped  <- function(prod) grepl("duck|geese|goose|palm", prod, ignore.case = TRUE)

pop[, group := fifelse(is_galliform(production), "galliform",
                       fifelse(is_palmiped(production),  "palmiped",  "other"))]

# Q3 simulator: records FOI attribution (expected infections) daily
simulate_seidr_q3 <- function(start_date, end_date,
                              pop, farm_ids, id_to_idx,
                              mov_by_date, prev_by_date,
                              w_size,
                              params,
                              edge_i, edge_j, edge_d, edge_w,
                              seed_farms = character(0),
                              init_state = NULL,
                              force_observed = FALSE,
                              cases = NULL,
                              rng_seed = NULL,
                              knockout_group = c("none","galliform","palmiped")) {
  
  knockout_group <- match.arg(knockout_group)
  if (!is.null(rng_seed)) set.seed(rng_seed)
  
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "day")
  TT <- length(dates)
  N <- length(farm_ids)
  
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
  empty_mv <- data.table(src_idx=integer(), dst_idx=integer(), volume=numeric())
  mv_queue <- replicate(L_mv, empty_mv, simplify = FALSE)
  
  gate_progression_by_activity <- isTRUE(params$gate_progression_by_activity)
  gate_detection_by_activity   <- isTRUE(params$gate_detection_by_activity)
  gate_removal_by_activity     <- isTRUE(params$gate_removal_by_activity)
  
  out <- data.table(
    date = dates,
    newE = integer(TT), newI = integer(TT), newD = integer(TT), newR = integer(TT),
    S = integer(TT), E = integer(TT), I = integer(TT), D = integer(TT), R = integer(TT)
  )
  
  attr <- data.table(
    date = dates,
    expInf_gal_sp = numeric(TT),
    expInf_pal_sp = numeric(TT),
    expInf_gal_mv = numeric(TT),
    expInf_pal_mv = numeric(TT),
    expInf_bg     = numeric(TT),
    expInf_total  = numeric(TT)
  )
  
  is_gal <- (pop$group == "galliform")
  is_pal <- (pop$group == "palmiped")
  
  inf_w  <- numeric(N)
  lambda_bg <- numeric(N)
  
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
    
    # infectiousness
    inf_w[] <- 0.0
    inf_w[X == I] <- 1.0
    inf_w[X == D] <- params$eta
    inf_w[A == 0L] <- 0.0
    
    # knockouts: stop that group from transmitting
    if (knockout_group == "galliform") inf_w[is_gal] <- 0.0
    if (knockout_group == "palmiped")  inf_w[is_pal] <- 0.0
    
    # background
    lambda_bg[] <- params$beta_hrz * pop$H_i + params$beta_env * pop$Z_i
    
    # spatial split by SOURCE group using edge list (i -> j)
    lambda_sp_gal <- numeric(N)
    lambda_sp_pal <- numeric(N)
    
    iw_edge <- inf_w[edge_i]
    valid <- iw_edge > 0
    if (any(valid)) {
      e_i <- edge_i[valid]
      e_j <- edge_j[valid]
      e_d <- edge_d[valid]
      e_w <- edge_w[valid]
      iw  <- iw_edge[valid]
      
      base_contrib <- iw * e_w * K_exp(e_d, params$alpha)
      
      gal_mask <- is_gal[e_i]
      pal_mask <- is_pal[e_i]
      
      if (any(gal_mask)) {
        rs_g <- rowsum(base_contrib[gal_mask], e_j[gal_mask], reorder = FALSE)
        lambda_sp_gal[as.integer(rownames(rs_g))] <- params$beta_sp * rs_g[,1]
      }
      if (any(pal_mask)) {
        rs_p <- rowsum(base_contrib[pal_mask], e_j[pal_mask], reorder = FALSE)
        lambda_sp_pal[as.integer(rownames(rs_p))] <- params$beta_sp * rs_p[,1]
      }
    }
    
    # movement split by SOURCE group (only broiler_2 recipients)
    lambda_mv_gal <- numeric(N)
    lambda_mv_pal <- numeric(N)
    
    if (nrow(mv_recent) > 0) {
      sus_mask <- (X == S & A == 1L)
      mv2 <- mv_recent[dst_idx %in% which(sus_mask)]
      if (nrow(mv2) > 0) {
        mv2[, src_inf := inf_w[src_idx]]
        mv2 <- mv2[src_inf > 0]
        if (nrow(mv2) > 0) {
          is_b2 <- (pop$production == "broiler_2")
          mv2 <- mv2[is_b2[dst_idx]]
          if (nrow(mv2) > 0) {
            mv2[, contrib := g_vol(volume) * src_inf]
            
            mv2_gal <- mv2[is_gal[src_idx]]
            mv2_pal <- mv2[is_pal[src_idx]]
            
            if (nrow(mv2_gal) > 0) {
              mv_sum_g <- mv2_gal[, .(mv_contrib = sum(contrib)), by = dst_idx]
              lambda_mv_gal[mv_sum_g$dst_idx] <- params$beta_mv * mv_sum_g$mv_contrib
            }
            if (nrow(mv2_pal) > 0) {
              mv_sum_p <- mv2_pal[, .(mv_contrib = sum(contrib)), by = dst_idx]
              lambda_mv_pal[mv_sum_p$dst_idx] <- params$beta_mv * mv_sum_p$mv_contrib
            }
          }
        }
      }
    }
    
    lambda_total <- pmax(lambda_bg, 0) + pmax(lambda_sp_gal, 0) + pmax(lambda_sp_pal, 0) +
      pmax(lambda_mv_gal, 0) + pmax(lambda_mv_pal, 0)
    
    # S->E
    p_inf <- 1 - exp(-pmax(lambda_total, 0))
    newE_vec <- rbinom(N, 1, p_inf)
    newE_vec[!(X == S & A == 1L)] <- 0L
    X[newE_vec == 1L] <- E
    
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
    
    # Expected attributable infections (sum over farms, by FOI fraction)
    lt <- lambda_total
    denom <- ifelse(lt > 0, lt, NA_real_)
    
    frac_gal_sp <- lambda_sp_gal / denom
    frac_pal_sp <- lambda_sp_pal / denom
    frac_gal_mv <- lambda_mv_gal / denom
    frac_pal_mv <- lambda_mv_pal / denom
    frac_bg     <- lambda_bg     / denom
    
    attr[tt, `:=`(
      expInf_gal_sp = sum(p_inf * ifelse(is.finite(frac_gal_sp), frac_gal_sp, 0), na.rm = TRUE),
      expInf_pal_sp = sum(p_inf * ifelse(is.finite(frac_pal_sp), frac_pal_sp, 0), na.rm = TRUE),
      expInf_gal_mv = sum(p_inf * ifelse(is.finite(frac_gal_mv), frac_gal_mv, 0), na.rm = TRUE),
      expInf_pal_mv = sum(p_inf * ifelse(is.finite(frac_pal_mv), frac_pal_mv, 0), na.rm = TRUE),
      expInf_bg     = sum(p_inf * ifelse(is.finite(frac_bg),     frac_bg,     0), na.rm = TRUE),
      expInf_total  = sum(p_inf, na.rm = TRUE)
    )]
  }
  
  list(out_daily = out, attr_daily = attr, final_state = X)
}

run_one_q3_forecast <- function(sim_id, rng_seed, knockout_group = "none") {
  assim <- simulate_seidr_q3(
    start_date = start_sim,
    end_date   = today,
    pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
    mov_by_date = mov_by_date, prev_by_date = prev_by_date,
    w_size = w_size,
    params = params,
    edge_i = edge_i, edge_j = edge_j, edge_d = edge_d, edge_w = edge_w,
    seed_farms = seed,
    init_state = NULL,
    force_observed = TRUE,
    cases = cases,
    rng_seed = rng_seed,
    knockout_group = "none"   # always none during assimilation
  )
  
  init_state <- assim$final_state
  
  fc <- simulate_seidr_q3(
    start_date = today + 1,
    end_date   = today + 28,
    pop = pop, farm_ids = farm_ids, id_to_idx = id_to_idx,
    mov_by_date = mov_by_date, prev_by_date = prev_by_date,
    w_size = w_size,
    params = params,
    edge_i = edge_i, edge_j = edge_j, edge_d = edge_d, edge_w = edge_w,
    seed_farms = character(0),
    init_state = init_state,
    force_observed = FALSE,
    cases = NULL,
    rng_seed = rng_seed + 9999,
    knockout_group = knockout_group
  )
  
  fc$out_daily[, sim := sim_id]
  fc$attr_daily[, sim := sim_id]
  list(traj = fc$out_daily, attr = fc$attr_daily)
}

# Choose Q3 sims (start smaller; 300 can be heavy)
n_sims_q3 <- 20

# baseline
set.seed(321)
res_base <- vector("list", n_sims_q3)
for (s in seq_len(n_sims_q3)) {
  res_base[[s]] <- run_one_q3_forecast(sim_id = s, rng_seed = 20000 + s, knockout_group = "none")
  if (s %% 25 == 0) cat("Q3 baseline sims:", s, " / ", n_sims_q3, "\n")
}
traj_base <- rbindlist(lapply(res_base, `[[`, "traj"))
attr_base <- rbindlist(lapply(res_base, `[[`, "attr"))

# knockouts
res_no_gal <- vector("list", n_sims_q3)
for (s in seq_len(n_sims_q3)) {
  res_no_gal[[s]] <- run_one_q3_forecast(sim_id = s, rng_seed = 30000 + s, knockout_group = "galliform")
  if (s %% 25 == 0) cat("Q3 no-galliform sims:", s, " / ", n_sims_q3, "\n")
}
traj_no_gal <- rbindlist(lapply(res_no_gal, `[[`, "traj"))

res_no_pal <- vector("list", n_sims_q3)
for (s in seq_len(n_sims_q3)) {
  res_no_pal[[s]] <- run_one_q3_forecast(sim_id = s, rng_seed = 40000 + s, knockout_group = "palmiped")
  if (s %% 25 == 0) cat("Q3 no-palmiped sims:", s, " / ", n_sims_q3, "\n")
}
traj_no_pal <- rbindlist(lapply(res_no_pal, `[[`, "traj"))

# Summarize attribution over 4-week forecast window
attr_tot <- attr_base[
  , .(
    gal_sp = sum(expInf_gal_sp),
    pal_sp = sum(expInf_pal_sp),
    gal_mv = sum(expInf_gal_mv),
    pal_mv = sum(expInf_pal_mv),
    bg     = sum(expInf_bg),
    total  = sum(expInf_total)
  ),
  by = sim
]

attr_tot[, `:=`(
  gal_total = gal_sp + gal_mv,
  pal_total = pal_sp + pal_mv,
  frac_gal  = fifelse(total > 0, (gal_sp + gal_mv)/total, NA_real_),
  frac_pal  = fifelse(total > 0, (pal_sp + pal_mv)/total, NA_real_)
)]

summary_attr <- attr_tot[, .(
  frac_gal_q05 = quantile(frac_gal, 0.05, na.rm=TRUE),
  frac_gal_q50 = quantile(frac_gal, 0.50, na.rm=TRUE),
  frac_gal_q95 = quantile(frac_gal, 0.95, na.rm=TRUE),
  frac_pal_q05 = quantile(frac_pal, 0.05, na.rm=TRUE),
  frac_pal_q50 = quantile(frac_pal, 0.50, na.rm=TRUE),
  frac_pal_q95 = quantile(frac_pal, 0.95, na.rm=TRUE)
)]
print(summary_attr)

fwrite(attr_base, file.path(OUT_Q3, "Q3_attr_daily_all_sims.csv"))
fwrite(attr_tot,  file.path(OUT_Q3, "Q3_attr_4week_totals_by_sim.csv"))
fwrite(summary_attr, file.path(OUT_Q3, "Q3_attr_4week_fraction_summary.csv"))

attr_path <- attr_tot[, .(
  frac_gal_from_spatial = gal_sp / pmax(gal_sp + gal_mv, 1e-12),
  frac_pal_from_spatial = pal_sp / pmax(pal_sp + pal_mv, 1e-12)
)]
fwrite(attr_path, file.path(OUT_Q3, "Q3_attr_pathway_shares_by_sim.csv"))

# Knockout impact
sum4 <- function(traj_dt) {
  traj_dt[, .(newE_4w = sum(newE), newD_4w = sum(newD)), by = sim]
}

base4  <- sum4(traj_base)
nogal4 <- sum4(traj_no_gal)
nopal4 <- sum4(traj_no_pal)

cmp <- Reduce(function(x,y) merge(x,y, by="sim", all=TRUE),
              list(
                base4,
                setnames(copy(nogal4), c("sim","newE_4w","newD_4w"), c("sim","newE_4w_no_gal","newD_4w_no_gal")),
                setnames(copy(nopal4), c("sim","newE_4w","newD_4w"), c("sim","newE_4w_no_pal","newD_4w_no_pal"))
              ))

cmp[, `:=`(
  red_newE_if_remove_gal = (newE_4w - newE_4w_no_gal) / pmax(newE_4w, 1),
  red_newE_if_remove_pal = (newE_4w - newE_4w_no_pal) / pmax(newE_4w, 1),
  red_newD_if_remove_gal = (newD_4w - newD_4w_no_gal) / pmax(newD_4w, 1),
  red_newD_if_remove_pal = (newD_4w - newD_4w_no_pal) / pmax(newD_4w, 1)
)]

summary_knock <- cmp[, .(
  red_newE_gal_q50 = quantile(red_newE_if_remove_gal, 0.50, na.rm=TRUE),
  red_newE_gal_q05 = quantile(red_newE_if_remove_gal, 0.05, na.rm=TRUE),
  red_newE_gal_q95 = quantile(red_newE_if_remove_gal, 0.95, na.rm=TRUE),
  
  red_newE_pal_q50 = quantile(red_newE_if_remove_pal, 0.50, na.rm=TRUE),
  red_newE_pal_q05 = quantile(red_newE_if_remove_pal, 0.05, na.rm=TRUE),
  red_newE_pal_q95 = quantile(red_newE_if_remove_pal, 0.95, na.rm=TRUE),
  
  red_newD_gal_q50 = quantile(red_newD_if_remove_gal, 0.50, na.rm=TRUE),
  red_newD_gal_q05 = quantile(red_newD_if_remove_gal, 0.05, na.rm=TRUE),
  red_newD_gal_q95 = quantile(red_newD_if_remove_gal, 0.95, na.rm=TRUE),
  
  red_newD_pal_q50 = quantile(red_newD_if_remove_pal, 0.50, na.rm=TRUE),
  red_newD_pal_q05 = quantile(red_newD_if_remove_pal, 0.05, na.rm=TRUE),
  red_newD_pal_q95 = quantile(red_newD_if_remove_pal, 0.95, na.rm=TRUE)
)]
print(summary_knock)

fwrite(cmp, file.path(OUT_Q3, "Q3_knockout_impacts_by_sim.csv"))
fwrite(summary_knock, file.path(OUT_Q3, "Q3_knockout_impacts_summary.csv"))

# Daily fraction plot (median, IQR) for galliform attribution
attr_daily_sum <- attr_base[
  , .(
    gal = expInf_gal_sp + expInf_gal_mv,
    pal = expInf_pal_sp + expInf_pal_mv,
    total = expInf_total
  ),
  by = .(sim, date)
]
attr_daily_sum[, `:=`(
  frac_gal = fifelse(total > 0, gal/total, NA_real_),
  frac_pal = fifelse(total > 0, pal/total, NA_real_)
)]

daily_med <- attr_daily_sum[
  , .(
    frac_gal_q50 = quantile(frac_gal, 0.50, na.rm=TRUE),
    frac_gal_q25 = quantile(frac_gal, 0.25, na.rm=TRUE),
    frac_gal_q75 = quantile(frac_gal, 0.75, na.rm=TRUE)
  ),
  by = date
]

p_frac <- ggplot(as.data.frame(daily_med), aes(x=date, y=frac_gal_q50)) +
  geom_ribbon(aes(ymin=frac_gal_q25, ymax=frac_gal_q75), alpha=0.25) +
  geom_line(linewidth=1) +
  labs(
    title="Q3: Fraction of infections attributable to galliform sources (median, IQR)",
    x="Date", y="Fraction attributable to galliform"
  ) +
  theme_minimal()

ggsave(file.path(OUT_Q3, "Q3_fraction_attributable_galliform_timeseries.png"),
       p_frac, width=8, height=4, dpi=300, bg="white")

cat("\n================ Q3 DONE ================\n")
cat("Outputs in: ", OUT_Q3, "\n")
cat("Attribution daily:   Q3_attr_daily_all_sims.csv\n")
cat("Attribution totals:  Q3_attr_4week_totals_by_sim.csv\n")
cat("Fraction summary:    Q3_attr_4week_fraction_summary.csv\n")
cat("Pathway shares:      Q3_attr_pathway_shares_by_sim.csv\n")
cat("Knockout impacts:    Q3_knockout_impacts_by_sim.csv\n")
cat("Knockout summary:    Q3_knockout_impacts_summary.csv\n")
cat("Figure:              Q3_fraction_attributable_galliform_timeseries.png\n")
cat("========================================\n")

