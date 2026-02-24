# ============================================================
# HPAI Phase 1 — Q1 Script (Descriptive epidemiology)
# Outputs:
#  1) Tables: outbreak distribution by species & production
#  2) Timeline: daily + weekly incidence (overall + by production)
#  3) Maps: farms + outbreak farms + HRZ overlay (and optional county choropleth if county name field matches)
#  4) CSV exports + PNG figures
#
# NOTE:
# - This script is OBSERVED-DATA ONLY (no simulation). Perfect for Q1.
# - Requires: population.csv, cases.csv, activity.csv, prev_culls.csv, movement.csv (not all used in Q1),
#            hrz_32626.geojson, clc_32626.geojson (optional for Q1), counties_32626.geojson
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(sf)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# 0) Helpers
# -----------------------------
stop_if_missing <- function(df, cols, df_name = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) stop(df_name, " missing columns: ", paste(missing, collapse = ", "))
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

guess_first_existing <- function(candidates) {
  candidates <- candidates[candidates %in% names(candidates)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

guess_col <- function(df, candidates) {
  candidates <- candidates[candidates %in% names(df)]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

# -----------------------------
# 1) Paths + read inputs
# -----------------------------
DATA_DIR <- "~/Desktop/hpai_challenge/Phase 1 documents"
OUT_DIR  <- file.path(DATA_DIR, "Q1_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

csv_files <- list(
  pop   = "population.csv",
  cases = "cases.csv",
  act   = "activity.csv",
  prev  = "prev_culls.csv",
  mov   = "movement.csv"
)

geo_files <- list(
  hrz      = "hrz_32626.geojson",
  clc      = "clc_32626.geojson",        # optional for Q1
  counties = "counties_32626.geojson"
)

# Read CSV
pop   <- fread(file.path(DATA_DIR, csv_files$pop))
cases <- fread(file.path(DATA_DIR, csv_files$cases))

# Optional CSVs (not strictly needed for Q1, but often useful to sanity-check)
act_path  <- file.path(DATA_DIR, csv_files$act)
prev_path <- file.path(DATA_DIR, csv_files$prev)
mov_path  <- file.path(DATA_DIR, csv_files$mov)

act  <- if (file.exists(act_path))  fread(act_path)  else NULL
prev <- if (file.exists(prev_path)) fread(prev_path) else NULL
mov  <- if (file.exists(mov_path))  fread(mov_path)  else NULL

# Read GeoJSON
hrz_sf      <- st_read(file.path(DATA_DIR, geo_files$hrz), quiet = TRUE)
counties_sf <- st_read(file.path(DATA_DIR, geo_files$counties), quiet = TRUE)

# Optional: CLC (only needed if you want an environmental layer for Q1 maps; not required)
clc_path <- file.path(DATA_DIR, geo_files$clc)
clc_sf   <- if (file.exists(clc_path)) st_read(clc_path, quiet = TRUE) else NULL

# -----------------------------
# 2) Standardize columns
# -----------------------------
# population.csv
stop_if_missing(pop, c("farm_id", "x", "y", "production"), "population.csv")
pop[, farm_id := as.character(farm_id)]
pop[, production := norm_str(production)]

# Identify species column if present
species_col <- guess_col(pop, c("species", "Species", "SPECIES", "sp", "animal", "Animal", "ANIMAL"))
if (!is.na(species_col)) {
  pop[, species := norm_str(get(species_col))]
} else {
  pop[, species := NA_character_]
}

# cases.csv
stop_if_missing(cases, c("farm_id", "date_confirmed"), "cases.csv")
cases[, farm_id := as.character(farm_id)]

date_cols_cases <- intersect(c("date_suspicious", "date_confirmed", "cull_start", "cull_end"), names(cases))
for (cc in date_cols_cases) cases[[cc]] <- as_date_safe(cases[[cc]])

# Keep only confirmed outbreaks with dates
cases_conf <- cases[!is.na(date_confirmed), .(farm_id, date_confirmed)]
if (nrow(cases_conf) == 0) stop("cases.csv has no non-NA date_confirmed; cannot build Q1 outputs.")

# -----------------------------
# 3) Build sf points for farms (EPSG:32626)
# -----------------------------
pop_sf <- st_as_sf(pop, coords = c("x", "y"), crs = 32626, remove = FALSE)

# Ensure polygons are in same CRS
if (!is.na(st_crs(hrz_sf)$epsg) && st_crs(hrz_sf)$epsg != 32626) hrz_sf <- st_transform(hrz_sf, 32626)
if (!is.na(st_crs(counties_sf)$epsg) && st_crs(counties_sf)$epsg != 32626) counties_sf <- st_transform(counties_sf, 32626)

# -----------------------------
# 4) Join counties to farms (to support county summaries/maps)
# -----------------------------
# 4) Join counties to farms (to support county summaries/maps)
# -----------------------------
# Spatial join farms -> counties (adds county polygon attributes onto points)
pop_join <- st_join(pop_sf, counties_sf, join = st_within, left = TRUE)
# rename county.y to county
pop_join <- pop_join %>%
  dplyr::rename(county = county.y)

# Drop geometry + convert to data.table for safe keyed merge
pop_lookup <- data.table::as.data.table(sf::st_drop_geometry(pop_join))
pop_lookup[, farm_id := as.character(farm_id)]

# IMPORTANT: in your counties file, the real county field is "county"
if (!("county" %in% names(pop_lookup))) {
  stop("After st_join, column 'county' not found. Available: ", paste(names(pop_lookup), collapse = ", "))
}

# farm_id -> county lookup
county_lookup <- unique(pop_lookup[, .(
  farm_id,
  county = as.character(county)  # <-- real county name (Fulton, Berks, etc.)
)])

# Merge onto pop by farm_id (robust, order-independent)
pop <- merge(pop, county_lookup, by = "farm_id", all.x = TRUE)

# Add HRZ membership (binary) for mapping/summaries
pop$H_i <- as.integer(lengths(st_intersects(pop_sf, hrz_sf)) > 0)
pop <- pop %>%
  dplyr::rename(county = county.y)
# -----------------------------
# 5) Q1-A Table: outbreaks by species & production
# -----------------------------
cases_meta <- merge(
  cases_conf,
  pop[, .(farm_id, production, species, county, H_i)],
  by = "farm_id",
  all.x = TRUE
)

# Unique outbreak farms (some farms can appear multiple times in cases)
outbreak_farms <- unique(cases_meta$farm_id)

tab_prod <- cases_meta[, .(n_outbreak_farms = uniqueN(farm_id)), by = .(production)][order(-n_outbreak_farms)]
tab_sp_prod <- cases_meta[, .(n_outbreak_farms = uniqueN(farm_id)), by = .(species, production)][order(-n_outbreak_farms)]

# Optional: HRZ vs non-HRZ outbreaks
tab_hrz <- cases_meta[, .(n_outbreak_farms = uniqueN(farm_id)), by = .(H_i)][order(H_i)]

# Save tables
fwrite(tab_prod,    file.path(OUT_DIR, "Q1_table_outbreaks_by_production.csv"))
fwrite(tab_sp_prod, file.path(OUT_DIR, "Q1_table_outbreaks_by_species_and_production.csv"))
fwrite(tab_hrz,     file.path(OUT_DIR, "Q1_table_outbreaks_by_HRZ.csv"))

# Print to console for quick view
cat("\n=== Outbreak farms by production ===\n")
print(tab_prod)
cat("\n=== Outbreak farms by species x production ===\n")
print(tab_sp_prod)
cat("\n=== Outbreak farms by HRZ membership (0/1) ===\n")
print(tab_hrz)

# -----------------------------
# 6) Q1-B Timeline: daily + weekly incidence (overall + by production)
# -----------------------------
# Overall daily
obs_daily <- cases_conf[, .N, by = .(date = as.Date(date_confirmed))][order(date)]

p_inc_daily <- ggplot(obs_daily, aes(x = date, y = N)) +
  geom_col() +
  labs(
    x = "Date",
    y = "Confirmed outbreaks (farms)",
    title = "Observed outbreak incidence (daily)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_timeline_daily_incidence.png"),
  plot = p_inc_daily,
  width = 8, height = 4, dpi = 300, bg = "white"
)

# Overall weekly
obs_weekly <- copy(obs_daily)
obs_weekly[, week := floor_date(date, unit = "week", week_start = 1)]
obs_weekly <- obs_weekly[, .(N = sum(N)), by = week][order(week)]

p_inc_weekly <- ggplot(obs_weekly, aes(x = week, y = N)) +
  geom_col() +
  labs(
    x = "Week (Mon-start)",
    y = "Confirmed outbreaks (farms)",
    title = "Observed outbreak incidence (weekly)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_timeline_weekly_incidence.png"),
  plot = p_inc_weekly,
  width = 8, height = 4, dpi = 300, bg = "white"
)

# By production (daily)
obs_daily_prod <- cases_meta[, .N, by = .(date = as.Date(date_confirmed), production)][order(date)]

p_inc_daily_prod <- ggplot(obs_daily_prod, aes(x = date, y = N)) +
  geom_col() +
  facet_wrap(~ production, scales = "free_y") +
  labs(
    x = "Date",
    y = "Confirmed outbreaks (farms)",
    title = "Observed outbreak incidence (daily) by production"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_timeline_daily_incidence_by_production.png"),
  plot = p_inc_daily_prod,
  width = 10, height = 7, dpi = 300, bg = "white"
)

# Export timeline data
fwrite(obs_daily,      file.path(OUT_DIR, "Q1_timeseries_observed_daily.csv"))
fwrite(obs_weekly,     file.path(OUT_DIR, "Q1_timeseries_observed_weekly.csv"))
fwrite(obs_daily_prod, file.path(OUT_DIR, "Q1_timeseries_observed_daily_by_production.csv"))

# -----------------------------
# 7) Q1-C Spatial distribution maps
# -----------------------------
# Build outbreak sf
obs_sf <- pop_sf[pop_sf$farm_id %in% outbreak_farms, ]

# Map 1: All farms (faint) + outbreak farms + HRZ + counties
p_map_outbreaks <- ggplot() +
  geom_sf(data = counties_sf, fill = NA, linewidth = 0.2) +
  geom_sf(data = hrz_sf, fill = NA, linewidth = 0.6) +
  geom_sf(data = pop_sf, alpha = 0.15, size = 0.6) +
  geom_sf(data = obs_sf, size = 2) +
  labs(title = "Spatial distribution of farms and confirmed outbreaks") +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_map_outbreak_farms.png"),
  plot = p_map_outbreaks,
  width = 8, height = 7, dpi = 300, bg = "white"
)

# Map 2 (optional): outbreaks colored by production type
obs_sf2 <- obs_sf
obs_sf2$production <- pop$production[match(obs_sf2$farm_id, pop$farm_id)]

p_map_prod <- ggplot() +
  geom_sf(data = counties_sf, fill = NA, linewidth = 0.2) +
  geom_sf(data = hrz_sf, fill = NA, linewidth = 0.6) +
  geom_sf(data = obs_sf2, aes(color = production), size = 2) +
  labs(title = "Confirmed outbreaks by production type", color = "Production") +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_map_outbreak_by_production.png"),
  plot = p_map_prod,
  width = 8, height = 7, dpi = 300, bg = "white"
)

# -----------------------------
# 8) Optional: County-level choropleth (outbreak farms per county)

# Count unique outbreak farms per county
county_out <- cases_meta[!is.na(county),
                         .(n_outbreak_farms = uniqueN(farm_id)),
                         by = county][order(-n_outbreak_farms)]

# Join counts onto county polygons by the TRUE county field
counties_map2 <- counties_sf %>%
  left_join(county_out, by = "county") %>%
  mutate(n_outbreak_farms = ifelse(is.na(n_outbreak_farms), 0, n_outbreak_farms))

p_choro <- ggplot(counties_map2) +
  geom_sf(aes(fill = n_outbreak_farms), color = "grey40", linewidth = 0.2) +
  labs(
    title = "Confirmed outbreak farms by county (cumulative)",
    fill = "Outbreak farms"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(OUT_DIR, "Q1_map_county_choropleth_outbreaks.png"),
  plot = p_choro,
  width = 8, height = 7, dpi = 300, bg = "white"
)

# Export county risk table
fwrite(county_out, file.path(OUT_DIR, "Q1_table_outbreak_farms_by_county.csv"))

# -----------------------------
# 9) Export a “farm metadata” file used in Q1 (handy for later Q2/Q3)
# -----------------------------
farm_meta <- pop[, .(farm_id, production, species, county, H_i, x, y)]
fwrite(farm_meta, file.path(OUT_DIR, "Q1_farm_metadata.csv"))

# -----------------------------
# 10) Quick console summary
# -----------------------------
cat("\n====================================\n")
cat("Q1 outputs written to:\n", OUT_DIR, "\n")
cat("Key files:\n")
cat(" - Q1_table_outbreaks_by_production.csv\n")
cat(" - Q1_table_outbreaks_by_species_and_production.csv\n")
cat(" - Q1_timeseries_observed_daily.csv\n")
cat(" - Q1_timeline_daily_incidence.png\n")
cat(" - Q1_map_outbreak_farms.png\n")
cat("====================================\n")

list.files(file.path(DATA_DIR, "Q1_outputs"), pattern = "Q1_table", full.names = TRUE)
list.files(file.path(DATA_DIR, "Q1_outputs"), pattern = "timeline", full.names = TRUE)
list.files(file.path(DATA_DIR, "Q1_outputs"), pattern = "Q1_map", full.names = TRUE)
