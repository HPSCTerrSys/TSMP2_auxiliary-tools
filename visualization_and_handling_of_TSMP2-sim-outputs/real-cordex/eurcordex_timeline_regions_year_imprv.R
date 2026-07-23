#!/usr/bin/env Rscript
#
# Regional time series of EURO-CORDEX / TSMP2 (ICON, eCLM, ParFlow) output.
#
# For each variable in `PLOTVARS` and each applicable experiment in `MODELS`
# (see select_models() — "saturation" excludes native ICON, which has no land
# model), reads hourly model output, computes regional (mean/25th/75th
# percentile) statistics per `REGIONS`, and writes three plots per
# region/variable: the full hourly series, a daily mean/min/max series, and a
# mean diurnal cycle. A separate section computes domain-wide relative soil
# saturation (ParFlow-based, daily cadence) from a smaller set of experiments.
#
# Usage:
#   Rscript eurcordex_timeline_regions_year_imprv.R
#
# All configuration (experiments, variables, date, paths, plot styling) is in
# the CONFIG section below.

suppressPackageStartupMessages({
  library(ncdf4)
  library(dplyr)
  library(abind)
  library(tibble)
})

########
### CONFIG
########

BASEDIR   <- "/p/scratch/cslts/poll1/sim/paper"
DIROUT    <- "/p/project1/cslts/poll1/tools/plots_tsmp2-year_tst/"
PLTTYPE   <- "png"                        # "png" or "X11"
DATE_SEL  <- as.Date("2018-01-01")        # only year is used to build file-search patterns
PLOTVARS  <- c("lhf") #c("t2m", "q2m", "lhf", "shf", "saturation") # add "ef" for evaporative fraction

# Overlay gridded observations (in black, with per-model RMSE) on the daily
# plots: E-OBS for t2m/q2m, GLEAM for lhf/shf, ESA CCI for saturation. Set to
# FALSE to disable.
PLOT_OBS <- TRUE

# Seasons (by calendar month) to additionally break the diurnal-cycle plot
# out by, on top of the whole-year diurnal plot. Set to list() to disable and
# only produce the whole-year diurnal plot.
DIURNAL_SEASONS <- list(
  DJF = c(12, 1, 2),
  MAM = c(3, 4, 5),
  JJA = c(6, 7, 8),
  SON = c(9, 10, 11)
)

# Soil layers (indices into eCLM's levsoi) to compute "saturation" for -- one
# full set of plots (regional timeseries and the ParFlow-native regional
# daily plot) is produced per level, with the level in each plot's filename.
SOIL_LEVELS <- c(1,3)

# Experiments used for the regional timeseries analysis (order matters: it
# fixes the column order of the result arrays and must match MODEL_STYLES).
MODELS <- c("eclm-pfl", "eclm", "icon", "icon-eclm", "icon-eclm-pfl") # c("icon", "eclm")
EXPERIMENTS <- file.path(BASEDIR, paste0("wfe_eur-11_revsetup_", MODELS))

STATIC <- list(
  icon_extpar_file = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/icon/static/external_parameter_icon_europe011_DOM01_tiles.nc",
  icon_grid_file  = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/icon/static/europe011_DOM01.nc",
  eclm_grid_file  = "/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/mkmapgrids/EUR-R13B05_189976_grid.nc",
  eclm_surf_file  = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/eclm/static/surfdata_ICON-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230302_gcvurb-pfsoil_halo.nc",
  parflow_indicator_file = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/parflow/static/EUR-11_TSMP_FZJ-IBG3_eCLMPFLDomain_444x432_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_alv.sa",
  parflow_grid_file = "/p/scratch/cslts/poll1/git/TSMP_EUR-11/static/clm/griddata_CLM_EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432.nc"
)

# E-OBS gridded observations (0.1 deg regular lat/lon), daily, 2011-2025.
# GLEAM (global, 0.1 deg, one file per year) for latent/sensible heat flux
# comparison. ESA CCI (EUR-11 rotated-pole grid, one file per year) for soil
# moisture / saturation comparison. "%s" is filled in with the DATE_SEL year.
OBS <- list(
  eobs_tg_file = "/p/data1/slts/shared_data/obs_THU_EOBS/o.data/tg_ens_mean_0.1deg_reg_2011-2025_v33.0e.nc",
  eobs_hu_file = "/p/data1/slts/shared_data/obs_THU_EOBS/o.data/hu_ens_mean_0.1deg_reg_2011-2025_v33.0e.nc",
  gleam_e_file = "/p/data1/slts/shared_data/obs_HLESM_GLEAM/o.data/E_%s_GLEAM_v4.3a.nc",
  gleam_h_file = "/p/data1/slts/shared_data/obs_HLESM_GLEAM/o.data/H_%s_GLEAM_v4.3a.nc",
  esacci_sm_file = "/p/data1/jibg31/DETECT/observations/ESACCI_SM/eu_cordex11/daily/merged/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-%s-fv08.1_EUCORDEX011.nc"
)

REGIONS <- tribble(
  ~region, ~xmin, ~xmax, ~ymin, ~ymax,
  "EU", -50, 70, 18, 75, # whole EUR-11
  "BI", -10,  2, 50, 59, # British Isles
  "IP", -10,  5, 36, 44, # Iberian Peninsula
  "FR",  -5,  5, 44, 50, # France
  "ME",   2, 16, 47, 55, # Mid-Europe
  "SC",   5, 30, 55, 71, # Scandinavia
  "AL",   5, 15, 44, 48, # Alps
  "MD",  -5, 25, 36, 44, # Mediterranean
  "EA",  15, 40, 44, 55, # Eastern Europe
  "CE", -50, 70, 36, 75  # EUR-11 without African continental Europe
)
RAD2DEG <- 180 / pi

# Label/color/line-type per model, looked up by name so plot styling always
# lines up with whichever subset/order of models an analysis actually uses.
MODEL_STYLES <- tribble(
  ~model,           ~label,               ~color,         ~linetype,
  "icon",           "ICON",               "dodgerblue3",  5,
  "eclm",           "eCLM",               "forestgreen",  3,
  "eclm-pfl",       "eCLM-ParFlow",       "olivedrab3",   2,
  "icon-eclm",      "ICON-eCLM",          "darkorange1",  4,
  "icon-eclm-pfl",  "ICON-eCLM-ParFlow",  "firebrick2",   6,
  "pfl",            "ParFlow",            "sienna",       1
)

########
### FUNCTIONS
########

# -- data loading -------------------------------------------------------

# collapse: pass FALSE for ParFlow output, which has exactly one timestep per
# file (a length-1 time dim) that ncvar_get() would otherwise silently drop,
# breaking the "combine along the last (time) dimension" abind below. Must
# stay TRUE (default) for ICON, which has a genuine singleton diagnostic
# level dimension on near-surface variables (t_2m etc.) that needs dropping.
#
# soil_level: for a variable with a "levsoi"/"levgrnd" dimension (H2OSOI),
# read only this one level via NetCDF start/count instead of loading every
# level into memory before subsetting -- H2OSOI has 20 levels, so loading
# them all for a full year of hourly output is ~20x more data than needed
# and will exhaust memory.
#
# eCLM's own "time" declares its units as "days since <this monthly chunk's
# own start date>", not one fixed reference for the whole run (confirmed via
# ncdump: a Jan-05 file and a Jun-05 file both report time = 4, each
# relative to its own month's 1st) -- read as a common absolute day count
# (days since 1970-01-01) instead, so concatenating time across chunks
# doesn't alias different months onto the same value (e.g. the day_tag
# grouping in run_soil_moisture_analysis below). ParFlow's "time" has no CF
# units attribute at all (a bare step index); passed through unchanged since
# it's never used to derive a calendar date. Native ICON's own "time" isn't
# "days since ..." either -- it's already an absolute encoded date/time
# (units "day as %Y%m%d.%f", e.g. 20180105.0417 for 01:00 on Jan 5), so this
# is a no-op for it too.
read_time_days <- function(ncfile) {
  raw <- ncvar_get(ncfile, "time")
  units_attr <- ncatt_get(ncfile, "time", "units")
  if (!units_attr$hasatt || !grepl("^days since", units_attr$value)) return(raw)
  origin <- as.POSIXct(sub("^days since\\s*", "", units_attr$value), tz = "UTC")
  raw + as.numeric(difftime(origin, as.POSIXct("1970-01-01", tz = "UTC"), units = "days"))
}

load_tsmp2_data <- function(vars, pattern, collapse = TRUE, soil_level = NULL) {
  files <- Sys.glob(pattern)
  if (length(files) == 0) stop("No matching files found for pattern: ", pattern, " var(s): ", paste(vars, collapse = ", "))

  data_lists <- setNames(vector("list", length(vars)), vars)
  time_list <- list()

  for (ii in seq_along(files)) {
    ncfile <- nc_open(files[ii])
    cat("Reading:", files[ii], "\n")

    for (vv in vars) {
      if (!vv %in% names(ncfile$var)) {
        warning("Variable '", vv, "' not found in: ", files[ii])
      } else {
        v <- ncfile$var[[vv]]
        # v$dim lists dimensions in ncdump's declared order, but ncvar_get()'s
        # start/count expect the reverse of that (the returned array's own,
        # fastest-varying-first order) -- e.g. dim 2 of 3 declared dims is
        # dim (3-2+1)=2 in that reversed order.
        lev_dim_declared <- which(vapply(v$dim, function(d) d$name, character(1)) %in% c("levsoi", "levgrnd"))
        if (!is.null(soil_level) && length(lev_dim_declared) == 1) {
          lev_dim <- v$ndims - lev_dim_declared + 1
          start <- rep(1, v$ndims); start[lev_dim] <- soil_level
          count <- rep(-1, v$ndims); count[lev_dim] <- 1
          data_lists[[vv]][[ii]] <- ncvar_get(ncfile, vv, start = start, count = count, collapse = collapse)
        } else {
          data_lists[[vv]][[ii]] <- ncvar_get(ncfile, vv, collapse = collapse)
        }
      }
    }
    time_list[[ii]] <- read_time_days(ncfile)
    nc_close(ncfile)
  }

  data_out <- lapply(data_lists, function(lst) {
    along_dim <- length(dim(lst[[length(lst)]]))
    do.call(abind, c(lst, list(along = along_dim)))
  })
  data_out$time <- unlist(time_list)
  data_out
}

# -- experiment/variable/file-pattern helpers ---------------------------

# eCLM-based runs (eclm, eclm-pfl, icon-eclm, icon-eclm-pfl) all have "eclm"
# in their directory name and share eCLM variable names/output layout.
# The ICON-only run is the sole experiment with "icon" but not "eclm".
classify_experiment <- function(experiment_dir) {
  name <- basename(experiment_dir)
  is_eclm <- grepl("eclm", name)
  is_icon_native <- !is_eclm && grepl("icon", name)
  if (!is_eclm && !is_icon_native) stop("Cannot classify experiment: ", experiment_dir)
  list(is_eclm = is_eclm, is_icon_native = is_icon_native)
}

# "saturation" expands (see PLOTVARS in CONFIG) into one "saturation_levN"
# plotvar per entry in SOIL_LEVELS; these two recognize/parse that.
is_saturation_pv <- function(pv) grepl("^saturation_lev", pv)
saturation_pv_level <- function(pv) as.integer(sub("^saturation_lev", "", pv))

get_var_names <- function(is_eclm, plotvar) {
  if (is_saturation_pv(plotvar)) {
    if (!is_eclm) stop("Unknown plotvar: ", plotvar) # native ICON has no soil column
    return("H2OSOI")
  }
  if (is_eclm) {
    switch(plotvar,
      t2m = "TSA", q2m = "Q2M", lhf = "EFLX_LH_TOT", shf = "FSH",
      ef  = c("FSH", "EFLX_LH_TOT"),
      stop("Unknown plotvar: ", plotvar))
  } else {
    switch(plotvar,
      t2m = "t_2m", q2m = "qv_2m", lhf = "lhfl_s", shf = "shfl_s",
      ef  = c("shfl_s", "lhfl_s"),
      stop("Unknown plotvar: ", plotvar))
  }
}

# "saturation" needs an eCLM soil column, which native ICON (no land model)
# doesn't have; every other model name contains "eclm" and has one.
select_models <- function(plotvar, models) {
  if (is_saturation_pv(plotvar)) models[models != "icon"] else models
}

# CLM/eCLM surface dataset -> porosity per gridcell/soil-level, used to turn
# volumetric soil water (H2OSOI) into relative saturation. This is a fallback
# for when eCLM's own WATSAT output isn't available (see
# get_eclm_porosity_at_level() below, which is what callers should actually
# use) -- it approximates CLM/CTSM's actual saturated water content, which
# blends mineral and organic porosity (see CTSM's
# SoilStateInitTimeConstMod.F90); the mineral (sand-only) term alone
# underestimates porosity in organic-rich soils (peatlands, boreal regions),
# which was pushing H2OSOI/porosity above 100% there.
soil_porosity <- function(surf_file) {
  ncfile <- nc_open(surf_file)
  on.exit(nc_close(ncfile))
  watsat_mineral <- 0.489 - 0.00126 * ncvar_get(ncfile, "PCT_SAND")
  om_frac <- pmin(ncvar_get(ncfile, "ORGANIC") / 130, 1)
  (1 - om_frac) * watsat_mineral + om_frac * 0.9
}

build_file_pattern <- function(experiment_dir, model, is_eclm, date_sel) {
  model_dirname <- gsub("pfl", "parflow", gsub("-", "", model))
  if (is_eclm) {
    datadir <- file.path(experiment_dir, "dta/simres", paste0(model_dirname, "_20*"), "out/eclm")
    file.path(datadir, paste0("eCLM_eur-12-iic.clm2.h0.", format(date_sel, "%Y"), "*0000.nc"))
  } else {
    datadir <- file.path(experiment_dir, "dta/simres", paste0(model_dirname, "_20*"), "out/icon")
    file.path(datadir, paste0("ICON_out_EU-R13B5_inst_DOM01_ML_", format(date_sel, "%Y"), "*T000000Z_1h.nc"))
  }
}

# eCLM writes its own exact saturated water content (WATSAT, on levgrnd) into
# the first history file of each monthly restart chunk -- prefer that over
# the PCT_SAND/ORGANIC approximation in soil_porosity() above, since it's the
# literal value CLM used internally, not an estimate of it, for this
# particular experiment. levsoi (used by H2OSOI/"level" everywhere else in
# this script) is always the top subset of levgrnd, so WATSAT[, level] lines
# up directly. Falls back to soil_porosity() if no available file has WATSAT.
get_eclm_porosity_at_level <- function(level, static, experiment_dir, model, date_sel) {
  pattern <- build_file_pattern(experiment_dir, model, is_eclm = TRUE, date_sel)
  for (f in sort(Sys.glob(pattern))) {
    ncfile <- nc_open(f)
    has_watsat <- "WATSAT" %in% names(ncfile$var)
    if (has_watsat) watsat <- ncvar_get(ncfile, "WATSAT")[, level]
    nc_close(ncfile)
    if (has_watsat) return(watsat)
  }
  soil_porosity(static$eclm_surf_file)[, level]
}

# -- regional masks -------------------------------------------------------

build_region_mask <- function(clon_deg, clat_deg, regions) {
  sapply(seq_len(nrow(regions)), function(rr) {
    clon_deg >= regions$xmin[rr] & clon_deg <= regions$xmax[rr] &
      clat_deg >= regions$ymin[rr] & clat_deg <= regions$ymax[rr]
  })
}

get_region_mask <- function(is_icon_native, regions, static) {
  grid_file <- if (is_icon_native) static$icon_grid_file else static$eclm_grid_file
  ncfile <- nc_open(grid_file)
  clon <- ncvar_get(ncfile, "clon")
  clat <- ncvar_get(ncfile, "clat")
  nc_close(ncfile)
  build_region_mask(clon * RAD2DEG, clat * RAD2DEG, regions)
}

# ParFlow's grid is a plain (x, y) raster; LONGXY/LATIXY are declared
# (lsmlat, lsmlon), which ncvar_get() reverses to R dim (lsmlon, lsmlat) =
# (x, y) -- the same x-fastest flatten order used elsewhere for ParFlow
# fields (e.g. dim(values) <- c(prod(orig_dim[1:2]), ...)), so no reordering
# is needed before comparing against REGIONS.
get_parflow_region_mask <- function(regions, static) {
  ncfile <- nc_open(static$parflow_grid_file)
  lon <- as.vector(ncvar_get(ncfile, "LONGXY"))
  lat <- as.vector(ncvar_get(ncfile, "LATIXY"))
  nc_close(ncfile)
  build_region_mask(lon, lat, regions)
}

get_land_mask <- function(static) {
  ncfile <- nc_open(static$icon_extpar_file)
  fr_land <- ncvar_get(ncfile, "FR_LAND")
  nc_close(ncfile)
  fr_land
}

# -- regional statistics -------------------------------------------------

compute_region_stats <- function(data_mat, mask) {
  ntime <- ncol(data_mat)
  nreg <- ncol(mask)
  mean_mat <- q25_mat <- q75_mat <- matrix(NA_real_, ntime, nreg)
  for (rr in seq_len(nreg)) {
    sub <- data_mat[mask[, rr], , drop = FALSE]
    mean_mat[, rr] <- colMeans(sub, na.rm = TRUE)
    q25_mat[, rr] <- apply(sub, 2, quantile, probs = 0.25, na.rm = TRUE)
    q75_mat[, rr] <- apply(sub, 2, quantile, probs = 0.75, na.rm = TRUE)
  }
  list(mean = mean_mat, q25 = q25_mat, q75 = q75_mat)
}

# -- observations (E-OBS / GLEAM) ------------------------------------------

# Standard-atmosphere barometric formula. Used to elevation-correct pressure
# per REGIONS box (from ICON's topography) when converting E-OBS relative
# humidity to specific humidity -- E-OBS has no pressure field of its own.
region_pressure_hpa <- function(elevation_m) {
  1013.25 * exp(-elevation_m / 8434.5)
}

# Mean ICON surface elevation per REGIONS box. Region-level (not per-gridcell)
# to match the level of precision used everywhere else in this script.
get_region_elevation <- function(regions, static) {
  ncfile <- nc_open(static$icon_extpar_file)
  clon <- ncvar_get(ncfile, "clon") * RAD2DEG
  clat <- ncvar_get(ncfile, "clat") * RAD2DEG
  elev <- ncvar_get(ncfile, "topography_c")
  nc_close(ncfile)
  mask <- build_region_mask(clon, clat, regions)
  vapply(seq_len(nrow(regions)), function(rr) mean(elev[mask[, rr]], na.rm = TRUE), numeric(1))
}

# Mean eCLM porosity (at the given soil level) per REGIONS box, used to
# convert ESA CCI's raw volumetric soil moisture into the same
# relative-saturation [%] units the "saturation" plotvar itself plots.
# Region-level, matching the rest of this script's level of precision (see
# get_region_elevation). Not tied to one particular experiment (ESA CCI is
# compared against every eclm-based model alike), so this just reads WATSAT
# from the first eclm-based experiment in MODELS/EXPERIMENTS that has it --
# a static soil property that doesn't vary between experiments sharing the
# same surface dataset.
get_region_porosity <- function(regions, static, level) {
  ncfile <- nc_open(static$eclm_grid_file)
  clon <- ncvar_get(ncfile, "clon") * RAD2DEG
  clat <- ncvar_get(ncfile, "clat") * RAD2DEG
  nc_close(ncfile)
  is_eclm_vec <- vapply(EXPERIMENTS, function(e) classify_experiment(e)$is_eclm, logical(1))
  porosity <- get_eclm_porosity_at_level(level, static, EXPERIMENTS[is_eclm_vec][1], MODELS[is_eclm_vec][1], DATE_SEL)
  mask <- build_region_mask(clon, clat, regions)
  vapply(seq_len(nrow(regions)), function(rr) mean(porosity[mask[, rr]], na.rm = TRUE), numeric(1))
}

# Magnus-formula saturation vapor pressure -> specific humidity (kg/kg) from
# relative humidity. t_celsius/rh_pct/p_hpa may be scalars or matching-shape
# vectors/matrices.
relhum_to_specific_humidity <- function(t_celsius, rh_pct, p_hpa) {
  es <- 6.1094 * exp(17.625 * t_celsius / (t_celsius + 243.04))
  e <- rh_pct / 100 * es
  0.622 * e / (p_hpa - 0.378 * e)
}

# GLEAM's actual evaporation (mm/day) -> latent heat flux (W/m^2), via the
# latent heat of vaporization (2.45 MJ/kg, the standard ~20 degC reference
# value also used e.g. in FAO-56 Penman-Monteith).
evap_to_latent_heat_flux <- function(e_mm_day) e_mm_day * 2.45e6 / 86400

# (start, count) NetCDF index range covering [lo, hi] for a monotonic
# (ascending or descending) coordinate vector.
grid_index_range <- function(coord, lo, hi) {
  in_range <- which(coord >= lo & coord <= hi)
  if (length(in_range) == 0) stop("No grid points in range [", lo, ", ", hi, "]")
  list(start = min(in_range), count = max(in_range) - min(in_range) + 1)
}

# Reads one regular-(lat,lon)-grid NetCDF variable for the given time index
# range, cropped to REGIONS' lon/lat extent (important for a global grid like
# GLEAM's -- avoids loading way more data than needed), and reduces it to a
# per-timestep, per-region mean.
load_gridded_region_means <- function(file, varname, lon_name, lat_name, regions, time_start, time_count) {
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, lon_name)
  lat <- ncvar_get(ncfile, lat_name)

  lon_idx <- grid_index_range(lon, min(regions$xmin), max(regions$xmax))
  lat_idx <- grid_index_range(lat, min(regions$ymin), max(regions$ymax))

  # ncvar_get()'s start/count are in the order of the *returned* R array
  # (fastest-varying dimension first), which is the reverse of ncdump's
  # declared order -- for a (time, lat, lon) variable that's (lon, lat, time).
  values <- ncvar_get(ncfile, varname,
                       start = c(lon_idx$start, lat_idx$start, time_start),
                       count = c(lon_idx$count, lat_idx$count, time_count))
  nc_close(ncfile)

  lon_sub <- lon[seq(lon_idx$start, length.out = lon_idx$count)]
  lat_sub <- lat[seq(lat_idx$start, length.out = lat_idx$count)]
  ntime <- dim(values)[3]

  dim(values) <- c(lon_idx$count * lat_idx$count, ntime)
  lon_grid <- rep(lon_sub, times = lat_idx$count)
  lat_grid <- rep(lat_sub, each = lon_idx$count)
  mask <- build_region_mask(lon_grid, lat_grid, regions)

  nreg <- nrow(regions)
  means <- matrix(NA_real_, ntime, nreg)
  for (rr in seq_len(nreg)) {
    means[, rr] <- colMeans(values[mask[, rr], , drop = FALSE], na.rm = TRUE)
  }
  means
}

# Reads one E-OBS variable for date_sel's year and reduces it to a per-day,
# per-region mean using the same REGIONS boxes as the model output.
load_eobs_region_means <- function(file, varname, regions, date_sel) {
  ncfile <- nc_open(file)
  time_days <- ncvar_get(ncfile, "time") # days since 1950-01-01
  nc_close(ncfile)
  # POSIXct, not Date, to match the time axes used everywhere else in this
  # script (mixing Date and POSIXct on the same plot mis-scales the x-axis).
  dates <- as.POSIXct(as.Date("1950-01-01") + time_days, tz = "UTC")

  year <- format(date_sel, "%Y")
  sel <- which(format(dates, "%Y") == year)
  if (length(sel) == 0) stop("No E-OBS data for year ", year, " in: ", file)

  means <- load_gridded_region_means(file, varname, "longitude", "latitude", regions,
                                      time_start = min(sel), time_count = length(sel))
  list(dates = dates[sel], mean = means)
}

# Reads one GLEAM variable (one file per year, so no year-filtering needed
# like E-OBS) and reduces it to a per-day, per-region mean.
load_gleam_region_means <- function(file, varname, regions) {
  ncfile <- nc_open(file)
  time_days <- ncvar_get(ncfile, "time") # days since 1900-01-01
  nc_close(ncfile)
  dates <- as.POSIXct(as.Date("1900-01-01") + time_days, tz = "UTC")

  means <- load_gridded_region_means(file, varname, "lon", "lat", regions,
                                      time_start = 1, time_count = length(dates))
  list(dates = dates, mean = means)
}

# Reads one ESA CCI variable (one file per year) and reduces it to a per-day,
# per-region mean. Unlike E-OBS/GLEAM, ESA CCI's lon/lat are 2D fields on a
# curvilinear rotated-pole grid (same EUR-11 444x432 grid as ParFlow) rather
# than 1D regular axes, and the file is already cropped to the EUR-11 domain
# (no further spatial subsetting needed), so this doesn't go through
# load_gridded_region_means.
load_esacci_region_means <- function(file, varname, regions) {
  ncfile <- nc_open(file)
  lon <- as.vector(ncvar_get(ncfile, "lon"))
  lat <- as.vector(ncvar_get(ncfile, "lat"))
  time_days <- ncvar_get(ncfile, "time") # days since 1970-01-01
  values <- ncvar_get(ncfile, varname)
  nc_close(ncfile)

  dates <- as.POSIXct(as.Date("1970-01-01") + time_days, tz = "UTC")
  ntime <- dim(values)[3]
  dim(values) <- c(length(lon), ntime)
  mask <- build_region_mask(lon, lat, regions)

  nreg <- nrow(regions)
  means <- matrix(NA_real_, ntime, nreg)
  for (rr in seq_len(nreg)) {
    means[, rr] <- colMeans(values[mask[, rr], , drop = FALSE], na.rm = TRUE)
  }
  list(dates = dates, mean = means)
}

# Loads/derives the observation series for a plotvar, region-reduced like the
# model output; NULL for plotvars with no obs counterpart, or when PLOT_OBS
# is off.
load_obs_for_plotvar <- function(pv, regions, static, date_sel) {
  if (!PLOT_OBS) return(NULL)

  if (pv == "t2m") {
    tg <- load_eobs_region_means(OBS$eobs_tg_file, "tg", regions, date_sel)
    return(list(dates = tg$dates, mean = tg$mean + 273.15, source = "E-OBS")) # Celsius -> Kelvin
  }

  if (pv == "q2m") {
    tg <- load_eobs_region_means(OBS$eobs_tg_file, "tg", regions, date_sel)
    hu <- load_eobs_region_means(OBS$eobs_hu_file, "hu", regions, date_sel)
    pressure <- region_pressure_hpa(get_region_elevation(regions, static))
    q_mean <- matrix(NA_real_, nrow(tg$mean), ncol(tg$mean))
    for (rr in seq_len(ncol(tg$mean))) {
      q_mean[, rr] <- relhum_to_specific_humidity(tg$mean[, rr], hu$mean[, rr], pressure[rr])
    }
    return(list(dates = tg$dates, mean = q_mean * 1000, source = "E-OBS")) # kg/kg -> g/kg, matches plotted q2m units
  }

  if (pv == "lhf") {
    e <- load_gleam_region_means(sprintf(OBS$gleam_e_file, format(date_sel, "%Y")), "E", regions)
    return(list(dates = e$dates, mean = evap_to_latent_heat_flux(e$mean), source = "GLEAM"))
  }

  if (pv == "shf") {
    h <- load_gleam_region_means(sprintf(OBS$gleam_h_file, format(date_sel, "%Y")), "H", regions)
    return(list(dates = h$dates, mean = h$mean, source = "GLEAM")) # already W/m^2
  }

  # ESA CCI is a near-surface satellite retrieval, only comparable to the
  # shallowest soil level.
  if (is_saturation_pv(pv) && saturation_pv_level(pv) == 1) {
    sm <- load_esacci_region_means(sprintf(OBS$esacci_sm_file, format(date_sel, "%Y")), "sm", regions)
    porosity <- get_region_porosity(regions, static, saturation_pv_level(pv))
    sat_mean <- matrix(NA_real_, nrow(sm$mean), ncol(sm$mean))
    for (rr in seq_len(ncol(sm$mean))) {
      sat_mean[, rr] <- sm$mean[, rr] / porosity[rr] * 100
    }
    return(list(dates = sm$dates, mean = sat_mean, source = "ESA CCI")) # volumetric m3/m3 -> relative saturation [%]
  }

  NULL
}

# For each plotvar, load its applicable models/experiments (see select_models;
# "saturation" excludes native ICON, which has no land model), mask
# non-land/spill-over months, and reduce to per-region mean/q25/q75. Returns a
# list keyed by plotvar, each holding mean/q25/q75 arrays [time, experiment,
# region] plus the models/is_eclm vectors actually used for that plotvar.
compute_all_stats <- function(experiments, models, plotvars, regions, date_sel, static) {
  nreg <- nrow(regions)
  results <- setNames(vector("list", length(plotvars)), plotvars)

  for (pv in plotvars) {
    cat("Processing plotvar:", pv, "\n")
    pv_models <- select_models(pv, models)
    pv_experiments <- experiments[match(pv_models, models)]
    nexp <- length(pv_experiments)
    is_eclm_vec <- vapply(pv_experiments, function(e) classify_experiment(e)$is_eclm, logical(1))
    exp_stats <- vector("list", nexp)

    for (iexp in seq_along(pv_experiments)) {
      cat("  Processing experiment:", pv_experiments[iexp], "\n")
      info <- classify_experiment(pv_experiments[iexp])
      mask <- get_region_mask(info$is_icon_native, regions, static)
      fr_land <- if (info$is_icon_native) get_land_mask(static) else NULL

      vars <- get_var_names(info$is_eclm, pv)
      pattern <- build_file_pattern(pv_experiments[iexp], pv_models[iexp], info$is_eclm, date_sel)
      soil_level <- if (is_saturation_pv(pv)) saturation_pv_level(pv) else NULL
      data <- load_tsmp2_data(vars = vars, pattern = pattern, soil_level = soil_level)

      values <- if (pv == "ef") {
        data[[vars[2]]] / (data[[vars[1]]] + data[[vars[2]]])
      } else if (is_saturation_pv(pv)) {
        porosity <- get_eclm_porosity_at_level(soil_level, static, pv_experiments[iexp], pv_models[iexp], date_sel)
        data[[vars[1]]] / porosity
      } else {
        data[[vars[1]]]
      }

      time_raw <- data$time
      if (info$is_icon_native) {
        values <- values[, -1, drop = FALSE] # drop spill-over month at the start of the ICON files
        values[fr_land < 0.5, ] <- NA
        time_raw <- time_raw[-1]
      }

      exp_stats[[iexp]] <- compute_region_stats(values, mask)
      # Real hour-of-day from each experiment's own raw timestamp (fractional
      # "days since ..."), not a synthetic assumed-hourly clock -- different
      # models' output can be offset from each other in ways a shared,
      # constructed time axis won't capture (e.g. native ICON).
      exp_stats[[iexp]]$hour <- round((time_raw %% 1) * 24) %% 24
    }

    # Experiments may not all have completed the same amount of simulated
    # time yet; align on the common (shortest) leading period rather than
    # assuming every experiment produced the same number of timesteps.
    ntime <- min(vapply(exp_stats, function(s) nrow(s$mean), integer(1)))
    mean_arr <- q25_arr <- q75_arr <- array(NA_real_, c(ntime, nexp, nreg))
    hour_mat <- matrix(NA_integer_, ntime, nexp)
    for (iexp in seq_along(exp_stats)) {
      mean_arr[, iexp, ] <- exp_stats[[iexp]]$mean[seq_len(ntime), ]
      q25_arr[, iexp, ] <- exp_stats[[iexp]]$q25[seq_len(ntime), ]
      q75_arr[, iexp, ] <- exp_stats[[iexp]]$q75[seq_len(ntime), ]
      hour_mat[, iexp] <- exp_stats[[iexp]]$hour[seq_len(ntime)]
    }

    results[[pv]] <- list(mean = mean_arr, q25 = q25_arr, q75 = q75_arr,
                           models = pv_models, is_eclm = is_eclm_vec, hour = hour_mat)
  }
  results
}

# -- plot value transforms ------------------------------------------------

# filename_prefix for "saturation_levN" bakes the soil level in (e.g.
# "rsat-eurocordex_year_lev3_") so different SOIL_LEVELS don't overwrite
# each other's output files.
saturation_plot_settings <- function(plotvar) {
  list(ylab = "relative saturation [%]", ylim = c(50, 100),
       filename_prefix = paste0("rsat-eurocordex_year_lev", saturation_pv_level(plotvar), "_"))
}

get_plot_settings <- function(plotvar) {
  if (is_saturation_pv(plotvar)) return(saturation_plot_settings(plotvar))
  switch(plotvar,
    t2m = list(ylab = "2m temperature [K]", ylim = c(255, 300),
               filename_prefix = "t2m-eurocordex_year_"),
    q2m = list(ylab = "2m specific humidity [g/kg]", ylim = c(1, 20),
               filename_prefix = "q2m-eurocordex_year_"),
    lhf = list(ylab = "latent heat flux [W/m2]", ylim = c(-50, 250),
               filename_prefix = "lhf-eurocordex_year_"),
    shf = list(ylab = "sensible heat fluxes [W/m2]", ylim = c(-50, 250),
               filename_prefix = "shf-eurocordex_year_"),
    ef  = list(ylab = "Evaporative Fraction [ ]", ylim = c(0, 1),
               filename_prefix = "ef-eurocordex_year_"),
    stop("Unknown plotvar: ", plotvar))
}

get_diurnal_plot_settings <- function(plotvar) {
  if (is_saturation_pv(plotvar)) return(saturation_plot_settings(plotvar))
  switch(plotvar,
    t2m = list(ylab = "2m temperature [K]", ylim = c(270, 295),
               filename_prefix = "t2m-eurocordex_year_"),
    q2m = list(ylab = "2m specific humidity [g/kg]", ylim = c(2, 16),
               filename_prefix = "q2m-eurocordex_year_"),
    lhf = list(ylab = "latent heat flux [W/m2]", ylim = c(-50, 200),
               filename_prefix = "lhf-eurocordex_year_"),
    shf = list(ylab = "sensible heat fluxes [W/m2]", ylim = c(-50, 200),
               filename_prefix = "shf-eurocordex_year_"),
    ef  = list(ylab = "Evaporative Fraction [ ]", ylim = c(0, 1),
               filename_prefix = "ef-eurocordex_year_"),
    stop("Unknown plotvar: ", plotvar))
}

axis_ticklabels <- function(ylim, pad, by) seq(ylim[1] - pad, ylim[2] + pad, by)

rmse <- function(model_vals, obs_vals) sqrt(mean((model_vals - obs_vals)^2, na.rm = TRUE))

# is_eclm_vec: logical vector, one per experiment column, matching mean_mat's column order
transform_plot_values <- function(pv, mean_mat, q25_mat, q75_mat, is_eclm_vec) {
  if (pv == "q2m") {
    mean_mat <- mean_mat * 1000
    q25_mat <- q25_mat * 1000
    q75_mat <- q75_mat * 1000
  }
  if (pv %in% c("shf", "lhf")) {
    flip <- !is_eclm_vec
    mean_mat[, flip, ] <- -mean_mat[, flip, ]
    q25_mat[, flip, ] <- -q25_mat[, flip, ]
    q75_mat[, flip, ] <- -q75_mat[, flip, ]
  }
  if (pv == "ef") {
    mean_mat[!is.finite(mean_mat)] <- NA
    q25_mat[!is.finite(q25_mat)] <- NA
    q75_mat[!is.finite(q75_mat)] <- NA
    mean_mat <- abs(mean_mat)
    q25_mat <- abs(q25_mat)
    q75_mat <- abs(q75_mat)
  }
  if (is_saturation_pv(pv)) {
    mean_mat <- mean_mat * 100
    q25_mat <- q25_mat * 100
    q75_mat <- q75_mat * 100
  }
  list(mean = mean_mat, q25 = q25_mat, q75 = q75_mat)
}

# -- daily / diurnal aggregation ------------------------------------------

compute_daily_stats <- function(time_vector, plot_mean, pv) {
  tag <- as.Date(time_vector)
  unique_days <- as.POSIXct(unique(tag), tz = "UTC")

  if (pv == "ef") {
    # Evaporative fraction is only meaningful around local midday.
    hour <- as.integer(format(time_vector, "%H", tz = "UTC"))
    sel <- hour >= 10 & hour <= 12
    daily_mean <- apply(plot_mean, c(2, 3), function(x) tapply(x[sel], tag[sel], mean, na.rm = TRUE))
    daily_max <- apply(plot_mean, c(2, 3), function(x) tapply(x[sel], tag[sel], max, na.rm = TRUE))
    daily_min <- apply(plot_mean, c(2, 3), function(x) tapply(x[sel], tag[sel], min, na.rm = TRUE))
    daily_mean_display <- apply(daily_mean, c(2, 3), function(x) stats::filter(x, rep(1 / 7, 7), sides = 2))
  } else {
    daily_mean <- apply(plot_mean, c(2, 3), function(x) tapply(x, tag, mean, na.rm = TRUE))
    daily_max <- apply(plot_mean, c(2, 3), function(x) tapply(x, tag, max, na.rm = TRUE))
    daily_min <- apply(plot_mean, c(2, 3), function(x) tapply(x, tag, min, na.rm = TRUE))
    daily_mean_display <- daily_mean
  }
  list(days = unique_days, mean = daily_mean_display, max = daily_max, min = daily_min)
}

# hour_mat: [time, experiment] real hour-of-day per experiment (see
# compute_all_stats), not a shared/assumed clock -- each experiment's own
# raw timestamps are used, so no per-model realignment (e.g. for native
# ICON) is needed here.
# rows: optional subset of timesteps to include (e.g. a season's months);
# defaults to every timestep, i.e. the whole-year diurnal cycle.
compute_diurnal_stats <- function(hour_mat, plot_mean, plot_q25, plot_q75, rows = NULL) {
  if (is.null(rows)) rows <- seq_len(nrow(hour_mat))
  nexp <- dim(plot_mean)[2]
  nreg <- dim(plot_mean)[3]
  day_mean <- day_q25 <- day_q75 <- array(NA_real_, c(24, nexp, nreg))

  for (ii in seq_len(nexp)) {
    hour_all <- factor(hour_mat[rows, ii], levels = 0:23)
    for (rr in seq_len(nreg)) {
      day_mean[, ii, rr] <- tapply(plot_mean[rows, ii, rr], hour_all, mean, na.rm = TRUE)
      day_q75[, ii, rr] <- tapply(plot_q75[rows, ii, rr], hour_all, mean, na.rm = TRUE)
      day_q25[, ii, rr] <- tapply(plot_q25[rows, ii, rr], hour_all, mean, na.rm = TRUE)
    }
  }
  list(hour = 0:23, mean = day_mean, q25 = day_q25, q75 = day_q75)
}

# -- plotting --------------------------------------------------------------

open_device <- function(filename, plttype, width, height) {
  if (plttype == "X11") {
    X11()
  } else {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    png(filename, width = width, height = height, units = "in", res = 600)
  }
}

close_device <- function(plttype) {
  if (plttype != "X11") dev.off()
}

set_plot_par <- function() {
  par(mgp = c(3.3, 0.8, 0), mar = c(4.8, 5.5, 0.5, 0.5), cex = 1.1, cex.axis = 1.5, cex.lab = 1.5)
}

draw_shaded_lines <- function(x, mean_mat, lo_mat, hi_mat, styles, use_lty) {
  for (ii in seq_len(nrow(styles))) {
    valid <- is.finite(mean_mat[, ii]) & is.finite(lo_mat[, ii]) & is.finite(hi_mat[, ii])
    if (!any(valid)) next
    polygon(c(x[valid], rev(x[valid])), c(hi_mat[valid, ii], rev(lo_mat[valid, ii])),
            col = styles$color_alpha[ii], border = NA)
    lty <- if (use_lty) styles$linetype[ii] else 1
    lines(x[valid], mean_mat[valid, ii], lwd = 2, col = styles$color[ii], lty = lty)
  }
}

# extra_label: optional single extra legend entry (e.g. "E-OBS"), always
# solid black, appended after the model entries.
draw_legend <- function(styles, use_lty, extra_label = NULL) {
  labels <- styles$label
  colors <- styles$color
  ltys <- if (use_lty) styles$linetype else rep(1, nrow(styles))
  if (!is.null(extra_label)) {
    labels <- c(labels, extra_label)
    colors <- c(colors, "black")
    ltys <- c(ltys, 1)
  }
  legend("topleft", legend = labels, col = colors, lwd = 3, bty = "n", lty = ltys)
}

styles_for <- function(models, model_styles) {
  styles <- model_styles[match(models, model_styles$model), ]
  if (anyNA(styles$model)) stop("No style defined for model(s): ", paste(models[is.na(styles$model)], collapse = ", "))
  styles$color_alpha <- sapply(styles$color, function(col) adjustcolor(col, alpha.f = 0.15))
  styles
}

plot_hourly_series <- function(time_vector, mean_mat, q25_mat, q75_mat, styles, settings, region, dirout, plttype) {
  filename <- file.path(dirout, paste0(settings$filename_prefix, region, "_hourly.png"))
  print(paste("erstelle", filename))
  open_device(filename, plttype, width = 9.5, height = 5.9)
  on.exit(close_device(plttype))
  set_plot_par()

  ticklabels <- axis_ticklabels(settings$ylim, 5, 5)
  month_ticks <- seq(as.POSIXct(format(time_vector[1], "%Y-01-01"), tz = "UTC"), by = "month", length.out = 13)

  plot(time_vector, mean_mat[, 1], type = "n", ylim = settings$ylim, xaxt = "n",
       xlab = format(time_vector[1], "%Y"), ylab = settings$ylab, las = 1, bty = "l")
  axis(2, at = ticklabels, tcl = par("tcl") * 0.5, labels = FALSE)
  axis.POSIXct(side = 1, at = month_ticks, format = "%b")
  draw_shaded_lines(time_vector, mean_mat, q25_mat, q75_mat, styles, use_lty = FALSE)
  draw_legend(styles, use_lty = FALSE)
}

# season: NULL for the whole-year daily plot, or a season label (e.g. "DJF")
# to tag the filename/title for a season-restricted daily plot.
# obs: optional list(days, values, source) of an observation series (E-OBS
# or GLEAM, already region-reduced), drawn as a solid black line with no
# min/max band; source (e.g. "E-OBS"/"GLEAM") labels it in the legend.
plot_daily_minmaxmean <- function(daily, styles, settings, region, dirout, plttype, season = NULL, obs = NULL) {
  suffix <- if (is.null(season)) "_daily.png" else paste0("_daily_", season, ".png")
  filename <- file.path(dirout, paste0(settings$filename_prefix, region, suffix))
  if (length(daily$days) == 0) {
    message("Skipping ", filename, ": no data available yet for this period")
    return(invisible(NULL))
  }
  print(paste("erstelle", filename))
  open_device(filename, plttype, width = 9.5, height = 5.9)
  on.exit(close_device(plttype))
  set_plot_par()

  ticklabels <- axis_ticklabels(settings$ylim, 5, 5)
  month_ticks <- seq(as.POSIXct(format(daily$days[1], "%Y-01-01"), tz = "UTC"), by = "month", length.out = 13)

  plot(daily$days, daily$mean[, 1], type = "n", ylim = settings$ylim, xaxt = "n",
       xlab = format(daily$days[1], "%Y"), ylab = settings$ylab, las = 1, bty = "l")
  axis(2, at = ticklabels, tcl = par("tcl") * 0.5, labels = FALSE)
  axis.POSIXct(side = 1, at = month_ticks, format = "%b")
  if (!is.null(season)) mtext(season, side = 3, line = 0.2, cex = 0.9)
  draw_shaded_lines(daily$days, daily$mean, daily$min, daily$max, styles, use_lty = TRUE)
  if (!is.null(obs)) {
    lines(obs$days, obs$values, col = "black", lwd = 2)

    # Per-model RMSE against the obs series, over the days both have; text
    # stacked in the lower center, one line per model in its own color.
    idx <- match(daily$days, obs$days)
    paired <- !is.na(idx)
    if (any(paired)) {
      rmse_vals <- vapply(seq_len(ncol(daily$mean)), function(ii) {
        rmse(daily$mean[paired, ii], obs$values[idx[paired]])
      }, numeric(1))
      xmid <- daily$days[1] + diff(range(daily$days)) / 2
      line_height <- diff(settings$ylim) * 0.045
      y0 <- settings$ylim[1] + diff(settings$ylim) * 0.02 + (nrow(styles) - 1) * line_height
      for (ii in seq_len(nrow(styles))) {
        text(xmid, y0 - (ii - 1) * line_height,
             labels = paste0(styles$label[ii], " RMSE = ", sprintf("%.2f", rmse_vals[ii])),
             col = styles$color[ii], cex = 0.8)
      }
    }
  }
  draw_legend(styles, use_lty = TRUE, extra_label = if (!is.null(obs)) obs$source else NULL)
}

# season: NULL for the whole-year diurnal cycle, or a season label (e.g.
# "DJF") to tag the filename and title for a season-restricted diurnal cycle.
plot_diurnal_cycle <- function(diurnal, styles, settings, region, dirout, plttype, season = NULL) {
  suffix <- if (is.null(season)) "_diurnal.png" else paste0("_diurnal_", season, ".png")
  filename <- file.path(dirout, paste0(settings$filename_prefix, region, suffix))
  print(paste("erstelle", filename))
  open_device(filename, plttype, width = 6.5, height = 5.9)
  on.exit(close_device(plttype))
  set_plot_par()

  ticklabels <- axis_ticklabels(settings$ylim, 5, 5)

  plot(diurnal$hour, diurnal$mean[, 1], type = "n", ylim = settings$ylim,
       xlab = "hour of day [UTC]", ylab = settings$ylab, las = 1, bty = "l")
  axis(2, at = ticklabels, tcl = par("tcl") * 0.5, labels = FALSE)
  axis(1, at = diurnal$hour, tcl = par("tcl") * 0.5, labels = FALSE)
  if (!is.null(season)) mtext(season, side = 3, line = 0.2, cex = 0.9)
  draw_shaded_lines(diurnal$hour, diurnal$mean, diurnal$q25, diurnal$q75, styles, use_lty = TRUE)
  draw_legend(styles, use_lty = TRUE)
}

# For the ParFlow-native soil-moisture analysis: one point per day (ParFlow's
# native cadence), so a plain numeric day axis is used instead of the
# POSIXct/month-tick axis the hourly plots above use.
# time_all: POSIXct, one point per day (ParFlow's native cadence -- no
# sub-day resolution exists here, so unlike the hourly plots this has no
# "hourly" counterpart). season: NULL for the whole-year plot, or a season
# label (e.g. "DJF") to tag the filename/title for a season-restricted plot.
plot_daily_native_series <- function(time_all, mean_mat, q25_mat, q75_mat, styles, settings, region, dirout, plttype, season = NULL) {
  suffix <- if (is.null(season)) "_native.png" else paste0("_native_", season, ".png")
  filename <- file.path(dirout, paste0(settings$filename_prefix, region, suffix))
  if (length(time_all) == 0) {
    message("Skipping ", filename, ": no data available yet for this period")
    return(invisible(NULL))
  }
  print(paste("erstelle", filename))
  open_device(filename, plttype, width = 9.5, height = 5.9)
  on.exit(close_device(plttype))
  set_plot_par()

  ticklabels <- axis_ticklabels(settings$ylim, 5, 5)
  month_ticks <- seq(as.POSIXct(format(time_all[1], "%Y-01-01"), tz = "UTC"), by = "month", length.out = 13)

  plot(time_all, mean_mat[, 1], type = "n", ylim = settings$ylim, xaxt = "n",
       xlab = format(time_all[1], "%Y"), ylab = settings$ylab, las = 1, bty = "l")
  axis(2, at = ticklabels, tcl = par("tcl") * 0.5, labels = FALSE)
  axis.POSIXct(side = 1, at = month_ticks, format = "%b")
  if (!is.null(season)) mtext(season, side = 3, line = 0.2, cex = 0.9)
  draw_shaded_lines(time_all, mean_mat, q25_mat, q75_mat, styles, use_lty = TRUE)
  draw_legend(styles, use_lty = TRUE)
}

########
### ANALYSIS: regional hourly / daily / diurnal timeseries
########

run_timeseries_analysis <- function() {
  dir.create(DIROUT, recursive = TRUE, showWarnings = FALSE)

  stats <- compute_all_stats(EXPERIMENTS, MODELS, PLOTVARS, REGIONS, DATE_SEL, STATIC)

  for (pv in PLOTVARS) {
    pv_stats <- stats[[pv]]
    pv_models <- pv_stats$models
    styles <- styles_for(pv_models, MODEL_STYLES)
    settings <- get_plot_settings(pv)
    diurnal_settings <- get_diurnal_plot_settings(pv)

    ntime <- dim(pv_stats$mean)[1]
    time_vector <- seq(as.POSIXct(paste0(format(DATE_SEL, "%Y"), "-01-01 00:00:00"), tz = "UTC"),
                        by = "hour", length.out = ntime)

    plot_vals <- transform_plot_values(pv, pv_stats$mean, pv_stats$q25, pv_stats$q75, pv_stats$is_eclm)

    # Whole-year diurnal cycle, plus one per season in DIURNAL_SEASONS (each
    # restricted to that season's calendar months); computed once here
    # (covers every region) and sliced per region below.
    month_of_year <- as.integer(format(time_vector, "%m"))
    diurnal_variants <- c(
      list(list(season = NULL,
                stats = compute_diurnal_stats(pv_stats$hour, plot_vals$mean, plot_vals$q25, plot_vals$q75))),
      lapply(names(DIURNAL_SEASONS), function(s) {
        rows <- which(month_of_year %in% DIURNAL_SEASONS[[s]])
        list(season = s,
             stats = compute_diurnal_stats(pv_stats$hour, plot_vals$mean, plot_vals$q25, plot_vals$q75, rows = rows))
      })
    )

    # Same idea for the daily mean/min/max plot: whole year plus one per
    # season, restricted to that season's days.
    daily <- compute_daily_stats(time_vector, plot_vals$mean, pv)
    daily_month <- as.integer(format(daily$days, "%m"))
    daily_variants <- c(
      list(list(season = NULL, rows = seq_along(daily$days))),
      lapply(names(DIURNAL_SEASONS), function(s) {
        list(season = s, rows = which(daily_month %in% DIURNAL_SEASONS[[s]]))
      })
    )

    # E-OBS observation overlay for the daily plot only (t2m/q2m only, and
    # only if PLOT_OBS is on) -- E-OBS has no sub-daily resolution to show on
    # the hourly/diurnal plots. Subset by season independently of the
    # model's own day indices, since E-OBS's date range may not match the
    # model's (e.g. if a model run hasn't completed the full year yet).
    obs <- load_obs_for_plotvar(pv, REGIONS, STATIC, DATE_SEL)
    obs_month <- if (is.null(obs)) NULL else as.integer(format(obs$dates, "%m"))

    for (rr in seq_along(REGIONS$region)) {
      region <- REGIONS$region[rr]
      plot_hourly_series(time_vector, plot_vals$mean[, , rr], plot_vals$q25[, , rr], plot_vals$q75[, , rr],
                          styles, settings, region, DIROUT, PLTTYPE)

      for (variant in daily_variants) {
        rows <- variant$rows
        obs_slice <- NULL
        if (!is.null(obs)) {
          orows <- if (is.null(variant$season)) seq_along(obs$dates) else which(obs_month %in% DIURNAL_SEASONS[[variant$season]])
          obs_slice <- list(days = obs$dates[orows], values = obs$mean[orows, rr], source = obs$source)
        }
        plot_daily_minmaxmean(list(days = daily$days[rows], mean = daily$mean[rows, , rr], min = daily$min[rows, , rr], max = daily$max[rows, , rr]),
                               styles, settings, region, DIROUT, PLTTYPE, season = variant$season, obs = obs_slice)
      }

      for (variant in diurnal_variants) {
        d <- variant$stats
        plot_diurnal_cycle(list(hour = d$hour, mean = d$mean[, , rr], q25 = d$q25[, , rr], q75 = d$q75[, , rr]),
                            styles, diurnal_settings, region, DIROUT, PLTTYPE, season = variant$season)
      }
    }
  }
}

# -- ParFlow-specific helpers ---------------------------------------------

# Reads a ParFlow "simple ASCII" (.sa) file: header line "nx ny nz" followed
# by nx*ny*nz values in ParFlow's native x-fastest, y-next, z-slowest order,
# which is exactly R's column-major array-fill order for dim = c(nx, ny, nz).
read_parflow_sa <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  header <- as.integer(strsplit(trimws(readLines(con, n = 1)), "\\s+")[[1]])
  nx <- header[1]; ny <- header[2]; nz <- header[3]
  values <- scan(con, what = numeric(), n = nx * ny * nz, quiet = TRUE)
  array(values, dim = c(nx, ny, nz))
}

# Parses a ParFlow .tcl namelist (coup_oas.tcl style) for the van Genuchten
# retention parameters (Alpha, N, SRes, SSat) of every subsurface unit named
# via `pfset GeomInput.<unit>.Value <indicator id>`. Units without their own
# `Geom.<unit>.Saturation.*` overrides (e.g. bedrock/lake/sea units in the
# EUR-11 TSMP setup) fall back to the `Geom.domain.Saturation.*` defaults.
# Returns a data.frame with one row per unit: id, alpha, n, sres, ssat.
parse_parflow_vg_table <- function(tcl_file) {
  lines <- readLines(tcl_file)

  value_pat <- "^\\s*pfset\\s+GeomInput\\.(\\S+)\\.Value\\s+(\\S+)"
  value_lines <- grep(value_pat, lines, value = TRUE)
  unit_names <- sub(value_pat, "\\1", value_lines)
  unit_ids <- as.numeric(sub(value_pat, "\\2", value_lines))

  get_param <- function(unit, param, default) {
    pat <- paste0("^\\s*pfset\\s+Geom\\.", unit, "\\.Saturation\\.", param, "\\s+(\\S+)")
    hit <- grep(pat, lines, value = TRUE)
    if (length(hit) == 0) return(default)
    as.numeric(sub(pat, "\\1", hit[1]))
  }

  domain <- list(alpha = get_param("domain", "Alpha", NA_real_),
                  n     = get_param("domain", "N", NA_real_),
                  sres  = get_param("domain", "SRes", NA_real_),
                  ssat  = get_param("domain", "SSat", NA_real_))

  data.frame(
    id    = unit_ids,
    alpha = vapply(unit_names, get_param, numeric(1), param = "Alpha", default = domain$alpha),
    n     = vapply(unit_names, get_param, numeric(1), param = "N", default = domain$n),
    sres  = vapply(unit_names, get_param, numeric(1), param = "SRes", default = domain$sres),
    ssat  = vapply(unit_names, get_param, numeric(1), param = "SSat", default = domain$ssat),
    row.names = unit_names
  )
}

# Van Genuchten (1980) effective-saturation relation. psi is pressure head
# (ParFlow convention: psi >= 0 means saturated/ponded).
van_genuchten_saturation <- function(psi, alpha, n, sres, ssat) {
  m <- 1 - 1 / n
  se <- ifelse(psi >= 0, 1, (1 + (alpha * abs(psi))^n)^(-m))
  sres + se * (ssat - sres)
}

########
### ANALYSIS: regional relative soil saturation (ParFlow-native, daily cadence)
########

# level: which SOIL_LEVELS entry to compute this for; baked into the output
# filenames via get_plot_settings() so different levels don't collide.
run_soil_moisture_analysis <- function(level) {
  # All ParFlow output (standalone "pfl" and every "*-pfl" coupled run) only
  # exposes "pressure", never "saturation" -- confirmed via ncdump against
  # actual output files for pfl, eclm-pfl and icon-eclm-pfl. Saturation is
  # derived from pressure via van Genuchten for all of them. Uses the same
  # model list as the "saturation" plotvar (select_models() drops native
  # ICON, which has no land model).
  pv <- paste0("saturation_lev", level)
  models <- select_models(pv, MODELS)
  experiments <- file.path(BASEDIR, paste0("wfe_eur-11_revsetup_", models))
  ilev <- level
  nreg <- nrow(REGIONS)

  parflow_indicator <- read_parflow_sa(STATIC$parflow_indicator_file)
  # ParFlow's grid has the same area coverage/resolution as the eCLM grid, so
  # REGIONS masking is done directly on each model's native grid -- no
  # remapping between the two is needed.
  eclm_mask <- get_region_mask(FALSE, REGIONS, STATIC)
  parflow_mask <- get_parflow_region_mask(REGIONS, STATIC)

  mean_list <- q25_list <- q75_list <- vector("list", length(experiments))

  for (i in seq_along(experiments)) {
    cat("Processing experiment:", experiments[i], "\n")
    model_dirname <- gsub("pfl", "parflow", gsub("-", "", models[i]))
    is_parflow <- grepl("pfl", models[i])

    if (is_parflow) {
      vars <- "pressure"
      datadir <- file.path(experiments[i], "dta/simres", paste0(model_dirname, "_20*"), "out/parflow")
      pattern <- file.path(datadir, "eur-12-iic.out.*.nc")
    } else {
      vars <- "H2OSOI"
      datadir <- file.path(experiments[i], "dta/simres", paste0(model_dirname, "_20*"), "out/eclm")
      pattern <- file.path(datadir, paste0("eCLM_eur-12-iic.clm2.h0.", format(DATE_SEL, "%Y"), "*0000.nc"))
    }

    soil_level <- if (is_parflow) NULL else ilev
    data <- load_tsmp2_data(vars = vars, pattern = pattern, collapse = !is_parflow, soil_level = soil_level)
    values <- data[[vars[1]]]

    if (is_parflow) {
      orig_dim <- dim(values)
      dim(values) <- c(prod(orig_dim[1:2]), orig_dim[3], orig_dim[4])
      values[values < -1e38] <- NaN
      zlev <- orig_dim[3] - ilev + 1
      values <- values[, zlev, ]

      stopifnot(dim(parflow_indicator)[1] == orig_dim[1],
                 dim(parflow_indicator)[2] == orig_dim[2],
                 dim(parflow_indicator)[3] == orig_dim[3])

      # Each monthly restart chunk has its own copy of the namelist; the van
      # Genuchten parameters don't change between restarts, so any one of
      # them (the first) is representative of the whole experiment.
      tcl_pattern <- file.path(experiments[i], "dta/simres", paste0(model_dirname, "_20*"), "nml", "coup_oas.tcl")
      tcl_files <- Sys.glob(tcl_pattern)
      if (length(tcl_files) < 1) stop("No ParFlow namelist found matching: ", tcl_pattern)
      vg_table <- parse_parflow_vg_table(tcl_files[1])

      indicator_xy <- matrix(parflow_indicator, nrow = prod(orig_dim[1:2]), ncol = orig_dim[3])[, zlev]
      unit_idx <- match(indicator_xy, vg_table$id)
      if (anyNA(unit_idx)) stop("Indicator value(s) with no matching Geom unit in the namelist: ",
                                 paste(unique(indicator_xy[is.na(unit_idx)]), collapse = ", "))

      values <- van_genuchten_saturation(values, vg_table$alpha[unit_idx], vg_table$n[unit_idx],
                                          vg_table$sres[unit_idx], vg_table$ssat[unit_idx])
      mask <- parflow_mask
    } else {
      porosity_lev <- get_eclm_porosity_at_level(ilev, STATIC, experiments[i], models[i], DATE_SEL)
      values <- values / porosity_lev
      # eCLM output here is hourly, but ParFlow only writes once a day; average
      # eCLM down to daily so every experiment in this plot shares one cadence.
      day_tag <- floor(data$time)
      days <- sort(unique(day_tag))
      values <- sapply(days, function(d) rowMeans(values[, day_tag == d, drop = FALSE], na.rm = TRUE))
      mask <- eclm_mask
    }

    stats <- compute_region_stats(values, mask)
    mean_list[[i]] <- stats$mean
    q25_list[[i]] <- stats$q25
    q75_list[[i]] <- stats$q75
  }

  # Experiments may not all have completed the same amount of simulated time
  # yet; align on the common (shortest) leading period rather than assuming
  # every experiment produced the same number of days. Use a synthetic,
  # index-based date axis (day 1 = DATE_SEL) rather than each experiment's
  # own raw time values, consistent with the rest of the script.
  ntime <- min(vapply(mean_list, nrow, integer(1)))
  time_all <- seq(as.POSIXct(paste0(format(DATE_SEL, "%Y"), "-01-01 00:00:00"), tz = "UTC"),
                   by = "day", length.out = ntime)

  nexp <- length(experiments)
  mean_arr <- q25_arr <- q75_arr <- array(NA_real_, c(ntime, nexp, nreg))
  for (i in seq_along(experiments)) {
    mean_arr[, i, ] <- mean_list[[i]][seq_len(ntime), ] * 100
    q25_arr[, i, ] <- q25_list[[i]][seq_len(ntime), ] * 100
    q75_arr[, i, ] <- q75_list[[i]][seq_len(ntime), ] * 100
  }

  styles <- styles_for(models, MODEL_STYLES)
  settings <- get_plot_settings(pv)

  # Whole-year plot, plus one per season in DIURNAL_SEASONS (restricted to
  # that season's calendar months). ParFlow's native cadence is daily, so
  # unlike the hourly plotvars there's no sub-day "hourly"/diurnal
  # counterpart for this analysis.
  month_of_year <- as.integer(format(time_all, "%m"))
  season_variants <- c(
    list(list(season = NULL, rows = seq_len(ntime))),
    lapply(names(DIURNAL_SEASONS), function(s) {
      list(season = s, rows = which(month_of_year %in% DIURNAL_SEASONS[[s]]))
    })
  )

  for (rr in seq_along(REGIONS$region)) {
    region <- REGIONS$region[rr]
    for (variant in season_variants) {
      rows <- variant$rows
      plot_daily_native_series(time_all[rows], mean_arr[rows, , rr], q25_arr[rows, , rr], q75_arr[rows, , rr],
                                styles, settings, region, DIROUT, PLTTYPE, season = variant$season)
    }
  }
}

########
### MAIN
########

# "saturation" in PLOTVARS expands to one "saturation_levN" entry per
# SOIL_LEVELS; the rest of the script treats each as its own plotvar (see
# is_saturation_pv()/saturation_pv_level()).
PLOTVARS <- unlist(lapply(PLOTVARS, function(pv) {
  if (pv == "saturation") paste0("saturation_lev", SOIL_LEVELS) else pv
}))

main <- function() {
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  run_timeseries_analysis()
  if (any(is_saturation_pv(PLOTVARS))) {
    for (level in SOIL_LEVELS) run_soil_moisture_analysis(level)
  }
}

main()
