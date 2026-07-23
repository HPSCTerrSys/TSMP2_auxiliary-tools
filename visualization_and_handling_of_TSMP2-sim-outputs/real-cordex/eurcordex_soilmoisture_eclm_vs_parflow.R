#!/usr/bin/env Rscript
#
# eCLM vs. ParFlow soil moisture coupling-consistency check (EURO-CORDEX /
# TSMP2, single "eclm-pfl" experiment).
#
# Focused variant of eurcordex_timeline_regions_year_imprv.R: within the one
# ParFlow-coupled experiment ("eclm-pfl"), compares the two land-surface
# components' own view of soil saturation --
#   - eCLM's diagnosed H2OSOI / porosity, on eCLM's unstructured grid
#   - ParFlow's own pressure state, converted to saturation via the van
#     Genuchten relation (using the subsurface indicator file and that
#     experiment's own coup_oas.tcl parameters), on ParFlow's native
#     444x432 raster grid
# This is a coupling-consistency check (do the two model components agree on
# the soil-moisture state at their shared interface?), not a comparison of
# two different simulations.
#
# For each soil level in SOIL_LEVELS, this produces:
#   - regional (mean/25th/75th percentile) daily timeseries per REGIONS, one
#     component against the other, with an optional ESA CCI overlay for the
#     shallowest level. Daily is ParFlow's native cadence (it writes once a
#     day); eCLM's hourly output is averaged down to match, as in the source
#     script's run_soil_moisture_analysis() -- so unlike the regular
#     regional-timeseries analysis, there is no hourly/diurnal-cycle plot
#     here.
#   - two 2D maps of annual-mean saturation (one per component, same colour
#     scale so they're visually comparable) -- eCLM's on its triangulated
#     mesh (styled after map_2d_trigrid_vis_modularized_LAM.py), ParFlow's on
#     its native raster grid. The two live on genuinely different
#     discretizations, so no cell-by-cell difference map is produced; only a
#     shared-scale side-by-side comparison.
#
# Usage:
#   Rscript eurcordex_soilmoisture_eclm_vs_parflow.R
#
# All configuration (experiment, date, paths, plot styling) is in the CONFIG
# section below.

suppressPackageStartupMessages({
  library(ncdf4)
  library(dplyr)
  library(abind)
  library(tibble)
  library(fields) # image.plot() colour-bar legend for the 2D maps only
})

########
### CONFIG
########

BASEDIR   <- "/p/scratch/cslts/poll1/sim/paper"
DIROUT    <- "/p/project1/cslts/poll1/tools/plots_tsmp2-soilmoisture_eclm-vs-parflow/"
PLTTYPE   <- "png"                        # "png" or "X11"
DATE_SEL  <- as.Date("2018-01-01")        # only year is used to build file-search patterns; simulation start

# Soil layers (indices into eCLM's levsoi / ParFlow's z, top = 1) to compare
# -- one full set of outputs (regional timeseries + 2D maps) is produced per
# level, with the level baked into every output filename so levels don't
# overwrite each other.
SOIL_LEVELS <- c(1, 3)

# The single experiment both components are read from.
EXPERIMENT_MODEL_NAME <- "eclm-pfl"
EXPERIMENT_DIR <- file.path(BASEDIR, paste0("wfe_eur-11_revsetup_", EXPERIMENT_MODEL_NAME))
SIMRES_DIRNAME <- gsub("pfl", "parflow", gsub("-", "", EXPERIMENT_MODEL_NAME)) # "eclmparflow"

STATIC <- list(
  eclm_grid_file = "/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/mkmapgrids/EUR-R13B05_189976_grid.nc",
  eclm_surf_file = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/eclm/static/surfdata_ICON-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230302_gcvurb-pfsoil_halo.nc",
  parflow_indicator_file = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/parflow/static/EUR-11_TSMP_FZJ-IBG3_eCLMPFLDomain_444x432_INDICATOR_regridded_rescaled_SoilGrids250-v2017_BGR3_alv.sa",
  parflow_grid_file = "/p/scratch/cslts/poll1/git/TSMP_EUR-11/static/clm/griddata_CLM_EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432.nc"
)

# ESA CCI gridded soil moisture (EUR-11 rotated-pole grid, one file per year),
# overlaid on the daily regional plot for the shallowest soil level only (a
# near-surface satellite retrieval, not comparable to deeper levels). Set
# PLOT_OBS to FALSE to disable.
PLOT_OBS <- TRUE
OBS_ESACCI_SM_FILE <- "/p/data1/jibg31/DETECT/observations/ESACCI_SM/eu_cordex11/daily/merged/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-%s-fv08.1_EUCORDEX011.nc"

# Seasons (by calendar month) to additionally break the daily regional plot
# out by, on top of the whole-year plot. Set to list() to disable and only
# produce the whole-year plot.
SEASONS <- list(
  DJF = c(12, 1, 2),
  MAM = c(3, 4, 5),
  JJA = c(6, 7, 8),
  SON = c(9, 10, 11)
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

COMPONENT_STYLES <- tribble(
  ~component, ~label,     ~color,        ~linetype,
  "eclm",     "eCLM",     "forestgreen", 3,
  "parflow",  "ParFlow",  "sienna",      1
)

# Rotated-pole parameters of the EUR-11 CORDEX domain, matching the companion
# Python triangulation tool (map_2d_trigrid_vis_modularized_LAM.py), so the 2D
# maps below show the same rotated-pole view.
LON_NORTHPOLE <- -162.0
LAT_NORTHPOLE <- 39.25

########
### FUNCTIONS
########

# -- data loading -------------------------------------------------------

# collapse: pass FALSE for ParFlow output, which has exactly one timestep per
# file (a length-1 time dim) that ncvar_get() would otherwise silently drop,
# breaking the "combine along the last (time) dimension" abind below. Must
# stay TRUE (default) for eCLM, whose H2OSOI has a genuine singleton
# dimension (the selected soil level) that needs dropping.
#
# soil_level: for H2OSOI's levsoi dimension, read only this one level via
# NetCDF start/count instead of loading every level into memory before
# subsetting -- H2OSOI has 20 levels, so loading them all for a full year of
# hourly output is ~20x more data than needed and will exhaust memory.
#
# dedupe_spillover: eCLM-only (see below); pass FALSE for ParFlow output,
# whose file names ("...out.00000.nc", ...00001.nc, ...) are a per-chunk
# step index that restarts at 0 every month -- every month's directory reuses
# the same basenames for genuinely different days, so basename-deduping them
# would wrongly discard eleven of every twelve months down to just one.
#
# eCLM's own "time" declares its units as "days since <this monthly chunk's
# own start date>", not one fixed reference for the whole run (confirmed via
# ncdump: a Jan-05 file and a Jun-05 file both report time = 4, each
# relative to its own month's 1st) -- read as a common absolute day count
# (days since 1970-01-01) instead, so concatenating time across chunks
# doesn't alias different months onto the same value (e.g. day_tag grouping
# in load_saturation_eclm_component below). ParFlow's "time" has no CF units
# attribute at all (a bare step index); passed through unchanged since it's
# never used to derive a calendar date.
read_time_days <- function(ncfile) {
  raw <- ncvar_get(ncfile, "time")
  units_attr <- ncatt_get(ncfile, "time", "units")
  if (!units_attr$hasatt || !grepl("^days since", units_attr$value)) return(raw)
  origin <- as.POSIXct(sub("^days since\\s*", "", units_attr$value), tz = "UTC")
  raw + as.numeric(difftime(origin, as.POSIXct("1970-01-01", tz = "UTC"), units = "days"))
}

load_tsmp2_data <- function(vars, pattern, collapse = TRUE, soil_level = NULL, dedupe_spillover = TRUE) {
  files <- Sys.glob(pattern)
  if (length(files) == 0) stop("No matching files found for pattern: ", pattern, " var(s): ", paste(vars, collapse = ", "))

  # Each monthly eCLM restart chunk's own out/eclm directory contains one
  # extra, spurious single-timestep snapshot for day 1 of the *following*
  # month (e.g. a ~60MB "...2018-02-01-00000.nc" inside the January chunk,
  # alongside the real, full 24-timestep ~1.25GB file of the same name
  # inside the February chunk) -- glob patterns spanning multiple months
  # pick up both, and abind() below fails on the shape mismatch. Keep only
  # the larger (complete) file wherever a basename is matched more than once.
  if (dedupe_spillover) {
    dupe_names <- unique(basename(files)[duplicated(basename(files))])
    for (bn in dupe_names) {
      idx <- which(basename(files) == bn)
      files <- files[-idx[-which.max(file.size(files[idx]))]]
    }
  }

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

# eCLM writes its own exact saturated water content (WATSAT, on levgrnd) into
# the first history file of each monthly restart chunk -- prefer that over
# the PCT_SAND/ORGANIC approximation in soil_porosity() above, since it's
# the literal value CLM used internally, not an estimate of it. levsoi (used
# by H2OSOI/"level" everywhere else in this script) is always the top subset
# of levgrnd, so WATSAT[, level] lines up directly. Falls back to
# soil_porosity() if no available file has WATSAT.
get_eclm_porosity_at_level <- function(level, static) {
  pattern <- file.path(EXPERIMENT_DIR, "dta/simres", paste0(SIMRES_DIRNAME, "_20*"), "out/eclm",
                        paste0("eCLM_eur-12-iic.clm2.h0.", format(DATE_SEL, "%Y"), "*0000.nc"))
  for (f in sort(Sys.glob(pattern))) {
    ncfile <- nc_open(f)
    has_watsat <- "WATSAT" %in% names(ncfile$var)
    if (has_watsat) watsat <- ncvar_get(ncfile, "WATSAT")[, level]
    nc_close(ncfile)
    if (has_watsat) return(watsat)
  }
  soil_porosity(static$eclm_surf_file)[, level]
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

# -- regional masks -------------------------------------------------------

build_region_mask <- function(clon_deg, clat_deg, regions) {
  sapply(seq_len(nrow(regions)), function(rr) {
    clon_deg >= regions$xmin[rr] & clon_deg <= regions$xmax[rr] &
      clat_deg >= regions$ymin[rr] & clat_deg <= regions$ymax[rr]
  })
}

get_eclm_region_mask <- function(regions, static) {
  ncfile <- nc_open(static$eclm_grid_file)
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

# Mean eCLM porosity (at the given soil level) per REGIONS box, used to
# convert ESA CCI's raw volumetric soil moisture into the same
# relative-saturation [%] units the model output is plotted in.
get_region_porosity <- function(regions, static, level) {
  ncfile <- nc_open(static$eclm_grid_file)
  clon <- ncvar_get(ncfile, "clon") * RAD2DEG
  clat <- ncvar_get(ncfile, "clat") * RAD2DEG
  nc_close(ncfile)
  porosity <- get_eclm_porosity_at_level(level, static)
  mask <- build_region_mask(clon, clat, regions)
  vapply(seq_len(nrow(regions)), function(rr) mean(porosity[mask[, rr]], na.rm = TRUE), numeric(1))
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

# -- observations (ESA CCI) ------------------------------------------------

# Reads one ESA CCI variable (one file per year) and reduces it to a per-day,
# per-region mean. ESA CCI's lon/lat are 2D fields on a curvilinear
# rotated-pole grid (same EUR-11 domain as the model output) rather than 1D
# regular axes, and the file is already cropped to the EUR-11 domain (no
# further spatial subsetting needed).
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

# Loads the ESA CCI observation series for the shallowest soil level,
# region-reduced like the model output; NULL if PLOT_OBS is off or the level
# isn't the shallowest one.
load_obs_for_level <- function(level, regions, static, date_sel) {
  if (!PLOT_OBS || level != 1) return(NULL)
  sm <- load_esacci_region_means(sprintf(OBS_ESACCI_SM_FILE, format(date_sel, "%Y")), "sm", regions)
  porosity <- get_region_porosity(regions, static, level)
  sat_mean <- matrix(NA_real_, nrow(sm$mean), ncol(sm$mean))
  for (rr in seq_len(ncol(sm$mean))) {
    sat_mean[, rr] <- sm$mean[, rr] / porosity[rr] * 100
  }
  list(dates = sm$dates, mean = sat_mean, source = "ESA CCI") # volumetric m3/m3 -> relative saturation [%]
}

# -- per-component loading (region stats + per-cell annual mean) ----------

# filename_prefix bakes the soil level in (e.g. "rsat-eclm_vs_parflow_lev3_")
# so different SOIL_LEVELS don't overwrite each other's output files.
plot_settings_for_level <- function(level) {
  list(ylab = "relative saturation [%]", ylim = c(50, 100),
       filename_prefix = paste0("rsat-eclm_vs_parflow_lev", level, "_"))
}

# Reads eCLM's own H2OSOI for one soil level (hourly cadence) and converts it
# to relative saturation via porosity. Returns region-reduced statistics
# already averaged down to ParFlow's native daily cadence, and the per-cell
# annual mean (for the 2D map).
load_saturation_eclm_component <- function(level, region_mask, static) {
  porosity_lev <- get_eclm_porosity_at_level(level, static)
  datadir <- file.path(EXPERIMENT_DIR, "dta/simres", paste0(SIMRES_DIRNAME, "_20*"), "out/eclm")
  pattern <- file.path(datadir, paste0("eCLM_eur-12-iic.clm2.h0.", format(DATE_SEL, "%Y"), "*0000.nc"))

  data <- load_tsmp2_data(vars = "H2OSOI", pattern = pattern, soil_level = level)
  values <- data$H2OSOI / porosity_lev

  # eCLM output here is hourly, but ParFlow only writes once a day; average
  # eCLM down to daily so both components share one cadence for comparison.
  day_tag <- floor(data$time)
  days <- sort(unique(day_tag))
  values <- sapply(days, function(d) rowMeans(values[, day_tag == d, drop = FALSE], na.rm = TRUE))

  list(region_stats = compute_region_stats(values, region_mask), cell_mean = rowMeans(values, na.rm = TRUE) * 100)
}

# Reads ParFlow's own pressure field for one soil level (one file per day)
# from the *same* eclm-pfl experiment, converts it to relative saturation via
# the van Genuchten relation using the indicator file and this experiment's
# own coup_oas.tcl parameters. Returns region-reduced daily statistics and
# the per-cell annual mean (for the 2D map).
load_saturation_parflow_component <- function(level, region_mask, static) {
  datadir <- file.path(EXPERIMENT_DIR, "dta/simres", paste0(SIMRES_DIRNAME, "_20*"), "out/parflow")
  pattern <- file.path(datadir, "eur-12-iic.out.*.nc")
  data <- load_tsmp2_data(vars = "pressure", pattern = pattern, collapse = FALSE, dedupe_spillover = FALSE)

  orig_dim <- dim(data$pressure) # (x, y, z, time)
  nx <- orig_dim[1]; ny <- orig_dim[2]; nz <- orig_dim[3]
  values <- data$pressure
  dim(values) <- c(nx * ny, nz, orig_dim[4])
  values[values < -1e38] <- NaN

  zlev <- nz - level + 1 # ParFlow's z increases upward; level 1 (shallowest) is the top layer
  values <- values[, zlev, ]

  parflow_indicator <- read_parflow_sa(static$parflow_indicator_file)
  stopifnot(dim(parflow_indicator)[1] == nx, dim(parflow_indicator)[2] == ny, dim(parflow_indicator)[3] == nz)

  # Each monthly restart chunk has its own copy of the namelist; the van
  # Genuchten parameters don't change between restarts, so any one of them
  # (the first) is representative of the whole experiment.
  tcl_pattern <- file.path(EXPERIMENT_DIR, "dta/simres", paste0(SIMRES_DIRNAME, "_20*"), "nml", "coup_oas.tcl")
  tcl_files <- Sys.glob(tcl_pattern)
  if (length(tcl_files) < 1) stop("No ParFlow namelist found matching: ", tcl_pattern)
  vg_table <- parse_parflow_vg_table(tcl_files[1])

  indicator_xy <- matrix(parflow_indicator, nrow = nx * ny, ncol = nz)[, zlev]
  unit_idx <- match(indicator_xy, vg_table$id)
  if (anyNA(unit_idx)) stop("Indicator value(s) with no matching Geom unit in the namelist: ",
                             paste(unique(indicator_xy[is.na(unit_idx)]), collapse = ", "))

  values <- van_genuchten_saturation(values, vg_table$alpha[unit_idx], vg_table$n[unit_idx],
                                      vg_table$sres[unit_idx], vg_table$ssat[unit_idx])

  list(region_stats = compute_region_stats(values, region_mask), cell_mean = rowMeans(values, na.rm = TRUE) * 100)
}

# -- regional plotting ------------------------------------------------------

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

# extra_label: optional single extra legend entry (e.g. "ESA CCI"), always
# solid black, appended after the component entries.
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

styles_for <- function(components, component_styles) {
  styles <- component_styles[match(components, component_styles$component), ]
  if (anyNA(styles$component)) stop("No style defined for component(s): ", paste(components[is.na(styles$component)], collapse = ", "))
  styles$color_alpha <- sapply(styles$color, function(col) adjustcolor(col, alpha.f = 0.15))
  styles
}

axis_ticklabels <- function(ylim, pad, by) seq(ylim[1] - pad, ylim[2] + pad, by)

rmse <- function(component_vals, obs_vals) sqrt(mean((component_vals - obs_vals)^2, na.rm = TRUE))

# ParFlow writes once a day; eCLM's hourly output is downsampled to match
# this cadence (see load_saturation_eclm_component), so both components
# share one native cadence here -- unlike the full regional-timeseries
# script, there is no meaningful hourly or diurnal-cycle comparison for a
# ParFlow variable.
# time_all: a synthetic daily POSIXct axis (day 1 = DATE_SEL), not either
# component's own raw timestamps -- consistent with how the source script's
# ParFlow-native analysis builds its time axis.
# season: NULL for the whole-year plot, or a season label (e.g. "DJF") to tag
# the filename/title for a season-restricted plot.
# obs: optional list(days, values, source) of an observation series (ESA
# CCI), drawn as a solid black line with a per-component RMSE readout.
plot_daily_native_series <- function(time_all, mean_mat, q25_mat, q75_mat, styles, settings, region, dirout, plttype, season = NULL, obs = NULL) {
  suffix <- if (is.null(season)) "_daily.png" else paste0("_daily_", season, ".png")
  filename <- file.path(dirout, paste0(settings$filename_prefix, region, suffix))
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
  if (!is.null(obs)) {
    lines(obs$days, obs$values, col = "black", lwd = 2)

    # Per-component RMSE against the obs series, over the days both have;
    # matched by calendar day since time_all is a synthetic axis while obs$days
    # are ESA CCI's real dates.
    idx <- match(as.Date(time_all), as.Date(obs$days))
    paired <- !is.na(idx)
    if (any(paired)) {
      rmse_vals <- vapply(seq_len(ncol(mean_mat)), function(ii) {
        rmse(mean_mat[paired, ii], obs$values[idx[paired]])
      }, numeric(1))
      xmid <- time_all[1] + diff(range(time_all)) / 2
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

# -- 2D maps -----------------------------------------------------------------
# eCLM's map is styled after map_2d_trigrid_vis_modularized_LAM.py's
# tripcolor plot: coastlines are skipped since soil moisture is only defined
# over land there, so the shape of the drawn (non-NA) cells already reads as
# the coastline. ParFlow's map uses image() instead of per-cell polygons
# (see plot_soil_moisture_map_grid for why that's valid here), with ocean
# cells masked out via LANDMASK for the same reason.

# eCLM's ICON-format grid file: clon_vertices/clat_vertices give each cell's
# 3 vertex lon/lat directly (unlike ICON's native vertex_of_cell scheme,
# which needs a separate shared-vertex table), so one flat-shaded triangle
# per cell can be drawn straight from these without any extra indexing.
load_eclm_triangulation <- function(grid_file) {
  ncfile <- nc_open(grid_file)
  vlon <- t(ncvar_get(ncfile, "clon_vertices")) * RAD2DEG # ncvar_get returns (nv, ncells); transpose to (ncells, nv)
  vlat <- t(ncvar_get(ncfile, "clat_vertices")) * RAD2DEG
  nc_close(ncfile)
  list(vlon = vlon, vlat = vlat)
}

# LONGXY/LATIXY are declared (lsmlat, lsmlon); ncvar_get() reverses this to R
# dim (lsmlon, lsmlat) = (x, y) -- the same orientation as ParFlow's own
# output fields, so no reordering is needed before reshaping the saturation
# values back onto this grid for the 2D map.
load_parflow_grid <- function(grid_file) {
  ncfile <- nc_open(grid_file)
  lon <- ncvar_get(ncfile, "LONGXY") # (nx, ny)
  lat <- ncvar_get(ncfile, "LATIXY")
  landmask <- ncvar_get(ncfile, "LANDMASK")
  nc_close(ncfile)
  list(lon = lon, lat = lat, landmask = landmask, nx = dim(lon)[1], ny = dim(lon)[2])
}

# Manual rotated-pole transform, matching the companion Python tool's
# rotate_coordinates() exactly, so the maps here share its rotated-pole view
# without depending on a GIS/cartopy equivalent in R.
rotate_to_pole <- function(lon, lat, lon_northpole, lat_northpole) {
  lon0 <- lon_northpole * pi / 180
  lat0 <- lat_northpole * pi / 180
  lon_r <- lon * pi / 180
  lat_r <- lat * pi / 180

  lon_rot <- atan2(
    -cos(lat_r) * sin(lon_r - lon0),
    -cos(lat_r) * sin(lat0) * cos(lon_r - lon0) + sin(lat_r) * cos(lat0)
  ) * 180 / pi
  lon_rot[lon_rot < -180] <- lon_rot[lon_rot < -180] + 360
  lon_rot[lon_rot > 180] <- lon_rot[lon_rot > 180] - 360

  lat_rot <- asin(
    sin(lat_r) * sin(lat0) + cos(lat_r) * cos(lat0) * cos(lon_r - lon0)
  ) * 180 / pi

  list(lon = lon_rot, lat = lat_rot)
}

# Draws one flat-shaded triangle per eCLM cell, coloured by `values` -- the
# same visual result as the Python tool's ax.tripcolor(shading='flat'), using
# base graphics polygon() (which accepts one colour per NA-separated
# sub-polygon) instead of matplotlib.
plot_soil_moisture_map_tri <- function(values, vlon_rot, vlat_rot, zlim, palette_colors, title, legend_lab, filename, plttype) {
  # Row i's 3 vertices followed by NA, flattened row-major so each
  # NA-separated chunk in x/y is one triangle in cell order.
  x <- as.vector(t(cbind(vlon_rot, NA_real_)))
  y <- as.vector(t(cbind(vlat_rot, NA_real_)))

  ncolors <- length(palette_colors)
  breaks <- seq(zlim[1], zlim[2], length.out = ncolors + 1)
  color_idx <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  cell_colors <- palette_colors[color_idx] # NA values (masked/ocean cells) stay unfilled

  open_device(filename, plttype, width = 6.5, height = 6.5)
  on.exit(close_device(plttype))
  # Right margin (mar[4]) must leave enough room for image.plot()'s colour
  # bar, its tick labels *and* legend.lab, or the rotated label text is
  # pushed past the device edge and silently clipped -- 7 lines is enough
  # for legend.mar = 5 plus the label itself at the default text size.
  par(mar = c(1, 1, 3, 7))

  plot(NA, xlim = range(vlon_rot), ylim = range(vlat_rot), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", main = title, cex.main = 1.1)
  polygon(x, y, col = cell_colors, border = NA)
  image.plot(zlim = zlim, col = palette_colors, legend.only = TRUE, add = TRUE,
             legend.mar = 5, legend.lab = legend_lab)
}

# ParFlow's native grid is regular once rotated into the same rotated-pole
# frame the eCLM triangles are drawn in -- LONGXY/LATIXY are, by
# construction, the geographic lon/lat of a plain regular rotated-pole
# raster (verified: rotating them gives a grid regular to ~1e-6 degrees) --
# so unlike the eCLM mesh, this map can just use image() directly, no
# per-cell polygon needed.
plot_soil_moisture_map_grid <- function(values_mat, rlon_axis, rlat_axis, zlim, palette_colors, title, legend_lab, filename, plttype) {
  open_device(filename, plttype, width = 6.5, height = 6.5)
  on.exit(close_device(plttype))
  par(mar = c(1, 1, 3, 7))

  image(rlon_axis, rlat_axis, values_mat, zlim = zlim, col = palette_colors, asp = 1,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = title, cex.main = 1.1)
  image.plot(zlim = zlim, col = palette_colors, legend.only = TRUE, add = TRUE,
             legend.mar = 5, legend.lab = legend_lab)
}

# Builds and saves the two 2D maps (per-component annual mean, shared colour
# scale) for one soil level, from the already loaded per-cell time-means --
# no data is re-read from disk here. No difference map: the two components
# live on different grids (eCLM's unstructured mesh vs. ParFlow's raster),
# so a cell-by-cell difference would need regridding first.
run_map_analysis <- function(cell_mean_eclm, cell_mean_parflow, level, tri, pfl_grid, styles, dirout, plttype) {
  rot_eclm <- rotate_to_pole(tri$vlon, tri$vlat, LON_NORTHPOLE, LAT_NORTHPOLE)
  vlon_rot <- matrix(rot_eclm$lon, ncol = 3)
  vlat_rot <- matrix(rot_eclm$lat, ncol = 3)

  rot_pfl <- rotate_to_pole(pfl_grid$lon, pfl_grid$lat, LON_NORTHPOLE, LAT_NORTHPOLE)
  rlon_axis <- rot_pfl$lon[, 1]
  rlat_axis <- rot_pfl$lat[1, ]
  parflow_vals_mat <- matrix(cell_mean_parflow, nrow = pfl_grid$nx, ncol = pfl_grid$ny)
  # ParFlow simulates ocean/lake cells too (with bedrock-like van Genuchten
  # units), but those aren't comparable to eCLM's land-only saturation, so
  # they're masked out of this visualisation (regional stats are left as-is,
  # matching the source script's ParFlow-native analysis).
  parflow_vals_mat[pfl_grid$landmask == 0] <- NA

  ncolors <- 100
  seq_colors <- hcl.colors(ncolors, "Viridis")
  zlim_abs <- quantile(c(cell_mean_eclm, parflow_vals_mat[pfl_grid$landmask == 1]),
                        c(0.01, 0.99), na.rm = TRUE)

  eclm_label <- styles$label[styles$component == "eclm"]
  filename_eclm <- file.path(dirout, paste0("map_soilmoisture_eclm_lev", level, ".png"))
  print(paste("erstelle", filename_eclm))
  plot_soil_moisture_map_tri(cell_mean_eclm, vlon_rot, vlat_rot, zlim_abs, seq_colors,
                              title = paste0(eclm_label, " mean relative saturation, soil level ", level),
                              legend_lab = "relative saturation [%]", filename = filename_eclm, plttype = plttype)

  pfl_label <- styles$label[styles$component == "parflow"]
  filename_pfl <- file.path(dirout, paste0("map_soilmoisture_parflow_lev", level, ".png"))
  print(paste("erstelle", filename_pfl))
  plot_soil_moisture_map_grid(parflow_vals_mat, rlon_axis, rlat_axis, zlim_abs, seq_colors,
                               title = paste0(pfl_label, " mean relative saturation, soil level ", level),
                               legend_lab = "relative saturation [%]", filename = filename_pfl, plttype = plttype)
}

########
### ANALYSIS
########

run_level_analysis <- function(level, eclm_mask, parflow_mask, tri, pfl_grid) {
  cat("Processing soil level:", level, "\n")

  eclm <- load_saturation_eclm_component(level, eclm_mask, STATIC)
  pfl <- load_saturation_parflow_component(level, parflow_mask, STATIC)

  styles <- styles_for(c("eclm", "parflow"), COMPONENT_STYLES)
  settings <- plot_settings_for_level(level)

  # Experiments may not all have completed the same amount of simulated time
  # yet, and eCLM/ParFlow could in principle differ by a day at the very end
  # of a run; align on the common (shortest) leading period. Uses a
  # synthetic, index-based date axis (day 1 = DATE_SEL) rather than either
  # component's own raw time values, matching the source script's
  # ParFlow-native analysis.
  ntime <- min(nrow(eclm$region_stats$mean), nrow(pfl$region_stats$mean))
  nreg <- ncol(eclm_mask)
  mean_arr <- q25_arr <- q75_arr <- array(NA_real_, c(ntime, 2, nreg))
  mean_arr[, 1, ] <- eclm$region_stats$mean[seq_len(ntime), ] * 100
  q25_arr[, 1, ] <- eclm$region_stats$q25[seq_len(ntime), ] * 100
  q75_arr[, 1, ] <- eclm$region_stats$q75[seq_len(ntime), ] * 100
  mean_arr[, 2, ] <- pfl$region_stats$mean[seq_len(ntime), ] * 100
  q25_arr[, 2, ] <- pfl$region_stats$q25[seq_len(ntime), ] * 100
  q75_arr[, 2, ] <- pfl$region_stats$q75[seq_len(ntime), ] * 100

  time_all <- seq(as.POSIXct(paste0(format(DATE_SEL, "%Y"), "-01-01 00:00:00"), tz = "UTC"),
                   by = "day", length.out = ntime)

  # Whole-year plot, plus one per season in SEASONS (restricted to that
  # season's calendar days).
  month_of_year <- as.integer(format(time_all, "%m"))
  season_variants <- c(
    list(list(season = NULL, rows = seq_len(ntime))),
    lapply(names(SEASONS), function(s) {
      list(season = s, rows = which(month_of_year %in% SEASONS[[s]]))
    })
  )

  obs <- load_obs_for_level(level, REGIONS, STATIC, DATE_SEL)

  for (rr in seq_along(REGIONS$region)) {
    region <- REGIONS$region[rr]
    for (variant in season_variants) {
      rows <- variant$rows
      if (length(rows) == 0) next # season not reached yet by the currently-available output
      obs_slice <- NULL
      if (!is.null(obs)) {
        orows <- if (is.null(variant$season)) seq_along(obs$dates) else which(as.integer(format(obs$dates, "%m")) %in% SEASONS[[variant$season]])
        obs_slice <- list(days = obs$dates[orows], values = obs$mean[orows, rr], source = obs$source)
      }
      plot_daily_native_series(time_all[rows], mean_arr[rows, , rr], q25_arr[rows, , rr], q75_arr[rows, , rr],
                                styles, settings, region, DIROUT, PLTTYPE, season = variant$season, obs = obs_slice)
    }
  }

  run_map_analysis(eclm$cell_mean, pfl$cell_mean, level, tri, pfl_grid, styles, DIROUT, PLTTYPE)
}

########
### MAIN
########

main <- function() {
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  dir.create(DIROUT, recursive = TRUE, showWarnings = FALSE)

  eclm_mask <- get_eclm_region_mask(REGIONS, STATIC)
  parflow_mask <- get_parflow_region_mask(REGIONS, STATIC)
  tri <- load_eclm_triangulation(STATIC$eclm_grid_file)
  pfl_grid <- load_parflow_grid(STATIC$parflow_grid_file)

  for (level in SOIL_LEVELS) run_level_analysis(level, eclm_mask, parflow_mask, tri, pfl_grid)
}

main()
