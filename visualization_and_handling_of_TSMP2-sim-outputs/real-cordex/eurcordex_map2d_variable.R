#!/usr/bin/env Rscript
#
# General-purpose 2D map plotter for a single EURO-CORDEX / TSMP2 output
# variable (ICON, eCLM, or ParFlow), styled exactly like the 2D maps in
# eurcordex_soilmoisture_eclm_vs_parflow.R:
#   - eCLM/ICON: one flat-shaded triangle per cell on the unstructured mesh
#     (rotated-pole view), drawn with base graphics polygon().
#   - ParFlow: image() on its native (rotated-pole-regular) raster grid.
# Both use the same viridis colour scale and fields::image.plot() legend.
#
# Unlike that script, this one is not a coupling-consistency comparison: it
# reads and maps exactly one variable from exactly one component/experiment,
# selected on the command line, so it can be pointed at any variable without
# editing the file.
#
# Usage:
#   Rscript eurcordex_map2d_variable.R <component> <variable> [options]
#
#   component   "eclm" | "icon" | "parflow"
#   variable    NetCDF variable name to read (e.g. H2OSOI, TSA, t_2m, pres,
#               pressure)
#
# Options:
#   --experiment=NAME    experiment name as in wfe_eur-11_revsetup_NAME
#                         (default: icon-eclm-pfl)
#   --date=YYYY-MM-DD     only the year is used, to build the file-search
#                         pattern / to build ParFlow's synthetic time axis
#                         (default: 2018-01-01)
#   --level=N             1-based index into the variable's vertical/soil
#                         dimension (levsoi, levgrnd, height, height_2, plev,
#                         ParFlow z, ...). Required if the variable has such
#                         a dimension; ignored (with a warning) otherwise.
#   --reduce=mean|min|max time reduction applied per cell to build the map
#                         (default: mean, i.e. an annual/period mean like the
#                         source script)
#   --season=DJF|MAM|JJA|SON  restrict the time reduction to that season's
#                         calendar months (default: the whole available
#                         period)
#   --zlim=lo,hi          fixed colour-scale range (default: 1st/99th
#                         percentile of the reduced field, as in the source
#                         script)
#   --landmask             ParFlow only: mask out ocean/lake cells via the
#                         static grid file's LANDMASK, as the source script
#                         does for soil saturation (default: off, since not
#                         every ParFlow variable is land-only)
#   --outdir=PATH          output directory (default: see DIROUT below)
#
# Examples:
#   Rscript eurcordex_map2d_variable.R eclm H2OSOI --level=1
#   Rscript eurcordex_map2d_variable.R icon t_2m --season=JJA
#   Rscript eurcordex_map2d_variable.R parflow pressure --level=1 --landmask

suppressPackageStartupMessages({
  library(ncdf4)
  library(abind)
  library(fields) # image.plot() colour-bar legend
})

########
### CONFIG
########

BASEDIR  <- "/p/scratch/cslts/poll1/sim/paper"
DIROUT   <- "/p/project1/cslts/poll1/tools/plots_tsmp2-map2d-variable/"
PLTTYPE  <- "png" # "png" or "X11"

STATIC <- list(
  eclm_grid_file    = "/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/mkmapgrids/EUR-R13B05_189976_grid.nc",
  icon_grid_file    = "/p/scratch/cslts/poll1/sim/paper/wfe_eur-11_icon-eclm-pfl/dta/geo/icon/static/europe011_DOM01.nc",
  parflow_grid_file = "/p/scratch/cslts/poll1/git/TSMP_EUR-11/static/clm/griddata_CLM_EUR-11_TSMP_FZJ-IBG3_CLMPFLDomain_444x432.nc"
)

# Per-component file layout. horiz_dim is the NetCDF dimension name that
# indexes gridcells, used to identify a variable's "extra" (vertical) axis by
# elimination -- everything that isn't "time" or horiz_dim. ParFlow isn't
# listed here: it has two horizontal dims (x, y) instead of one, so it's
# handled by its own loader below rather than this generic one.
COMPONENT_INFO <- list(
  eclm = list(subdir = "out/eclm", file_glob = "eCLM_eur-12-iic.clm2.h0.%s*0000.nc",
              horiz_dim = "lndgrid", grid_file = STATIC$eclm_grid_file),
  icon = list(subdir = "out/icon", file_glob = "ICON_out_EU-R13B5_inst_DOM01_ML_%s*T000000Z_1h.nc",
              horiz_dim = "ncells", grid_file = STATIC$icon_grid_file)
)
COMPONENT_LABELS <- list(eclm = "eCLM", icon = "ICON", parflow = "ParFlow")

SEASONS <- list(DJF = c(12, 1, 2), MAM = c(3, 4, 5), JJA = c(6, 7, 8), SON = c(9, 10, 11))

RAD2DEG <- 180 / pi

# Rotated-pole parameters of the EUR-11 CORDEX domain, matching the companion
# Python triangulation tool (map_2d_trigrid_vis_modularized_LAM.py) and the
# eCLM-vs-ParFlow script, so this map shows the same rotated-pole view.
LON_NORTHPOLE <- -162.0
LAT_NORTHPOLE <- 39.25

PALETTE_COLORS <- hcl.colors(100, "Viridis")

########
### ARGUMENT PARSING
########

USAGE <- paste(
  "Usage: Rscript eurcordex_map2d_variable.R <component: eclm|icon|parflow> <variable> [options]",
  "Options: --experiment=NAME --date=YYYY-MM-DD --level=N --reduce=mean|min|max",
  "         --season=DJF|MAM|JJA|SON --zlim=lo,hi --landmask --outdir=PATH",
  sep = "\n"
)

parse_args <- function(argv) {
  if (length(argv) < 2) stop(USAGE, call. = FALSE)
  component <- argv[1]
  variable <- argv[2]
  if (!component %in% names(COMPONENT_LABELS)) {
    stop("component must be one of: ", paste(names(COMPONENT_LABELS), collapse = ", "), " (got '", component, "')\n", USAGE, call. = FALSE)
  }

  opts <- list(experiment = "icon-eclm-pfl", date = "2018-01-01", level = NULL,
               reduce = "mean", season = NULL, zlim = NULL, landmask = FALSE, outdir = DIROUT)
  for (arg in argv[-(1:2)]) {
    if (arg == "--landmask") { opts$landmask <- TRUE; next }
    m <- regmatches(arg, regexec("^--([a-zA-Z]+)=(.*)$", arg))[[1]]
    if (length(m) != 3) stop("Unrecognized argument: ", arg, "\n", USAGE, call. = FALSE)
    key <- m[2]; val <- m[3]
    if (key == "experiment") opts$experiment <- val
    else if (key == "date") opts$date <- val
    else if (key == "level") opts$level <- as.integer(val)
    else if (key == "reduce") opts$reduce <- val
    else if (key == "season") opts$season <- val
    else if (key == "zlim") opts$zlim <- as.numeric(strsplit(val, ",")[[1]])
    else if (key == "outdir") opts$outdir <- val
    else stop("Unrecognized option: --", key, "\n", USAGE, call. = FALSE)
  }

  if (!opts$reduce %in% c("mean", "min", "max")) stop("--reduce must be one of: mean, min, max (got '", opts$reduce, "')", call. = FALSE)
  if (!is.null(opts$season) && !opts$season %in% names(SEASONS)) stop("--season must be one of: ", paste(names(SEASONS), collapse = ", "), call. = FALSE)
  if (!is.null(opts$zlim) && length(opts$zlim) != 2) stop("--zlim needs exactly two comma-separated values, e.g. --zlim=0,100", call. = FALSE)

  list(component = component, variable = variable, experiment = opts$experiment,
       date_sel = as.Date(opts$date), level = opts$level, reduce = opts$reduce,
       season = opts$season, zlim = opts$zlim, landmask = opts$landmask, outdir = opts$outdir)
}

########
### DATA LOADING
########

# eCLM's own "time" declares its units as "days since <this monthly chunk's
# own start date>", not one fixed reference for the whole run -- read as a
# common absolute day count (days since 1970-01-01) instead, so concatenating
# time across chunks doesn't alias different months onto the same value.
# Native ICON's "time" isn't "days since ..." either -- it's an absolute
# encoded date/time (units "day as %Y%m%d.%f", e.g. 20180105.0417), decoded
# separately in as_dates() below, so it's passed through unchanged here.
read_time_raw <- function(ncfile, component) {
  raw <- ncvar_get(ncfile, "time")
  if (component == "icon") return(raw)
  units_attr <- ncatt_get(ncfile, "time", "units")
  if (!units_attr$hasatt || !grepl("^days since", units_attr$value)) return(raw)
  origin <- as.POSIXct(sub("^days since\\s*", "", units_attr$value), tz = "UTC")
  raw + as.numeric(difftime(origin, as.POSIXct("1970-01-01", tz = "UTC"), units = "days"))
}

as_dates <- function(component, time_raw, date_sel) {
  if (component == "icon") return(as.Date(sprintf("%08d", floor(time_raw)), format = "%Y%m%d"))
  as.Date("1970-01-01") + floor(time_raw) # eclm
}

# Reads one variable for eCLM or native ICON (single-horizontal-dim
# components) across every matching monthly/period file, concatenated along
# time. If the variable has one extra (non-time, non-horizontal) dimension --
# levsoi/levgrnd for eCLM, height/height_2/plev/... for ICON -- `level`
# selects a single index via NetCDF start/count instead of loading every
# level into memory before subsetting.
load_unstructured_variable <- function(component, variable, experiment_dir, simres_dirname, date_sel, level) {
  info <- COMPONENT_INFO[[component]]
  pattern <- file.path(experiment_dir, "dta/simres", paste0(simres_dirname, "_20*"), info$subdir,
                        sprintf(info$file_glob, format(date_sel, "%Y")))
  files <- Sys.glob(pattern)
  if (length(files) == 0) stop("No matching files found for pattern: ", pattern)

  # eCLM-only: each monthly restart chunk's own out/eclm directory contains
  # one extra, spurious single-timestep snapshot for day 1 of the *following*
  # month, alongside the real, full file of the same name inside that next
  # chunk -- keep only the larger (complete) file wherever a basename repeats
  # (see eurcordex_soilmoisture_eclm_vs_parflow.R's load_tsmp2_data for the
  # full story).
  if (component == "eclm") {
    dupe_names <- unique(basename(files)[duplicated(basename(files))])
    for (bn in dupe_names) {
      idx <- which(basename(files) == bn)
      files <- files[-idx[-which.max(file.size(files[idx]))]]
    }
  }

  value_list <- vector("list", length(files))
  time_list <- vector("list", length(files))
  units <- NULL

  for (ii in seq_along(files)) {
    ncfile <- nc_open(files[ii])
    cat("Reading:", files[ii], "\n")
    if (!variable %in% names(ncfile$var)) stop("Variable '", variable, "' not found in: ", files[ii])
    v <- ncfile$var[[variable]]
    if (is.null(units)) {
      units_attr <- ncatt_get(ncfile, variable, "units")
      if (units_attr$hasatt) units <- units_attr$value
    }

    # v$dim lists dimensions in ncdump's declared order, but ncvar_get()'s
    # start/count expect the reverse of that (the returned array's own,
    # fastest-varying-first order) -- e.g. dim 2 of 3 declared dims is dim
    # (3-2+1)=2 in that reversed order.
    dim_names <- vapply(v$dim, function(d) d$name, character(1))
    extra <- which(!dim_names %in% c("time", info$horiz_dim))
    if (length(extra) > 1) {
      stop("Variable '", variable, "' has more than one non-horizontal, non-time dimension (",
           paste(dim_names[extra], collapse = ", "), "); this script only supports variables with at most one extra dimension.")
    }
    if (length(extra) == 1) {
      if (is.null(level)) {
        stop("Variable '", variable, "' has dimension '", dim_names[extra], "' (length ", v$dim[[extra]]$len,
             ") -- pass --level=<1..", v$dim[[extra]]$len, "> to select one.")
      }
      lev_dim <- v$ndims - extra + 1
      start <- rep(1, v$ndims); start[lev_dim] <- level
      count <- rep(-1, v$ndims); count[lev_dim] <- 1
      value_list[[ii]] <- ncvar_get(ncfile, variable, start = start, count = count, collapse = TRUE)
    } else {
      if (!is.null(level)) warning("--level was given but variable '", variable, "' has no extra dimension; ignoring.")
      value_list[[ii]] <- ncvar_get(ncfile, variable, collapse = TRUE)
    }
    time_list[[ii]] <- read_time_raw(ncfile, component)
    nc_close(ncfile)
  }

  along_dim <- length(dim(value_list[[length(value_list)]]))
  values <- do.call(abind, c(value_list, list(along = along_dim))) # (horiz, time)
  list(values = values, time_raw = unlist(time_list), units = units)
}

# Reads one variable from ParFlow's own output (one file per day, from the
# *same* experiment) for the ParFlow-native raster grid. ParFlow files have
# exactly one timestep each (a length-1 time dim that ncvar_get() would
# otherwise silently drop), so collapse must stay FALSE, unlike the
# eCLM/ICON loader above -- see eurcordex_soilmoisture_eclm_vs_parflow.R's
# load_tsmp2_data for the full story. `level` (if the variable has a z
# dimension) is applied by direct indexing after loading rather than via
# start/count, since ParFlow's two horizontal dims (x, y) don't fit the
# single-extra-dimension scheme used above.
load_parflow_variable <- function(variable, experiment_dir, simres_dirname, level) {
  datadir <- file.path(experiment_dir, "dta/simres", paste0(simres_dirname, "_20*"), "out/parflow")
  pattern <- file.path(datadir, "eur-12-iic.out.*.nc")
  files <- Sys.glob(pattern)
  if (length(files) == 0) stop("No matching files found for pattern: ", pattern)

  value_list <- vector("list", length(files))
  units <- NULL
  for (ii in seq_along(files)) {
    ncfile <- nc_open(files[ii])
    cat("Reading:", files[ii], "\n")
    if (!variable %in% names(ncfile$var)) stop("Variable '", variable, "' not found in: ", files[ii])
    if (is.null(units)) {
      units_attr <- ncatt_get(ncfile, variable, "units")
      if (units_attr$hasatt) units <- units_attr$value
    }
    value_list[[ii]] <- ncvar_get(ncfile, variable, collapse = FALSE)
    nc_close(ncfile)
  }

  orig_dim <- dim(value_list[[1]]) # (x, y[, z], time)
  values <- do.call(abind, c(value_list, list(along = length(orig_dim))))
  values[values < -1e38] <- NaN # ParFlow's fill-value convention

  has_z <- length(orig_dim) == 4
  if (has_z) {
    nx <- orig_dim[1]; ny <- orig_dim[2]; nz <- orig_dim[3]
    dim(values) <- c(nx * ny, nz, dim(values)[length(dim(values))])
    if (is.null(level)) stop("Variable '", variable, "' has a z dimension (length ", nz, ") -- pass --level=<1..", nz, "> to select one.")
    zlev <- nz - level + 1 # ParFlow's z increases upward; level 1 (shallowest) is the top layer
    values <- values[, zlev, ]
  } else {
    nx <- orig_dim[1]; ny <- orig_dim[2]
    dim(values) <- c(nx * ny, dim(values)[length(dim(values))])
    if (!is.null(level)) warning("--level was given but variable '", variable, "' has no z dimension; ignoring.")
  }

  list(values = values, nx = nx, ny = ny, units = units)
}

########
### GRIDS
########

# eCLM's and native ICON's grid files are both in the same "ICON grid"
# format, and both carry clon_vertices/clat_vertices(cell, nv) -- each cell's
# 3 vertex lon/lat directly, in radians -- so this one loader (and the
# tripcolor-style plot below) works for either component without needing
# ICON's separate vertex_of_cell shared-vertex table.
load_triangulation <- function(grid_file) {
  ncfile <- nc_open(grid_file)
  vlon <- t(ncvar_get(ncfile, "clon_vertices")) * RAD2DEG # ncvar_get returns (nv, ncells); transpose to (ncells, nv)
  vlat <- t(ncvar_get(ncfile, "clat_vertices")) * RAD2DEG
  nc_close(ncfile)
  list(vlon = vlon, vlat = vlat)
}

# LONGXY/LATIXY are declared (lsmlat, lsmlon); ncvar_get() reverses this to R
# dim (lsmlon, lsmlat) = (x, y) -- the same orientation as ParFlow's own
# output fields, so no reordering is needed before reshaping values back onto
# this grid.
load_parflow_grid <- function(grid_file) {
  ncfile <- nc_open(grid_file)
  lon <- ncvar_get(ncfile, "LONGXY") # (nx, ny)
  lat <- ncvar_get(ncfile, "LATIXY")
  landmask <- ncvar_get(ncfile, "LANDMASK")
  nc_close(ncfile)
  list(lon = lon, lat = lat, landmask = landmask, nx = dim(lon)[1], ny = dim(lon)[2])
}

# Manual rotated-pole transform, matching the companion Python tool's
# rotate_coordinates() and eurcordex_soilmoisture_eclm_vs_parflow.R exactly.
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

########
### PLOTTING
########

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

# Draws one flat-shaded triangle per eCLM/ICON cell, coloured by `values` --
# the same visual result as the companion Python tool's
# ax.tripcolor(shading='flat'), using base graphics polygon() (which accepts
# one colour per NA-separated sub-polygon) instead of matplotlib.
plot_map_tri <- function(values, vlon_rot, vlat_rot, zlim, palette_colors, title, legend_lab, filename, plttype) {
  # Row i's 3 vertices followed by NA, flattened row-major so each
  # NA-separated chunk in x/y is one triangle in cell order.
  x <- as.vector(t(cbind(vlon_rot, NA_real_)))
  y <- as.vector(t(cbind(vlat_rot, NA_real_)))

  ncolors <- length(palette_colors)
  breaks <- seq(zlim[1], zlim[2], length.out = ncolors + 1)
  color_idx <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  cell_colors <- palette_colors[color_idx] # NA values (masked/missing cells) stay unfilled

  open_device(filename, plttype, width = 6.5, height = 6.5)
  on.exit(close_device(plttype))
  # Right margin (mar[4]) must leave enough room for image.plot()'s colour
  # bar, its tick labels *and* legend.lab, or the rotated label text is
  # pushed past the device edge and silently clipped -- 7 lines is enough for
  # legend.mar = 5 plus the label itself at the default text size.
  par(mar = c(1, 1, 3, 7))

  plot(NA, xlim = range(vlon_rot), ylim = range(vlat_rot), asp = 1,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", main = title, cex.main = 1.1)
  polygon(x, y, col = cell_colors, border = NA)
  image.plot(zlim = zlim, col = palette_colors, legend.only = TRUE, add = TRUE,
             legend.mar = 5, legend.lab = legend_lab)
}

# ParFlow's native grid is regular once rotated into the same rotated-pole
# frame the eCLM/ICON triangles are drawn in, so unlike those meshes this map
# can just use image() directly, no per-cell polygon needed.
plot_map_grid <- function(values_mat, rlon_axis, rlat_axis, zlim, palette_colors, title, legend_lab, filename, plttype) {
  open_device(filename, plttype, width = 6.5, height = 6.5)
  on.exit(close_device(plttype))
  par(mar = c(1, 1, 3, 7))

  image(rlon_axis, rlat_axis, values_mat, zlim = zlim, col = palette_colors, asp = 1,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = title, cex.main = 1.1)
  image.plot(zlim = zlim, col = palette_colors, legend.only = TRUE, add = TRUE,
             legend.mar = 5, legend.lab = legend_lab)
}

########
### HELPERS
########

season_columns <- function(dates, season) {
  if (is.null(season)) return(seq_along(dates))
  which(as.integer(format(dates, "%m")) %in% SEASONS[[season]])
}

# rowMeans() is used for the common "mean" case since it's much faster than
# apply() over what can be hundreds of thousands of rows; min/max fall back
# to apply() since base R has no built-in rowMins/rowMaxs.
reduce_rows <- function(m, reduce) {
  if (reduce == "mean") return(rowMeans(m, na.rm = TRUE))
  fun <- if (reduce == "min") min else max
  apply(m, 1, function(x) {
    v <- suppressWarnings(fun(x, na.rm = TRUE))
    if (is.infinite(v)) NA_real_ else v
  })
}

build_filename <- function(args) {
  parts <- c(args$component, args$variable)
  if (!is.null(args$level)) parts <- c(parts, paste0("lev", args$level))
  parts <- c(parts, args$reduce)
  if (!is.null(args$season)) parts <- c(parts, args$season)
  file.path(args$outdir, paste0("map2d_", paste(parts, collapse = "_"), ".png"))
}

build_title <- function(component_label, args) {
  reduce_label <- c(mean = "mean", min = "minimum", max = "maximum")[[args$reduce]]
  level_txt <- if (!is.null(args$level)) paste0(" (level ", args$level, ")") else ""
  season_txt <- if (!is.null(args$season)) paste0(" ", args$season) else ""
  paste0(component_label, " ", args$variable, level_txt, " -- ", reduce_label, " ", format(args$date_sel, "%Y"), season_txt)
}

legend_label <- function(units, variable) {
  if (!is.null(units) && nzchar(units)) paste0(variable, " [", units, "]") else variable
}

########
### MAIN
########

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  simres_dirname <- gsub("pfl", "parflow", gsub("-", "", args$experiment))
  experiment_dir <- file.path(BASEDIR, paste0("wfe_eur-11_revsetup_", args$experiment))

  if (args$component == "parflow") {
    d <- load_parflow_variable(args$variable, experiment_dir, simres_dirname, args$level)
    # ParFlow's "time" is a bare per-chunk step index with no calendar
    # meaning; Sys.glob's lexical file order is chronological (month-chunk
    # dirs and same-width step suffixes both sort correctly), so a synthetic
    # daily axis starting at --date works, as in the source script.
    dates <- args$date_sel + (seq_len(ncol(d$values)) - 1)
    cols <- season_columns(dates, args$season)
    cell_stat <- reduce_rows(d$values[, cols, drop = FALSE], args$reduce)

    grid <- load_parflow_grid(STATIC$parflow_grid_file)
    if (args$landmask) cell_stat[as.vector(grid$landmask) == 0] <- NA
    zlim <- if (!is.null(args$zlim)) args$zlim else quantile(cell_stat, c(0.01, 0.99), na.rm = TRUE)

    rot <- rotate_to_pole(grid$lon, grid$lat, LON_NORTHPOLE, LAT_NORTHPOLE)
    values_mat <- matrix(cell_stat, nrow = grid$nx, ncol = grid$ny)

    filename <- build_filename(args)
    plot_map_grid(values_mat, rot$lon[, 1], rot$lat[1, ], zlim, PALETTE_COLORS,
                  title = build_title("ParFlow", args), legend_lab = legend_label(d$units, args$variable),
                  filename = filename, plttype = PLTTYPE)
  } else {
    info <- COMPONENT_INFO[[args$component]]
    d <- load_unstructured_variable(args$component, args$variable, experiment_dir, simres_dirname, args$date_sel, args$level)
    dates <- as_dates(args$component, d$time_raw, args$date_sel)
    cols <- season_columns(dates, args$season)
    cell_stat <- reduce_rows(d$values[, cols, drop = FALSE], args$reduce)

    tri <- load_triangulation(info$grid_file)
    zlim <- if (!is.null(args$zlim)) args$zlim else quantile(cell_stat, c(0.01, 0.99), na.rm = TRUE)

    rot <- rotate_to_pole(tri$vlon, tri$vlat, LON_NORTHPOLE, LAT_NORTHPOLE)
    vlon_rot <- matrix(rot$lon, ncol = 3)
    vlat_rot <- matrix(rot$lat, ncol = 3)

    filename <- build_filename(args)
    plot_map_tri(cell_stat, vlon_rot, vlat_rot, zlim, PALETTE_COLORS,
                title = build_title(COMPONENT_LABELS[[args$component]], args), legend_lab = legend_label(d$units, args$variable),
                filename = filename, plttype = PLTTYPE)
  }

  cat("Wrote:", filename, "\n")
}

main()
