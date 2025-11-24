#!/usr/bin/env python
# coding: utf-8

'''
Visualisation of global ParFlow model input and output on a map.

Purpose:
This script combines the plotting scripts for the different figures in 
Kollet et al. "Global groundwater modeling: Proof-of-concept of 3D 
variably saturated flow simulation at km-resolution" in JHX for further
reuse and adjustment. 

Inputs:
- Locally stored compilation of global datasets.
- ParFlow simulation results and diagnostics. Differ in netCDF structure and attributes.
- ParFlow external parameter files.
- Additional data from external sources: glacier mask, permafrost mask.
- All data is stored here: `/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources`

Ouputs:
- Image files
- Examples: `/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources/example_plots`

Usage:
`python3 map_2D_ParFlow-global-highres.py`

Notes:
Plot size, data set size in grid elements, map scale, dpi, image file format are
adjusted with each other. Also the viewer makes a difference.
The following combinations are implemented:
- Fig1 hydrofacies + global
- Fig2 Reff + global
- Fig3 Sr + global
- Fig4a WTD + global
- Fig4b WTD + southAmerica1
- Fig poster WTD + southAmerica2 -- some poster vis, less map content
- Fig EoCoE WTD + europe -- some website vis, less map content
Total runtime about 4min.

Limitations:
- May need manual code edits, especially in `main`. 
- Shape-file is not yet part of the data object.
- "fullres" size, works, but not yet implemented as an option.
- Adjustable resampling not yet implemented.

Operating environment:
- locally with conda forge env (tested), or
- JSC/JURECA JSC/JUWELS frontnodes (tested), source this environment initialisation file:
  https://icg4geo.icg.kfa-juelich.de/Organisation/hpc_it_dev_data_howtos-and-bestpractice/hpc_scientific_operating_environments/-/blob/main/loadenv.JURECA-DC_JUWELS_2023_GCC-OpenMPI_Std-Py-AI.ini`
'''

import time, sys, os
#os.environ["CARTOPY_DATA_DIR"] = "/p/home/jusers/goergen1/jureca/.local/share/cartopy"

import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FixedLocator
from matplotlib import colormaps
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
import matplotlib.cm as cm

#import cartopy
#cartopy.config['data_dir']='/p/home/jusers/goergen1/jureca/.local/share/cartopy'
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
import cartopy.io.img_tiles as cimgt

from dataclasses import dataclass, field


__authors__ = "Klaus GOERGEN"
__orcid__ = "https://orcid.org/0000-0002-4208-3444"
__email__ = ""
__maintainer__ = "Klaus GOERGEN"
__copyright__ = "Copyright (c) 2025, https://www.fz-juelich.de (FZJ)"
__license__ = "MIT"
__version__ = "1.0.0"
__date__ = "2025-11-16"
__status__ = "Production"
__credits__ = [
    "Stefan KOLLET",
    "Bibi NAZ",
    "Alexandre BELLEFLAMME"
    ]
__acknowledgements__ = "" 
__references__ = ""


@dataclass
class simAuxData:
    """
    Container class for simulation datasets.

    This class loads multiple datasets from disk (NumPy arrays or xarray datasets)
    and provides basic access to them for further analysis or visualization.

    Parameters
    ----------
    *_dir : str
        Path to the first dataset file.
    *_file : str
        Filename of dataset file.
    *_varname : str
        Name of netCDF variable to import. 
    *_resample : int, {0, 1}
        Flag indicating whether the imported data array resolution shall be reduced.
        0 disables resampling, data stays as is; 1 enables it, reduction factor of 4 is hardcoded.                       
        Not all datasets need this.

    Attributes
    ----------
    *_data : np.ndarray
        Loaded data.
    *_lon, lon : np.ndarray
        Loaded longitude data; decimal degrees; center point.
        Not all datasets provide usable lon/lat information, then it is not read.
    *_lat, lat : np.ndarray
        Loaded latitude data; decimal degrees; center point.
        Not all datasets provide usable lon/lat information, then it is not read.
    """

    sim_var1_dir: str
    sim_var1_file: str
    sim_var1_varname: str
    sim_var1_resample: int
    sim_var1_data: np.ndarray = field(init=False)
    sim_var2_dir: str
    sim_var2_file: str
    sim_var2_varname: str
    sim_var2_resample: int
    sim_var2_data: np.ndarray = field(init=False)
    mask_landsea_dir: str
    mask_landsea_file: str
    mask_landsea_varname: str
    mask_landsea_resample: int
    mask_landsea_data: np.ndarray = field(init=False)
    mask_permafrost_dir: str
    mask_permafrost_file: str
    mask_permafrost_varname: str
    mask_permafrost_resample: int
    mask_permafrost_data: np.ndarray = field(init=False)
    mask_permafrost_lon: np.ndarray = field(init=False)
    mask_permafrost_lat: np.ndarray = field(init=False)
    mask_karst_dir: str
    mask_karst_file: str
    mask_karst_varname: str
    mask_karst_resample: int
    mask_karst_data: np.ndarray = field(init=False)
    sim_lon: np.ndarray = field(init=False)
    sim_lat: np.ndarray = field(init=False)
    forcing_var1_dir: str
    forcing_var1_file: str
    forcing_var1_varname: str
    forcing_var1_resample: int
    forcing_var1_data: np.ndarray = field(init=False)
    extpar_var1_dir: str
    extpar_var1_file: str
    extpar_var1_varname: str
    extpar_var1_resample: int
    extpar_var1_data: np.ndarray = field(init=False)

    def __post_init__(self):
        self.sim_var1_data = self._read_sim_var1(self.sim_var1_dir, self.sim_var1_file, self.sim_var1_varname, self.sim_var1_resample)
        self.sim_var2_data = self._read_sim_var2(self.sim_var2_dir, self.sim_var2_file, self.sim_var2_varname, self.sim_var2_resample)
        self.mask_landsea_data = self._read_mask_landsea(self.mask_landsea_dir, self.mask_landsea_file, self.mask_landsea_varname, self.mask_landsea_resample)
        self.mask_permafrost_data, self.mask_permafrost_lon, self.mask_permafrost_lat = self._read_mask_permafrost(self.mask_permafrost_dir, self.mask_permafrost_file, self.mask_permafrost_varname, self.mask_permafrost_resample)
        self.mask_karst_data, self.sim_lon, self.sim_lat = self._read_mask_karst(self.mask_karst_dir, self.mask_karst_file, self.mask_karst_varname, self.mask_karst_resample)
        self.forcing_var1_data = self._read_forcing_var1(self.forcing_var1_dir, self.forcing_var1_file, self.forcing_var1_varname, self.forcing_var1_resample)
        self.extpar_var1_data = self._read_extpar_var1(self.extpar_var1_dir, self.extpar_var1_file, self.extpar_var1_varname, self.extpar_var1_resample)
  
    def _read_sim_var1(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        ds = xr.open_dataset(pn+"/"+fn)    
        if regr == 0:
            #data = ds[varname].values
            data = ds[varname].isel(z=0, time=0).values
        if regr == 1:
            # resample, factor 4 along x and y axes
            # align the number of datapoints with the nnumber of pixels to print
            # data contains NaN
            data = ds[varname].isel(z=0, time=0)
            datac = data.coarsen(x=4, y=4, boundary='exact', side='left').mean()
            del data
            data = datac.values

        #print('min', np.nanmin(data))
        #print('max', np.nanmax(data))

        return data

    def _read_sim_var2(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        ds = xr.open_dataset(pn+"/"+fn)    
        if regr == 0:
            data = ds[varname].values
        if regr == 1:
            # resample, factor 4 along x and y axes
            # align the number of datapoints with the nnumber of pixels to print
            # data contains NaN
            data = ds[varname]
            datac = data.coarsen(nx=4, ny=4, boundary='exact', side='left').mean()
            del data
            data = datac.values

        return data

    def _read_mask_landsea(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        # this is a normal 0 1 mask, nothing with masked array and nothing with 0-1, but 0,1
        ds_mask = xr.open_dataset(pn+"/"+fn)
        if regr == 0:
            mask = ds_mask[varname].values
            #print('mask lsm: ', mask.shape)
            #print(isinstance(mask, np.ma.MaskedArray))
            #print(mask)
            #print(mask.min())
            #print(mask.max())
            #print(np.unique(mask))
        if regr == 1:
            # resample, factor 4 along x and y axes
            # align the number of datapoints with the nnumber of pixels to print
            mask = ds_mask[varname]
            maskc = mask.coarsen(x=4, y=4, boundary='exact', side='left').mean()
            #print('mask: ', mask.shape)
            #print('maskc: ', maskc.shape)
            del mask
            mask = maskc.values
            # averaging a binary mask where, >=0.75, i.e. 3 of 4 grid elements == 1
            mask[mask < 0.75] = 0

        return mask

    def _read_mask_permafrost(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        #print('----------------------------------------')
        ds_mask = xr.open_dataset(pn+"/"+fn, decode_times=False)
        mask = ds_mask[varname].isel(time=0) #.values
        lon =  ds_mask['lon'] #.values
        lat =  ds_mask['lat'] #.values
        #print('mask isimip lon:', lon.shape)
        #print('mask isimip lat:', lat.shape)

        #if mtype == 'permafrost':
        #    print('permafrost')
        #    mask = np.ma.masked_where(mask == 4, mask)

        if regr == 0:
            mask = mask.values
            lon = lon.values
            lat = lat.values

        # xarray, considers NaNs, but mean is from valid data only
        if regr == 1:
            maskc = mask.coarsen(lat=4, lon=4, boundary='exact', side='left').mean()
            lonc = lon.coarsen(lon=4, boundary='exact', side='left').mean()
            latc = lat.coarsen(lat=4, boundary='exact', side='left').mean()
            #print('mask: ', mask.shape)
            #print('maskc: ', maskc.shape)
            #print('lon: ', lon.shape)
            #print('lonc: ', lonc.shape)
            #print('lat: ', lat.shape)
            #print('latc: ', latc.shape)
            del mask
            del lon
            del lat
            mask = maskc.values
            lon = lonc.values
            lat = latc.values

        mask = np.ma.masked_where(mask <= 25, mask)

        return mask, lon, lat

    def _read_mask_karst(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        #print('----------------------------------------')
        ds_mask = xr.open_dataset(pn+"/"+fn, decode_times=False)
        mask = ds_mask[varname].isel(time=0) #.values
        lon =  ds_mask['lon']
        lat =  ds_mask['lat']

        #print('karst')
        # the NaNs are ignored
        # 0.06 porosity matches with https://doi.org/10.1002/2014GL059856 Gleeson et al. (2014) "carbonate" class
        mask[:] = np.where(mask == 0.06, 1, 0)
        #mask[:] = np.where((mask == 0.06) | (mask == 0.12), 1, 0)

        if regr == 0:
            mask = mask.values
            lon = lon.values
            lat = lat.values

        if regr == 1:
            maskc = mask.coarsen(lat=4, lon=4, boundary='exact', side='left').mean()
            lonc = lon.coarsen(lon=4, boundary='exact', side='left').mean()
            latc = lat.coarsen(lat=4, boundary='exact', side='left').mean()
            #print('mask: ', mask.shape)
            #print('maskc: ', maskc.shape)
            #print('lon: ', lon.shape)
            #print('lonc: ', lonc.shape)
            #print('lat: ', lat.shape)
            #print('latc: ', latc.shape)
            del mask
            del lon
            del lat
            mask = maskc.values
            lon = lonc.values
            lat = latc.values
            mask[mask < 0.01] = 0   # 0.75

        return mask, lon, lat

    def _read_forcing_var1(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        ds = xr.open_dataset(pn+"/"+fn, decode_times=False)
        if regr == 0:
            data = ds[varname].values
            print('data basic:', data.shape)
        if regr == 1:
            # resample, factor 4 along x and y axes
            # align the number of datapoints with the nnumber of pixels to print
            # data contains no NaN
            # min -1e-10
            # max 0.028552477762256483
            # must be m / day
            data = ds[varname].isel(nz=14, time=0)
            datac = data.coarsen(nx=4, ny=4, boundary='exact', side='left').mean()
            #print('data: ', data.shape)
            #print('datac: ', datac.shape)
            del data
            data = datac.values

            #print('min', np.nanmin(data))
            #print('max', np.nanmax(data))
            #print('mean', np.nanmean(data))

        # h-1 -> m/h  -> m/year ->  mm / year
        return data *  0.02 * 365. * 24.  * 1000.

    def _read_extpar_var1(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        ds = xr.open_dataset(pn+"/"+fn)    
        if regr == 0:
            #data = ds[varname].values
            data = ds[varname].isel(z=0, time=0).values
        if regr == 1:
            # resample, factor 4 along x and y axes
            # align the number of datapoints with the nnumber of pixels to print
            # data contains NaN
            data = ds[varname].isel(z=0, time=0)
            datac = data.coarsen(x=4, y=4, boundary='exact', side='left').mean()
            del data
            data = datac.values

        #print('min', np.nanmin(data))
        #print('max', np.nanmax(data))

        return data

def plot_2D_map(data_obj, *, plottype, mapfocus, size, pn_out, fn_out, fileformat):
    """
    Visualize data, plot on map, interactive display and/or write to image file.
    On purpose many hardcoded settings, adjusted to the use case (global ParFlow IHM).

    Parameters
    ----------
    data_obj : simAuxData
        Instance of the simAuxData dataclass containing `sim_data`, `mask_data`, `permafrost_data`, 
        `permafrost_lon`, `permafrost_lat`, `karst_data`, `lon`, `lat` as NumPy arrays.
    plottype : {'Sr', 'WTD', 'Reff', 'hydrofacies'}
        Type of plot to produce / variable to visualize. Affects automatically the type of custom 
        legend, colormap, labelling, etc. as this is not a generic plotting routine. 
        Options are:

        - 'Sr': for relative saturation [-]
        - 'WTD': for water table depth [m]
        - 'Reff': for effective recharge [mm year-1] 
        - 'hydrofacies': for hydrofacies distribution [hydrofacies class]

    mapfocus : {'global', 'europe', 'southAmerica1', 'southAmerica2'}
        Area of the globe to be shown. This is purely done by adjusting the map extent. The extent
        can be arbitrary, but some good looking subsets have been predefined.
        Options are:
    
        - 'global': Selfexplanatory.
        - 'europe': Selfexplanatory.
        - 'southAmerica1': Large parts of Africa included. 
        - 'southAmerica2': Very little of Africa included. 

    size: {'pagewidth', 'fullres'}
        The size of the plot (PDF or PNG file).
        Options are:
    
        - 'pagewidth': Produces an image file which fits onto a letter or A4 prtrait page, full
          linewidth. Normal use case, in this case also resampling is needed. 
        - 'fullres': Creates a huge plot, to be printed on a plotter height is matching the with of 
          the plotter paper. 

    pn_out : str
        Pathname of the output image file.
    fn_out : str
        Filename of the output image file.
    fileformat : {'pdf', 'png'}
        Matplotlib file format, recommended is pdf, related to resolution, plot size, etc. All 
        available Matplotlib file formats are OK.
        Options are:
    
        - 'pdf': PDF.
        - 'png': PNG.

    Returns
    -------
    Path
        Path object pointing to the saved image file.
    """

    # sim data, ParFlow sim data
    match plottype:
        case "Sr":
            data = data_obj.sim_var1_data[:, :]
        case "WTD":
            data = data_obj.sim_var2_data[:, :]
        case "Reff":
            data = data_obj.forcing_var1_data[:, :]
        case "hydrofacies":
            data = data_obj.extpar_var1_data[:, :]
        case _:
            print("plottype does not exist, exiting script")
            sys.exit()
    # coordinates
    data_lon = mask2_lon = mask4_lon = data_obj.sim_lon[:]
    data_lat = mask2_lat = mask4_lat = data_obj.sim_lat[:]
    # permafrost, ISIMIP data, coarse resolution
    mask1 = data_obj.mask_permafrost_data[:, :]
    mask1_lon = data_obj.mask_permafrost_lon[:] 
    mask1_lat = data_obj.mask_permafrost_lat[:] 
    # land ocean lake, derived from ParFlow external paremeter files
    mask3 = data_obj.mask_landsea_data[:, :]
    # glacier, derived from shapefile, here for legacy and consistency
    mask2 = mask3[:, :]
    # karst rock 
    mask4 = data_obj.mask_karst_data[:, :]
    #sys.exit()

    # interactive plotting
    plt.switch_backend('TkAgg')

    # regular rectilinear lon lat grid
    crs_data = crs_map = ccrs.PlateCarree()

    # ends up with plot just below 20cm in width
    fig1 = plt.figure(figsize=(9.5, 9.5))

    # higher resolution map also shows islands, scales 10m 50m 110m
    ax1 = plt.subplot(111, projection=crs_map)
    ax1.set_aspect('equal')
    ax1.set_global()
    ax1.coastlines(resolution="10m", color='gray', linewidth=0.075)  # gray
    match mapfocus:
        case "southAmerica2" | "europe":
            ax1.add_feature(cfeature.OCEAN, facecolor='midnightblue')
    #lakes = cfeature.LAKES.with_scale('10m')  # options: '110m', '50m', '10m'
    #ax1.add_feature(lakes, edgecolor='darkgreen', linewidth=0.075)
    #ax1.add_feature(cfeature.RIVERS.with_scale('10m'), edgecolor='blue', linewidth=0.2)
 
    # minimum number of gridlines
    match mapfocus:
        case "global" | "southAmerica1":
            gl_main = ax1.gridlines(draw_labels=False, linewidth=0.1, color="dimgray", linestyle="-", alpha=0.5, zorder=50)
            gl_main.xlocator = FixedLocator([0])  # Greenwich
            gl_main.ylocator = FixedLocator([0])  # Equator
            gl_extra = ax1.gridlines(draw_labels=False, linewidth=0.1, color="dimgray", linestyle="--", zorder=50)
            gl_extra.xlocator = FixedLocator([-180, -90, 90, 180])
            gl_extra.ylocator = FixedLocator([-45, 45])

    # thin map border
    ax1.spines['geo'].set_linewidth(0.5)

    #ax1.set_title('ParFlow IHM (GPU), 43200x17400@0.00833deg (approx. 1km), resampled to 4x4km^2, karst orig (aggr>0.01) Gleeson et al. (2014) porosity=0.06 + WHYMAP WOKAM 1+2', fontsize=6)

    match plottype:
        case "Sr":
            cmapDiscr = plt.get_cmap('Spectral', 50)
            levelsVals = np.linspace(0, 1, 51) # linear / logarithmic
            norm = BoundaryNorm(levelsVals, ncolors=cmapDiscr.N, clip=True)
        case "WTD":
            cmapDiscr = plt.get_cmap('Spectral_r', 50)
            levelsVals = [0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 
                          1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0, 
                          3.5,4.5,5.0,5.5,6,
                          7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                          25,30,35,40,45,50]  # non-linear custom, need to align with ParFlow dz
            norm = BoundaryNorm(levelsVals, ncolors=cmapDiscr.N, clip=True)
        case "Reff":
            levelsVals = np.logspace(-1, 3, num=50, base=10.0)
            N = len(levelsVals) - 1
            cmapDiscr = plt.get_cmap('turbo_r', N+2)
            norm = BoundaryNorm(levelsVals, ncolors=N+2, clip=False)
        case "hydrofacies":
            cmapDiscr = plt.get_cmap('gist_ncar_r', 22)
            levelsVals = np.arange(23)
            norm = BoundaryNorm(levelsVals, ncolors=cmapDiscr.N, clip=True)
        case _:
            print("plottype does not exist, exiting script")
            sys.exit()
    #print(len(levelsVals))

    # sim results plotting
    # use lon lat from indicator, this is precise 
    # masked with land sea mask
    data[mask3 == 0.] = np.nan
    # contourf and pcolormesh are no option, take way too long
    # values below / above vmin / vmax are assigned to the lowest / highest color
    plt_data_grid = ax1.imshow(data, 
        extent=[mask2_lon.min(), mask2_lon.max(), 
                mask2_lat.min(), mask2_lat.max()], 
        origin="lower", transform=crs_data,
        cmap=cmapDiscr,
        visible=False,
        norm=norm)  # arbitrary  
    #     #vmin=levelsVals[0], vmax=levelsVals[-1])  # linear
    #     #norm=LogNorm(vmin=levelsVals[0]+0.01, vmax=levelsVals[-1]))  # logarithmic
    plt_data_grid.set_rasterized(True) # does not make it much better but also not worse

    match plottype:
        case "Sr" | "WTD":
            # mask permafrost form ISIMIP
            #mask1_gt50 = np.ma.masked_where(mask1 <= 50, mask1)
            plt_mask_permafrost_grid = ax1.imshow(mask1, extent=[mask1_lon.min(), mask1_lon.max(), mask1_lat.min(), mask1_lat.max()], 
                                       origin="lower", transform=crs_data, alpha=0.5, cmap=ListedColormap(['none', (0.37, 0.25, 0.48)]), zorder=79, vmin=0.7, vmax=100)  # cadetblue
            ##plt.rcParams['hatch.linewidth'] = 0.10  # 0.2
            ##plt_mask_permafrost_cont = ax1.contourf(mask1_lon, mask1_lat, mask1_gt50, levels=[50, 100], #levels=[0.5, 1.5], 
            ##                           origin="upper", transform=crs_data, zorder=71, colors="none", hatches=["XXXXXXXXXXXXXXX"], alpha=0) #, 7 10 / linewidth=0.01)
            ##                           #origin="upper", transform=crs_data, zorder=71, colors="none", hatches=["/////////////////"], alpha=0) #, 7 10 / linewidth=0.01)
            ##for collection in plt_mask_permafrost_cont.collections:
            ##    collection.set_edgecolor("black")  # this sets hatch line color # dimgrey
            ##    #collection.set_linewidth(0.05)   # hatch thickness
            #plt_mask_permafrost_cont = ax1.contourf(mask1_lon, mask1_lat, mask1, levels=[50, 100],
            #                           origin="upper", transform=crs_data, zorder=71, colors="#c7dbe6", alpha=0.7)

            # mask glacier from somewhere, coarse resolution perhaps, fit in to ISIMIP
            shapes = shapereader.Reader('/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources/masking_glaciers_NaturalEarthData/o.data/ne_10m_glaciated_areas.shp')
            shape_feature = cfeature.ShapelyFeature(shapes.geometries(), crs_data, 
                            edgecolor='deepskyblue', linewidth=0.01, facecolor='deepskyblue', alpha=0.9, zorder=80)    # snow
            ax1.add_feature(shape_feature)

            # mask carstic regions from indicator
            # either overplot with transparency -> large areas
            # or plot as scattered points, but then need thin out much further -> too many dots otherwise
            ###plt_mask_karst_grid = ax1.imshow(mask4, extent=[mask4_lon.min(), mask4_lon.max(), mask4_lat.min(), mask4_lat.max()], 
            ###                      origin="lower", transform=crs_data, alpha=0.8, cmap=ListedColormap(['none', 'slategrey']), zorder=80)
            # ax.scatter(lon, lat, color='red', edgecolor='black', s=100, marker='*', transform=crs_data, zorder=5)  # needs more handling of lon and lat

            # add karstifiable rocks from WHYMAP WOKAM shapefile
            # use the same color coding for the polygons as on the original WHYMAP WOKAM map
            shapes_whymap = shapereader.Reader('/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources/masking_additional-data-WHYMAP-karst_BGR/o.data/shp/whymap_karst__v1_poly.shp')
            # rock type = 1 carbonate rock continuous
            filtered_geometries1 = [rec.geometry for rec in shapes_whymap.records() if rec.attributes.get("rock_type") == 1]
            custom_feature1 = cfeature.ShapelyFeature(filtered_geometries1, crs_data, edgecolor=[(0.42, 0.60, 0.83)], linewidth=0.01, facecolor=[(0.42, 0.60, 0.83)], alpha=1.0, zorder=78)
            ax1.add_feature(custom_feature1)
            # rock type = 2 carbonate rock discontinuous
            filtered_geometries2 = [rec.geometry for rec in shapes_whymap.records() if rec.attributes.get("rock_type") == 2]
            custom_feature2 = cfeature.ShapelyFeature(filtered_geometries2, crs_data, edgecolor=[(0.67, 0.85, 0.89)], linewidth=0.01, facecolor=[(0.67, 0.85, 0.89)], alpha=1.0, zorder=78)
            ax1.add_feature(custom_feature2)
            # rock type = 3 continuous evaporite rocks 
            filtered_geometries3 = [rec.geometry for rec in shapes_whymap.records() if rec.attributes.get("rock_type") == 3]
            custom_feature3 = cfeature.ShapelyFeature(filtered_geometries3, crs_data, edgecolor=[(0.90, 0.55, 0.71)], linewidth=0.01, facecolor=[(0.90, 0.55, 0.71)], alpha=1.0, zorder=78)
            ax1.add_feature(custom_feature3)
            # rock type = 4 discontinuous evporite rocks
            filtered_geometries4 = [rec.geometry for rec in shapes_whymap.records() if rec.attributes.get("rock_type") == 4]
            custom_feature4 = cfeature.ShapelyFeature(filtered_geometries4, crs_data, edgecolor=[(0.60, 0.85, 0.55)], linewidth=0.01, facecolor=[(0.60, 0.85, 0.55)], alpha=1.0, zorder=78)
            ax1.add_feature(custom_feature4)
            # rock type = 5 mixed carbonate and evaporite rocks
            filtered_geometries5 = [rec.geometry for rec in shapes_whymap.records() if rec.attributes.get("rock_type") == 5]
            custom_feature5 = cfeature.ShapelyFeature(filtered_geometries5, crs_data, edgecolor=[(0.98, 0.86, 0.49)], linewidth=0.01, facecolor=[(0.98, 0.86, 0.49)], alpha=1.0, zorder=78)
            ax1.add_feature(custom_feature5)

            # mask glacier from indicator file
            # does not look OK, too small glaciated area, overlay rather needs to be a bit bigger to show on the world map 
            # fits wel when overlaid with the Natural Earth data though
            #plt_mask_glacier_grid = ax1.imshow(mask2, extent=[mask2_lon.min(), mask2_lon.max(), mask2_lat.min(), mask2_lat.max()], 
            #                        origin="lower", transform=crs_data, alpha=0.8, cmap=ListedColormap(['none', 'purple']), zorder=99)

            # check land ocean lakes mask, can always be overplotted
            # mask is binary [0,1]
            #plt_mask_ocean_grid = ax1.imshow(mask3, extent=[mask2_lon.min(), mask2_lon.max(), mask2_lat.min(), mask2_lat.max()], origin="lower", transform=crs_data, alpha=0.5, cmap=ListedColormap(['none', 'red']) ) #, vmin=0.0, vmax=0.01)

            #stamen_terrain = cimgt.Stamen('terrain-background')
            #ax1.add_image(stamen_terrain, zoom=1, extent=[-180, 180, -60, 85])
            #img_relief = ax1.background_img(name='natural-earth', resolution='low')
            ##img_relief = ax1.stock_img()
            ##img_relief.set_cmap(cm.gray)
            #ax1.stock_img()

    # map extend
    match mapfocus:
        case "global":
            plt.xlim(-180, 180)
            plt.ylim(-60, 85)
        case "europe":
            plt.xlim(-25, 40)
            plt.ylim(30, 72)
        case "southAmerica1":
            plt.xlim(-95, 32)
            plt.ylim(-51, 15)
        case "southAmerica2":
            plt.xlim(-85, -9)
            plt.ylim(-36, 14)

    # thin colorbar
#    match plottype:
#        case "Sr":
#            cb = plt.colorbar(plt_data_grid, ax=ax1, extend='neither', pad=0.007, shrink=0.2, drawedges=False, 
#                            orientation='horizontal', ticks=np.linspace(0, 1, 11))
#            cb.set_ticklabels(['0.0','\ndry','','','','0.5','','','','\nwet','1.0'])
#            cb.set_label('Relative saturation [-]', fontsize=6, fontweight='bold')
#        case "WTD":
#            cb = plt.colorbar(plt_data_grid, ax=ax1, extend='max', pad=0.007, shrink=0.2, drawedges=False, 
#                             orientation='horizontal', ticks=[1, 3, 10, 50])
#            cb.set_ticklabels(['1', '3', '<10\nshallow', '50\ndeep'])
#            match mapfocus:
#                case "global" | "southAmerica1":
#                    cb.set_label('Water table depth [m]', fontsize=6, fontweight='bold')
#                case "southAmerica2" | "europe":
#                    cb.ax.tick_params(color="whitesmoke", labelcolor="whitesmoke")
#                    cb.set_label('Water table depth [m]\nParFlow IHM, 1km global', fontsize=6, fontweight='bold', color='whitesmoke')
#                    cb.outline.set_edgecolor("whitesmoke")
#        case "Reff":
#            cb = plt.colorbar(plt_data_grid, ax=ax1, extend='both', pad=0.007, shrink=0.2, boundaries=levelsVals, spacing='uniform',
#                            orientation='horizontal') #, format=LogFormatterSciNotation(base=10, labelOnlyBase=True)) 
#            cb.set_ticks([1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03])
#            cb.set_ticklabels(['10$^{-1}$','10$^{0}$','10$^{1}$','10$^{2}$','10$^{3}$'])
#            cb.set_label('Effective recharge [mm year$^{-1}$]', fontsize=6, fontweight='bold')
#        case "hydrofacies":
#            cb = plt.colorbar(plt_data_grid, ax=ax1, extend='neither', pad=0.007, shrink=0.2, drawedges=True, 
#                           orientation='horizontal', ticks=[0.5, 3.5, 6.5, 9.5, 12.5, 15.5, 18.5, 21.5])
#            cb.set_ticklabels(['1', '4', '7', '10', '13', '16', '19', '22'])
#            cb.set_label('Hydrofacies, bottom layer', fontsize=6, fontweight='bold')
#        case _:
#            print("plottype does not exist, exiting script")
#            sys.exit()
#    cb.outline.set_linewidth(0.3)
#    cb.ax.tick_params(labelsize=6)
#
#    # re-position the color bar
#    pos = cb.ax.get_position()
#    match mapfocus:
#        case "global":
#            match plottype:
#                case "Sr" | "WTD":
#                    new_pos = [pos.x0 - 0.3, pos.y0 + 0.075, pos.width, pos.height]
#                case "Reff" | "hydrofacies":
#                    new_pos = [pos.x0 - 0.3, pos.y0 + 0.09, pos.width, pos.height]
#        case "southAmerica1":
#            new_pos = [pos.x0 + 0.05, pos.y0 + 0.098, pos.width, pos.height]
#        case "southAmerica2":
#            new_pos = [pos.x0 + 0.25, pos.y0 + 0.098, pos.width, pos.height]
#        case "europe":
#            new_pos = [pos.x0 - 0.3, pos.y0 + 0.125, pos.width, pos.height]
#    cb.ax.set_position(new_pos)

    # add sub-figure label
    #ax1.text(0.015, 0.96,'a)',transform=ax1.transAxes,fontsize=8,fontweight="bold",va="top", ha="left")

    # add custom legends, permafrost, karst, glacier
    # width/height in degrees
    match mapfocus:
        case "global":
            lonstart = -175.
            latstart =  -6. #-10.
            boxheight = 3.5
            boxwidth = 7.
            txtlonstart = 7.8
            txtlatstart = 1.6
            latoffset = -6
        case "southAmerica1":
            lonstart = -36.
            latstart = -28.
            boxheight = 1.8
            boxwidth = 4.5
            txtlonstart = 5
            txtlatstart = 0.8
            latoffset = -3

    match mapfocus:
        case "global" | "southAmerica1":
            match plottype:
                case "Sr" | "WTD":

                    print("adding custom legend for permafrost, karst, glacier")

                    # glaciers
                    rect0 = Rectangle((lonstart, latstart + latoffset * 0.), boxwidth, boxheight, facecolor="deepskyblue", alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect0)
                    rect0.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 0. + txtlatstart, "Glaciated regions", va="center", ha="left", fontsize=6)

                    # permafrost
                    rect1 = Rectangle((lonstart, latstart + latoffset * 1.), boxwidth, boxheight, facecolor=(0.37, 0.25, 0.48), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect1)
                    rect1.set_zorder(100)
                #    #rect11 = Rectangle((lonstart, latstart + latoffset * 1.), boxwidth, boxheight, facecolor="dodgerblue", alpha=0.3, edgecolor="black", linewidth=0.3)
                #    #ax1.add_patch(rect11)
                #    #rect11.set_zorder(101)
                #    rect12 = Rectangle((lonstart, latstart + latoffset * 1.), boxwidth, boxheight, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.2, hatch='XXXXXXXXX') # dimgrey
                #    #rect12 = Rectangle((lonstart, latstart + latoffset * 1.), boxwidth, boxheight, facecolor="none", alpha=1.0, edgecolor="dimgrey", linewidth=0.3, hatch='///////')
                #    ax1.add_patch(rect12)
                #    rect12.set_zorder(102)
                #    rect13 = Rectangle((lonstart, latstart + latoffset * 1.), boxwidth, boxheight, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.3)
                #    ax1.add_patch(rect13)
                #    rect13.set_zorder(103)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 1. + txtlatstart, "Permafrost", va="center", ha="left", fontsize=6)

                    # karst
                #    rect2 = Rectangle((lonstart, latstart + latoffset * 2.), boxwidth, boxheight, facecolor="white", alpha=1.0, edgecolor="black", linewidth=0.3)
                #    ax1.add_patch(rect2)
                #    rect2.set_zorder(100)
                #    rect21 = Rectangle((lonstart, latstart + latoffset * 2.), boxwidth, boxheight, facecolor="slategrey", alpha=0.8, edgecolor="black", linewidth=0.3)
                #    ax1.add_patch(rect21)
                #    rect21.set_zorder(101)
                #    rect22 = Rectangle((lonstart, latstart + latoffset * 2.), boxwidth, boxheight, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.3)
                #    ax1.add_patch(rect22)
                #    rect22.set_zorder(102)
                #    ax1.text(lonstart + txtlonstart, latstart + latoffset * 2. + txtlatstart, "Karst rock (Gleeson et al., 2014; porosity=0.06)", va="center", ha="left", fontsize=6)
                    rect2 = Rectangle((lonstart, latstart + latoffset * 2.), boxwidth, boxheight, facecolor=(0.42, 0.60, 0.83), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect2)
                    rect2.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 2. + txtlatstart, "Carbonate rocks, continuous", va="center", ha="left", fontsize=6)
                    
                    rect3 = Rectangle((lonstart, latstart + latoffset * 3.), boxwidth, boxheight, facecolor=(0.67, 0.85, 0.89), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect3)
                    rect3.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 3. + txtlatstart, "Carbonate rocks, discontinuous", va="center", ha="left", fontsize=6)

                    rect4 = Rectangle((lonstart, latstart + latoffset * 4.), boxwidth, boxheight, facecolor=(0.90, 0.55, 0.71), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect4)
                    rect4.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 4. + txtlatstart, "Evaporite rocks, continuous", va="center", ha="left", fontsize=6)

                    rect5 = Rectangle((lonstart, latstart + latoffset * 5.), boxwidth, boxheight, facecolor=(0.60, 0.85, 0.55), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect5)
                    rect5.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 5. + txtlatstart, "Evaporite rocks, discontinuous", va="center", ha="left", fontsize=6)

                    rect6 = Rectangle((lonstart, latstart + latoffset * 6.), boxwidth, boxheight, facecolor=(0.98, 0.86, 0.49), alpha=1.0, edgecolor="black", linewidth=0.3)
                    ax1.add_patch(rect6)
                    rect6.set_zorder(100)
                    ax1.text(lonstart + txtlonstart, latstart + latoffset * 6. + txtlatstart, "Mixed carbonate and evaporite rocks", va="center", ha="left", fontsize=6)

    # 300 600 1200 2400
    #os.makedirs(pn_out, exist_ok=True)
    fig1.savefig(pn_out+"/"+fn_out+"."+fileformat, bbox_inches='tight', pad_inches=0.1, dpi=1200)

    # for interactive checking and zooming into the data down to full resolution
    #plt.show()

    return pn_out+"/"+fn_out+"."+fileformat

def main():

    time_start = time.time()

    dir_data = "/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources"

    # create dataclass holding all relevant datasets, for code 
    data_obj = simAuxData(dir_data+"/"+"plotdata_manuscript_SKo_links", "satur1x1.nc", "saturation", 1,
                      dir_data+"/"+"plotdata_manuscript_SKo_links", "wtd1x1.nc", "wtd", 1,
                      dir_data+"/"+"plotdata_manuscript_SKo_links", "output_mask1x1.nc", "mask", 1,
                      dir_data+"/"+"masking_permafrost-ESA-CCI_CEDA/p.data.mean_permafrost_distribution", "ESACCI-PERMAFROST-L4-PFR-1997to2023_float_mv_cat_compr4_timmean_compr4.nc", "PFR", 1,
                      dir_data+"/"+"masking_carst-glaciated_indicator-gleeson", "global_gleeson_porosity_008333.nc", "Gleeson", 1,
                      dir_data+"/"+"plotdata_manuscript_SKo_links", "meanforcing.nc", "evaptrans", 1,
                      dir_data+"/"+"plotdata_manuscript_SKo_links", "indicator.nc", "indicator", 1)

    time_intermediate = time.time()
    print('exec wallclock time reading and processing =  %0.3f s' % (time_intermediate - time_start)) 

    #sys.exit()

    # add: , interactiveplot=False
    # manually activate / deactivate
    #pn_fn_image_output = plot_2D_map(data_obj, plottype="hydrofacies", mapfocus="global", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_Fig1_hydrofacies", fileformat="pdf")

    #pn_fn_image_output = plot_2D_map(data_obj, plottype="Reff", mapfocus="global", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_Fig2_Reff", fileformat="pdf")
    
    pn_fn_image_output = plot_2D_map(data_obj, plottype="Sr", mapfocus="global", size="pagewidth",
        pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_MaterialMethods-AddMapLayers_r2", fileformat="pdf")
    
    #pn_fn_image_output = plot_2D_map(data_obj, plottype="WTD", mapfocus="global", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_Fig4a_WTD", fileformat="pdf")
    
    #pn_fn_image_output = plot_2D_map(data_obj, plottype="WTD", mapfocus="southAmerica1", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_Fig4b_WTD_SouthAmerica", fileformat="pdf")
    
    #pn_fn_image_output = plot_2D_map(data_obj, plottype="WTD", mapfocus="southAmerica2", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_poster_WTD_SouthAmerica", fileformat="pdf")
    
    #pn_fn_image_output = plot_2D_map(data_obj, plottype="WTD", mapfocus="europe", size="pagewidth",
    #    pn_out="./", fn_out="JHX_KolletEtAl_ParFlowGlobal_website_WTD_Europe", fileformat="pdf")

    print('exec wallclock time plotting  =  %0.3f s' % (time.time() - time_intermediate)) 

if __name__ == '__main__':

    main()
