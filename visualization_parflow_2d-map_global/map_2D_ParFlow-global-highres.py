#!/usr/bin/env python
# coding: utf-8

'''
Visualisation of global ParFlow model input and output on a map.

Purpose:
This script combines the visualisation for the different figures in 
Kollet et al. "Global groundwater modeling: Proof-of-concept of 3D 
variably saturated flow simulation at km-resolution" in JHX for further
reuse and adjustment. 

Inputs:
- Locally stored compilation of global datasets.
- ParFlow simulation results.
- ParFlow external parameter files.
- Additional data from external sources: glacier mask, permafrost mask.

Ouputs:
- Image files

Notes:
Plot size, data set size in grid elements, map scale, dpi, image file format are
adjusted with each other. Also the viewer makes a difference.

Limitations:
Currently only Fig3 from Kollet et al., plottype='Sr', is implemented, other
plots will follow. 

Operating environment:
- locally with conda forge env (tested), or
- JSC/JURECA frontnode (tested) (JSC/JUWELS should also work), source this environment initialisation file:
  https://icg4geo.icg.kfa-juelich.de/Organisation/hpc_it_dev_data_howtos-and-bestpractice/hpc_scientific_operating_environments/-/blob/main/loadenv.JURECA-DC_JUWELS_2023_GCC-OpenMPI_Std-Py-AI.ini`
'''

import time, sys

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

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader

from dataclasses import dataclass, field


__authors__ = "Klaus GOERGEN"
__orcid__ = "https://orcid.org/0000-0002-4208-3444"
__email__ = ""
__maintainer__ = "Klaus GOERGEN"
__copyright__ = "Copyright (c) 2025, https://www.fz-juelich.de (FZJ)"
__license__ = "MIT"
__version__ = "1.0.0"
__date__ = "2025-10-09"
__status__ = "Beta"
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

    sim_dir: str
    sim_file: str
    sim_varname: str
    sim_resample: int
    sim_data: np.ndarray = field(init=False)
    mask_dir: str
    mask_file: str
    mask_varname: str
    mask_resample: int
    mask_data: np.ndarray = field(init=False)
    permafrost_dir: str
    permafrost_file: str
    permafrost_varname: str
    permafrost_data: np.ndarray = field(init=False)
    permafrost_lon: np.ndarray = field(init=False)
    permafrost_lat: np.ndarray = field(init=False)
    karst_dir: str
    karst_file: str
    karst_varname: str
    karst_resample: int
    karst_data: np.ndarray = field(init=False)
    lon: np.ndarray = field(init=False)
    lat: np.ndarray = field(init=False)

    def __post_init__(self):
        self.sim_data = self._read_sim_data_proc(self.sim_dir, self.sim_file, self.sim_varname, self.sim_resample)
        self.mask_data = self._read_mask(self.mask_dir, self.mask_file, self.mask_varname, self.mask_resample)
        self.permafrost_data, self.permafrost_lon, self.permafrost_lat = self._read_permafrost(self.permafrost_dir, self.permafrost_file, self.permafrost_varname)
        self.karst_data, self.lon, self.lat = self._read_karst(self.karst_dir, self.karst_file, self.karst_varname, self.karst_resample)
  
    def _read_sim_data_proc(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        #print('----------------------------------------')
        ds = xr.open_dataset(pn+"/"+fn)    
        if regr == 0:
            data = ds[varname].values
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

    def _read_mask(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

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

    def _read_permafrost(self, pn: str, fn: str, varname: str) -> np.ndarray:

        #print('----------------------------------------')
        ds_mask = xr.open_dataset(pn+"/"+fn, decode_times=False)
        mask = ds_mask[varname].isel(time=0).values
        #print('mask isimip: ', mask.shape)

        lon =  ds_mask['lon'].values
        lat =  ds_mask['lat'].values
        #print('mask isimip lon:', lon.shape)
        #print('mask isimip lat:', lat.shape)

        #if mtype == 'permafrost':
        #    print('permafrost')
        #    mask = np.ma.masked_where(mask == 4, mask)

        return mask, lon, lat

    def _read_karst(self, pn: str, fn: str, varname: str, regr: int) -> np.ndarray:

        #print('----------------------------------------')
        ds_mask = xr.open_dataset(pn+"/"+fn, decode_times=False)
        mask = ds_mask[varname].isel(time=0) #.values
        lon =  ds_mask['lon']
        lat =  ds_mask['lat']

        #print('karst')
        # the NaNs are ignored
        mask[:] = np.where(mask == 0.06, 1, 0)

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
            mask[mask < 0.75] = 0

        return mask, lon, lat

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
    data = data_obj.sim_data[:, :]
    data_lon = mask2_lon = mask4_lon = data_obj.lon[:]
    data_lat = mask2_lat = mask4_lat = data_obj.lat[:]
    # permafrost, ISIMIP data, coarse resolution
    mask1 = data_obj.permafrost_data[:, :]
    mask1_lon = data_obj.permafrost_lon[:] 
    mask1_lat = data_obj.permafrost_lat[:] 
    # land ocean lake, derived from ParFlow external paremeter files
    mask3 = data_obj.mask_data[:, :]
    # glacier, derived from shapefile, here for legacy and consistency
    mask2 = mask3[:, :]
    # karst rock 
    mask4 = data_obj.karst_data[:, :]
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
    ax1.coastlines(resolution="10m", color='gray', linewidth=0.075)
 
    # minimum number of gridlines
    gl_main = ax1.gridlines(draw_labels=False, linewidth=0.1, color="dimgray", linestyle="-", alpha=0.5)
    gl_main.xlocator = FixedLocator([0])  # Greenwich
    gl_main.ylocator = FixedLocator([0])  # Equator
    gl_extra = ax1.gridlines(draw_labels=False, linewidth=0.1, color="dimgray", linestyle="--")
    gl_extra.xlocator = FixedLocator([-180, -90, 90, 180])
    gl_extra.ylocator = FixedLocator([-45, 45])

    # thin map border
    ax1.spines['geo'].set_linewidth(0.5)

    ax1.set_title('ParFlow IHM (GPU), 43200x17400@0.00833deg (approx. 1km), resampled to 4x4km^2', fontsize=7)

    cmapDiscr = plt.get_cmap('Spectral', 50)
    levelsVals = np.linspace(0, 1, 51) # linear / logarithmic
    #levelsVals = [0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0, 
    #              1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0, 
    #              3.5,4.5,5.0,5.5,6,
    #              7,8,9,10,11,12,13,14,15,16,17,18,19,20,
    #              25,30,35,40,45,50]  # non-linear custom, need to align with ParFlow dz
    #print(len(levelsVals))
    norm = BoundaryNorm(levelsVals, ncolors=cmapDiscr.N, clip=True)

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
        norm=norm)  # arbitrary  
    #     #vmin=levelsVals[0], vmax=levelsVals[-1])  # linear
    #     #norm=LogNorm(vmin=levelsVals[0]+0.01, vmax=levelsVals[-1]))  # logarithmic
    plt_data_grid.set_rasterized(True) # does not make it much better but also not worse

    # mask permafrost form ISIMIP
    #plt_mask_permafrost_grid = ax1.imshow(mask1, extent=[mask1_lon.min(), mask1_lon.max(), mask1_lat.min(), mask1_lat.max()], 
    #                           origin="upper", transform=crs_data, alpha=0.3, cmap=ListedColormap(['none', 'dodgerblue']), zorder=70)  #, vmin=0.5, vmax=1)  # better
    plt.rcParams['hatch.linewidth'] = 0.10  # 0.2
    plt_mask_permafrost_cont = ax1.contourf(mask1_lon, mask1_lat, mask1, levels=[0.5, 1.5], 
                               origin="upper", transform=crs_data, zorder=71, colors="none", hatches=["XXXXXXXXXXXXXXX"], alpha=0) #, 7 10 / linewidth=0.01)
                               #origin="upper", transform=crs_data, zorder=71, colors="none", hatches=["/////////////////"], alpha=0) #, 7 10 / linewidth=0.01)
    for collection in plt_mask_permafrost_cont.collections:
        collection.set_edgecolor("black")  # this sets hatch line color # dimgrey
        #collection.set_linewidth(0.05)   # hatch thickness

    # mask glacier from somewhere, coarse resolution perhaps, fit in to ISIMIP
    shapes = shapereader.Reader('/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources/masking_glaciers_NaturalEarthData/ne_10m_glaciated_areas/ne_10m_glaciated_areas.shp')
    shape_feature = cfeature.ShapelyFeature(shapes.geometries(), crs_data, 
                    edgecolor='snow', linewidth=0.01, facecolor='snow', alpha=1.0, zorder=98)
    ax1.add_feature(shape_feature)

    # mask carstic regions from indicator
    # either overplot with transparency -> large areas
    # or plot as scattered points, but then need thin out much further -> too many dots otherwise
    plt_mask_karst_grid = ax1.imshow(mask4, extent=[mask4_lon.min(), mask4_lon.max(), mask4_lat.min(), mask4_lat.max()], 
                          origin="lower", transform=crs_data, alpha=0.6, cmap=ListedColormap(['none', 'slategrey']), zorder=80)
    # ax.scatter(lon, lat, color='red', edgecolor='black', s=100, marker='*', transform=crs_data, zorder=5)  # needs more handling of lon and lat

    # mask glacier from indicator file
    # does not look OK, too small glaciated area, overlay rather needs to be a bit bigger to show on the world map 
    # fits wel when overlaid with the Natural Earth data though
    #plt_mask_glacier_grid = ax1.imshow(mask2, extent=[mask2_lon.min(), mask2_lon.max(), mask2_lat.min(), mask2_lat.max()], 
    #                        origin="lower", transform=crs_data, alpha=0.8, cmap=ListedColormap(['none', 'purple']), zorder=99)

    # check land ocean lakes mask, can always be overplotted
    # mask is binary [0,1]
    #plt_mask_ocean_grid = ax1.imshow(mask3, extent=[mask2_lon.min(), mask2_lon.max(), mask2_lat.min(), mask2_lat.max()], origin="lower", transform=crs_data, alpha=0.5, cmap=ListedColormap(['none', 'red']) ) #, vmin=0.0, vmax=0.01)

    # map extend
    plt.xlim(-180, 180)
    plt.ylim(-60, 85)

    # thin colorbar
    cb = plt.colorbar(plt_data_grid, ax=ax1, extend='neither', pad=0.007, shrink=0.2, drawedges=False, 
                       orientation='horizontal', ticks=np.linspace(0, 1, 11))
    cb.ax.tick_params(labelsize=6)
    cb.set_ticklabels(['0.0','\ndry','','','','0.5','','','','\nwet','1.0'])
    cb.set_label('Relative saturation [-]', fontsize=6, fontweight='bold')
    cb.outline.set_linewidth(0.3)
    # re-position the color bar
    pos = cb.ax.get_position()
    new_pos = [pos.x0 - 0.3, pos.y0 + 0.075, pos.width, pos.height]
    cb.ax.set_position(new_pos)

    # add sub-figure label
    #ax1.text(0.015, 0.96,'a)',transform=ax1.transAxes,fontsize=8,fontweight="bold",va="top", ha="left")

    # add custom legends, permafrost, karst, glacier
    # width/height in degrees
    rect0 = Rectangle((-175, -10), 7, 3.5, facecolor="snow", alpha=1.0, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect0)
    rect0.set_zorder(100)
    ax1.text(-175 + 7.8, -10 + 1.6, "Glaciated regions", va="center", ha="left", fontsize=6)

    rect1 = Rectangle((-175, -16), 7, 3.5, facecolor="white", alpha=1.0, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect1)
    rect1.set_zorder(100)
    #rect11 = Rectangle((-175, -16), 7, 3.5, facecolor="dodgerblue", alpha=0.3, edgecolor="black", linewidth=0.3)
    #ax1.add_patch(rect11)
    #rect11.set_zorder(101)
    rect12 = Rectangle((-175, -16), 7, 3.5, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.2, hatch='XXXXXXXXX') # dimgrey
    #rect12 = Rectangle((-175, -16), 7, 3.5, facecolor="none", alpha=1.0, edgecolor="dimgrey", linewidth=0.3, hatch='///////')
    ax1.add_patch(rect12)
    rect12.set_zorder(102)
    rect13 = Rectangle((-175, -16), 7, 3.5, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect13)
    rect13.set_zorder(103)
    ax1.text(-175 + 7.8, -16 + 1.6, "Permafrost regions", va="center", ha="left", fontsize=6)

    rect2 = Rectangle((-175, -22), 7, 3.5, facecolor="white", alpha=1.0, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect2)
    rect2.set_zorder(100)
    rect21 = Rectangle((-175, -22), 7, 3.5, facecolor="slategrey", alpha=0.6, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect21)
    rect21.set_zorder(101)
    rect22 = Rectangle((-175, -22), 7, 3.5, facecolor="none", alpha=1.0, edgecolor="black", linewidth=0.3)
    ax1.add_patch(rect22)
    rect22.set_zorder(102)
    ax1.text(-175 + 7.8, -22 + 1.6, "Karst rock", va="center", ha="left", fontsize=6)

    # 300 600 1200 2400
    fig1.savefig(pn_out+"/"+fn_out+"."+fileformat, bbox_inches='tight', pad_inches=0.1, dpi=1200)

    # for interactive checking and zooming into the data down to full resolution
    #plt.show()

    return pn_out+"/"+fn_out+"."+fileformat

def main():

    time_start = time.time()

    dir_data = "/p/data1/slts/shared_data/collection_ParFlow-global_misc-sources"

    # create dataclass holding all relevant datasets, for code 
    data_obj = simAuxData(dir_data+"/"+"plotdata_manuscript_SKo_links", "satur1x1.nc", "saturation", 1,
                      dir_data+"/"+"plotdata_manuscript_SKo_links", "output_mask1x1.nc", "mask", 1,
                      dir_data+"/"+"masking_permafrost_isimip", "isimip2b_clm45_permafrost_mask_2005_3m.nc", "mask",
                      dir_data+"/"+"masking_carst-glaciated_indicator-gleeson", "global_gleeson_porosity_008333.nc", "Gleeson", 1)

    time_intermediate = time.time()
    print('exec wallclock time reading and processing =  %0.3f s' % (time_intermediate - time_start)) 

    #sys.exit()

    # add: , interactiveplot=False
    pn_fn_image_output = plot_2D_map(data_obj, plottype="Sr", mapfocus="global", size="pagewidth",
        pn_out="/p/project/cjjsc39/jjsc3900/tmp.visualization_parflow_2d-map_global", fn_out="test", fileformat="pdf")

    print('exec wallclock time plotting  =  %0.3f s' % (time.time() - time_intermediate)) 

if __name__ == '__main__':

    main()