#!/usr/bin/env python
# coding: utf-8

'''
Visualisation testing nd demo tool to plot data on triangular grids

Purpose:
Use this as a basis ofr further vis tools, when, e.g. ICON data is to be plotted
on its native grid. Deliberately hardcoded and many comments, need to uncomment.
Many explanations added during exploration of the fucntionalities.
'''

import time, os, sys
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature


__authors__ = "Klaus GOERGEN"
__orcid__ = "https://orcid.org/0000-0002-4208-3444"
__email__ = ""
__maintainer__ = "Klaus GOERGEN"
__copyright__ = "Copyright (c) 2025, https://www.fz-juelich.de (FZJ)"
__license__ = "MIT"
__version__ = "1.0.0"
__date__ = "2025-08-11"
__status__ = "Beta"
__credits__ = [
    "DKRZ"
    ]
__acknowledgements__ = [
    "The authors acknowledge the use of the Zonda service (https://zonda.ethz.ch/) for the preparation of external parameter data for the ICON model. Zonda has been developed and is maintained by the Center for Climate Systems Modeling (C2SM), ETH Zurich, in cooperation with the Federal Office for Meteorology and Climatology MeteoSwiss and Deutscher Wetterdienst (DWD)."
    ]
__references__ = [
    "https://docs.dkrz.de/doc/visualization/sw/python/source_code/python-matplotlib-example-unstructured-icon-triangles-plot-python-3.html",
    "https://docs.dkrz.de/doc/visualization/sw/python/source_code/python-matplotlib-triangular-grid-with-tripcolor-ICON.html",
    "https://easy.gems.dkrz.de/Processing/playing_with_triangles/tripcolor.html",
    "https://www.dwd.de/DE/leistungen/nwv_icon_tutorial/pdf_einzelbaende/icon_tutorial2024.pdf",
    "https://zonda.ethz.ch"
    ]

def main():

    # Agg: file-output
    # TkAgg: on-screen, interactive with macOS
    #plt.switch_backend('TkAgg')

    t1 = time.time()

    dsPnFn = '../data_test_geo.EUR-12_R13B05_ICON_via-zonda_v20250731orig/europe011_DOM01_external_parameter.nc'
    gridPnFn = '../data_test_geo.EUR-12_R13B05_ICON_via-zonda_v20250731orig/europe011_DOM01.nc'
    gridRotClPnFn = '../data_test_geo.EUR-12_R13B05_ICON_via-zonda_v20250731orig/europe011_DOM01_latlon_rotated.nc'
    
    ds = xr.open_dataset(dsPnFn)
    dsGrid = xr.open_dataset(gridPnFn)
    dsRotClGrid = xr.open_dataset(gridRotClPnFn)
    #print(ds)
    #print(dsGrid)
    #print(dsRotClGrid)

    data = ds['topography_c'][:].values
    print('data shape: ', data.shape)

    # get various grid specs, needed for different plotting routines
    # triangle (=grid cells) center coords, triangle vertex coords, indiv triangle vertices coords
    # determine also number of cells (= triangles) and number of vertices
    # most ICON variables (e.g., wind!) are calculated for the clon and clat loacations, at the center of each triangle face
    # in case of the CORDEX-based pan-EU grids (i.e., EUR-11/12) the icosahedron 
    # is already rotated
    clon = np.rad2deg(dsGrid.clon.values)
    clat = np.rad2deg(dsGrid.clat.values)
    vlon = np.rad2deg(dsGrid.vlon.values)
    vlat = np.rad2deg(dsGrid.vlat.values)
    clon_vertices = np.rad2deg(dsGrid.clon_vertices.values)
    clat_vertices = np.rad2deg(dsGrid.clat_vertices.values)
    ncells, nv = clon_vertices.shape[0], clon_vertices.shape[1]
    print('Cells: %6d ' % clon.size)
    print('clat shape: ', clat.shape)
    print('vlat shape: ', vlat.shape)
    print('clat_vertices shape: ', clat_vertices.shape)
    lon = dsRotClGrid.rlon.values
    lat = dsRotClGrid.rlat.values

    # define a mask, e.g., western Alps
    mask = (
    (clat > 40)
    & (clat < 50)
    & (clon > 4)
    & (clon < 9)
    )

    # vocs (vertex of cell, corner indices for each triangular grid cell)
    # counting needs to be reset from starting at 1 to 0, and transposed
    # without mask: use the complete grid
    voc = dsGrid.vertex_of_cell.T.values - 1
    #voc = dsGrid.vertex_of_cell.T[mask].values - 1
    print('voc orig shape: ', dsGrid.vertex_of_cell.values.shape)
    print('voc.shape: ', voc.shape)
    print(voc)

    # get the extend of the geolocated data (model domain) 
    # when using the rotated lon/lat, this has to come from rlon, rlat
    used_vertices = np.unique(voc)
    print('used_vertices: ', used_vertices.shape)
    #lat_min = vlat[used_vertices].min()
    #lat_max = vlat[used_vertices].max()
    #lon_min = vlon[used_vertices].min()
    #lon_max = vlon[used_vertices].max()
    lat_min = lat.min()-2
    lat_max = lat.max()+2
    lon_min = lon.min()-2
    lon_max = lon.max()+2
    print(lat_min, lat_max, lon_min, lon_max)

    # triangles option 1 
    # create the triangles for plotting, always the same with DKRZ examples
    # used with the `polycollection`
    ##clon_vertices = np.where(clon_vertices < -180., clon_vertices + 360., clon_vertices)
    ##clon_vertices = np.where(clon_vertices >  180., clon_vertices - 360., clon_vertices)
    triangles = np.zeros((ncells, nv, 2), np.float32)
    for i in range(0, ncells, 1):
        triangles[i,:,0] = np.array(clon_vertices[i,:])
        triangles[i,:,1] = np.array(clat_vertices[i,:])
    print('triangles, self-made from in main grid file infos: ', triangles.shape)

    # triangles option 2 
    # used with `tripcolor`
    # create a triangulation object "triang" to describe the triangular grid
    # then triangulation does not have to be done on-the-fly with the plotting
    # once the triangulation has been specified other functions are available
    # the triangulation object is described in the ICON main grid file
    # here the cell center points are provided and a Delauney triangulation is then used
    # results in many more triangles, does not make sense in case
    # the triangles are fully defined
    triang2 = tri.Triangulation(clon, clat)  
    print(triang2)
    print(triang2.triangles.shape)
    print(triang2.is_delaunay)
    # here triangles are specified, this is the one to be used
    triang3 = tri.Triangulation(vlon, vlat, voc)
    print(triang3)
    print(triang3.triangles.shape)
    print(triang3.is_delaunay)

    #sys.exit()

    # we need the data projection / "grid" (input to transform keyword)
    # plus the plotting map projection
    # describes the ICON grid best in the Cartopy crs selection
    # from ICON "main grid file"
    # :grid_mapping_name = "lat_long_on_sphere" ;
    # :crs_id = "urn:ogc:def:cs:EPSG:6.0:6422" ;
    # :crs_name = "Spherical 2D Coordinate System" ;
    # combine plate-carree with rotated pole to get the alignment of map and 
    # data points
    data_crs = ccrs.PlateCarree()  
    map_projection = ccrs.RotatedPole(pole_longitude=dsRotClGrid.rotated_pole.grid_north_pole_longitude, pole_latitude=dsRotClGrid.rotated_pole.grid_north_pole_latitude)
    # for testing:
    #map_projection = ccrs.PlateCarree()
    #map_projection = ccrs.Orthographic()

    fig1 = plt.figure()  #figsize=(5.75, 5.75))

    ax1 = plt.subplot(111, projection=map_projection)
    ax1.set_aspect('equal')
    ax1.coastlines()
    ax1.set_global()  # overridden by xlim and ylim, if not set, limit automatically
    gl = ax1.gridlines()  # just plot the labels for double-checking the coords used
    gl.bottom_labels = True
    gl.left_labels   = True
    gl.top_labels    = False
    gl.right_labels  = False
    #ax.coastlines(resolution='50m',linewidth=1.0)
    #ax.add_feature(cfeature.BORDERS, linestyle='-', alpha=.5, linewidth=1.0)
    #ax1.add_feature(cfeature.OCEAN, color='azure')

    ax1.set_title('ICON plotting, tests with icosahedral grid')

    # option 1
    # do this on the fly
    # super fast, no hairline artifacts between triangle plot objects
    # too simplified, but OK for global plots
    # because it relates neighbouring triangles, this plotting routine results
    # in artifacts along the LAM edges, OK inside if x/ylim are used
    #pdo = ax1.tricontourf(clon, clat, data, transform=data_crs)  
    # ? does not work:
    #pdo = ax1.tricontourf(triang3, data, transform=data_crs) 
    #pdo = ax1.tricontourf(vlon, vlat, data, triangles=voc, transform=data_crs)

    # option 2
    # tripcolor
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tripcolor.html
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/tripcolor_demo.html
    # if the triangulation is not provided -> calc on the fly based on the definition of the triangles
    # results in a polycollection from above
    # this mapping takes most of the time
    # the data is passed as a positional color parameter or by the facecolors keyword
    # this controls how the values is displayed, as a point information for the cell center / face or as the mean of the triangles corner points
    # option 2a
    # use the pre-processed traingulation information, which is based on
    # different data, trinag3 is be more accurate, trinag2 does not make sense
    # then one might as well regrid
    #pdo = ax1.tripcolor(triang2, data, transform=data_crs) 
    #pdo = ax1.tripcolor(triang2, data, transform=data_crs, shading='gouraud')
    # option 2b
    #pdo = ax1.tripcolor(triang3, data, transform=data_crs)
    #pdo = ax1.tripcolor(triang3, facecolors=data, transform=data_crs)
    #pdo = ax1.tripcolor(triang3, facecolors=data, transform=data_crs, shading='flat')  # optimum!! -> exactly showing orig data, but thin hairlines possible
    pdo = ax1.tripcolor(triang3, facecolors=data, transform=data_crs, shading='flat', rasterized=True)  # optimum!!! test RE hairline -> no hairlines, all other map elements vectorized, only the triangles are rasterized but high resolution -> control via the dpis in the output
    #pdo = ax1.tripcolor(triang3, facecolors=data, transform=data_crs, shading='flat', edgecolors='none', antialiased=False)  # test RE hairline -> no effect
    #pdo = ax1.tripcolor(triang3, data, transform=data_crs, edgecolors='k', linewidth=0.05) 
    # option 2c
    # complete dataset, w/o mask, some optimization with the plotting
    # or with a geographic mask / subset
    # works OK, same result as 4.0 with triang3
    #pdo = ax1.tripcolor(vlon, vlat, facecolors=data, triangles=voc, transform=data_crs, shading='flat', edgecolors='none')  # optimum!
    #pdo = ax1.tripcolor(vlon, vlat, ds['topography_c'][:].isel(cell=mask).values, triangles=voc, transform=data_crs)

    # option 3
    # polycollection
    # not yet working
    # https://matplotlib.org/stable/api/collections_api.html#matplotlib.collections.PolyCollection
    # one needs to set the colors right
    # option 3a
    # works, as optimum 2b and 2c
    #pdo = PolyCollection(triangles, array=data, transform=data_crs)  # optimum
    #ax1.add_collection(pdo)
    # option 3b
    # not yet implemented:
    #colors   = np.ones([ncells,4], np.float32)
    #pdo = PolyCollection(triangles, array=None, facecolors=colors, edgecolor='black', linewidth=0.05, transform=data_crs, zorder=0)
    #pdo = PolyCollection(triangles, array=None, facecolors=colors, transform=data_crs)
    
    # option 4
    # datashader
    # https://docs.dkrz.de/doc/visualization/sw/python/source_code/python-matplotlib-example-ICON-R02B09-datashader-plot.html#dkrz-python-matplotlib-example-icon-r02b09-datashader-plot
    # to test-implement later
    # this is especially relevant for very large global data plots at high
    # resolution

    # despite the global keyword, this is constraining to the real extend
    # with "plate-caree"; it does not work with always-global "orthographic"
    # these values in case of the pan-EU ICON CORDEX based grids this needs
    # to be reotated coordinates according to the setting from the rotated pole,
    # either rotattion is done here based on the settings of rotated pole
    # or it is done
    plt.xlim(lon_min, lon_max)
    plt.ylim(lat_min, lat_max)
    cb = plt.colorbar(pdo)

    # a last minute rasterization for the visualisation is also possible
    # beforehand; works OK, but only when used as above
    #pdo.set_rasterized(True)

    # issues:
    # - the triangular mesh is only retained when writing the plot to pdf or svg
    # - with png and gridded formats the triangular mesh is just resterized
    # - when using pdf and svg though:
    #   * there are still thin hairlines in the PDF betwen the triangles, even if 
    #     edges are not plotted, this is not ideal
    #   * does not show interactively and 
    #   * does not show in Adobe Acrobat (i.e., barely visible)
    #   * no hairlines in the print-outs though
    # - solution, see 2b: 
    #   set rasterized-keyword only to the data layer plot, output at 
    #   high dpi, pdf format -> no hairline, exact plotting, filesize OK, 
    #   printing OK, works OK with different pdf viewers
    #plt.show()
    # adjust the dpis to truly control the precision of the triangles in
    # the plot in case it is rasterized
    #fig1.savefig('./test_allVecOption.pdf')  # all-vector version
    fig1.savefig('./test.pdf', dpi=3000)  # this seems optimum
    #fig1.savefig('./test.eps', dpi=600)  # very large filesize 
    #fig1.savefig('./test.svg', dpi=600)  # same filesize as pdf
    fig1.savefig('./test.png', dpi=1500)  # also works, larger filesize than pdf
    #fig1.savefig('./test_3a_pc_data.pdf', bbox_inches='tight', pad_inches=0.1, dpi=600)
    #fig1.savefig('./test.png', bbox_inches='tight', pad_inches=0.1, dpi=600)
 
    t2 = time.time()
    print('Wallclock time:  %0.3f seconds' % (t2-t1))


if __name__ == '__main__':

    main()
