#!/usr/bin/env python
# coding: utf-8

'''
Visualisation template for plotting model data on its icosahedral grid on a map

Inputs:
- ICON/eCLM grid file
- External parameter file / Surface file
- Rotated pole longitude and latitude

Outputs:
- 2D plot of eCLM/ICON variables
'''

import time
import xarray as xr
import numpy as np
import getpass
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

__authors__ = "Stefan POLL"
__references__ = [
    "https://github.com/HPSCTerrSys/TSMP2_auxiliary-tools",
    "https://docs.dkrz.de/doc/visualization/sw/python/source_code/python-matplotlib-example-unstructured-icon-triangles-plot-python-3.html",
    "https://docs.dkrz.de/doc/visualization/sw/python/source_code/python-matplotlib-triangular-grid-with-tripcolor-ICON.html",
    ]

def load_datasets(dsPnFn: str, gridPnFn: str):
    ds = xr.open_dataset(dsPnFn)
    dsGrid = xr.open_dataset(gridPnFn)
    return ds, dsGrid


def extract_coordinates(dsGrid, model):
    if (model=='eclm'):
        clon = np.squeeze(dsGrid.xc.values) # ncells
        clat = np.squeeze(dsGrid.yc.values) # ncells
        vlon = np.squeeze(dsGrid.xv.values) # ncells, 3
        vlat = np.squeeze(dsGrid.yv.values) # ncells, 3
    else:
        clon = np.rad2deg(dsGrid.clon.values) # ncells
        clat = np.rad2deg(dsGrid.clat.values) # ncells
        vlon = np.rad2deg(dsGrid.vlon.values) # nvortex
        vlat = np.rad2deg(dsGrid.vlat.values) # nvortex
    return clon, clat, vlon, vlat

def create_mask(clon, clat, model):

    if (model=='icon') or (model=='eclm2'):
        mask = (
            (clat > 0)
            & (clat < 90)
            & (clon > -180)
            & (clon < 180)
        )
    
    if (model=='eclm'):
        mask = (
            (clat > 0)
            & (clat < 90)
            & (clon > 0)
            & (clon < 360)
        )

    return mask


def make_triangulation(vlon, vlat, mask=None):
    """
    Convert (ncells, 3) lon/lat arrays into a matplotlib Triangulation.
    """
    if mask is not None:
        vlon = vlon[mask]
        vlat = vlat[mask]

#    points = np.column_stack([vlon.ravel(), vlat.ravel()])
#
#    coords, inverse_idx = np.unique(points, axis=0, return_inverse=True)
#    x = coords[:, 0]
#    y = coords[:, 1]
#
#    triangles = inverse_idx.reshape(vlon.shape[0], vlon.shape[1])  # (ncells, 3)
#    triang = tri.Triangulation(x, y, triangles)

#    triangles = np.arange(vlon.size).reshape(-1, 3)
#    points = np.column_stack([vlon.ravel(), vlat.ravel()])
#    triang = tri.Triangulation(points[:,0], points[:,1], triangles)

    points = np.column_stack([vlon.ravel(), vlat.ravel()])
    coords, inverse_idx = np.unique(points, axis=0, return_inverse=True)
    triangles = inverse_idx.reshape(vlon.shape[0], 3)
    triang = tri.Triangulation(coords[:, 0], coords[:, 1], triangles)

#    ncells = vlon.shape[0]
#    triang = np.zeros((ncells,vlon.shape[1] , 2), np.float32)
#    for i in range(0, ncells, 1):
#        triang[i,:,0] = np.array(vlon[i,:])
#        triang[i,:,1] = np.array(vlat[i,:])

    return triang


def get_triangulation(dsGrid, mask, vlon, vlat, model):
    if (model=='eclm'):
        # check mask for eCLM
        triang = make_triangulation(vlon, vlat, mask)
        used_vertices = np.arange(vlon.shape[0])

#        print("Shape of triang 2:", triang.shape)
        print("Shape of triangles array:", triang.triangles.shape)
#        print("Shape of vlon:", vlon.shape)

    else:
        voc = dsGrid.vertex_of_cell.T[mask].values - 1
        triang = tri.Triangulation(vlon, vlat, voc)
        used_vertices = np.unique(voc)
    return triang, used_vertices


def rotate_coordinates(vlon, vlat, lon_northpole: float, lat_northpole: float):
    lon_northpole_rot = np.radians(lon_northpole)
    lat_northpole_rot = np.radians(lat_northpole)

    vlon_rad = np.deg2rad(vlon)
    vlat_rad = np.deg2rad(vlat)

    vlon_rot = np.degrees(np.arctan2(
        -np.cos(vlat_rad) * np.sin(vlon_rad - lon_northpole_rot),
        -np.cos(vlat_rad) * np.sin(lat_northpole_rot) * np.cos(vlon_rad - lon_northpole_rot)
        + np.sin(vlat_rad) * np.cos(lat_northpole_rot)
    ))

    vlon_rot[vlon_rot < -180] += 360
    vlon_rot[vlon_rot > 180] -= 360

    vlat_rot = np.degrees(np.arcsin(
        np.sin(vlat_rad) * np.sin(lat_northpole_rot)
        + np.cos(vlat_rad) * np.cos(lat_northpole_rot) * np.cos(vlon_rad - lon_northpole_rot)
    ))
    return vlon_rot, vlat_rot


def get_rotated_bounds(vlon_rot, vlat_rot, used_vertices):
    vlat_rot_min = vlat_rot[used_vertices].min()
    vlat_rot_max = vlat_rot[used_vertices].max()
    vlon_rot_min = vlon_rot[used_vertices].min()
    vlon_rot_max = vlon_rot[used_vertices].max()
    return vlon_rot_min, vlon_rot_max, vlat_rot_min, vlat_rot_max


def select_plotting_var(ds, variable, mask):

    if ( variable == 'terrain' ):
        var = np.ma.masked_where(
               ds['FR_LAND'][:].isel(cell=mask).values <= 0.5,
               ds['topography_c'][:].isel(cell=mask).values
              )
        cmap_used = build_colormap_terrain()

    elif ( variable == 'natpft' ):
        pct_dompft=ds['PCT_NAT_PFT'][:].isel(gridcell=mask).idxmax(dim="natpft",skipna=True).values
#        var = np.ma.masked_where(
#               ds['LANDFRAC_PFT'][:].isel(gridcell=mask).values <= 0.5,
#               pct_dompft
#              )
        var = pct_dompft
        cmap_used = build_colormap_pftnat()

    elif ( variable == 'pft' ):

        pct_nat_pft_cor = ds["PCT_NATVEG"] * ds["PCT_NAT_PFT"] / 100
        pct_cft_cor = ds["PCT_CROP"] * ds["PCT_CFT"] / 100
        
        pct_urban   = ds["PCT_URBAN"].sum(dim="numurbl")
        pct_wetland = ds["PCT_WETLAND"]
        pct_lake    = (ds["PCT_LAKE"])
        pct_glacier = (ds["PCT_GLACIER"])

        nat_arr   = pct_nat_pft_cor.to_numpy()   # shape: (gridcell, new_class_nat)
        cft_arr   = pct_cft_cor.to_numpy()   # shape: (gridcell, new_class_cft)
        urban_arr = pct_urban.to_numpy()[:, np.newaxis] # shape: (gridcell, 1)
        wetland_arr = pct_wetland.to_numpy()[:, np.newaxis]
        lake_arr = pct_lake.to_numpy()[:, np.newaxis]
        glacier_arr = pct_glacier.to_numpy()[:, np.newaxis]

        combined_array = np.concatenate([nat_arr, cft_arr, urban_arr,wetland_arr,lake_arr,glacier_arr], axis=1)

        print(combined_array.shape)
        
        class_labels = [ "BARE","NETTe","NETBo","NDTBo","BETTr","BETTe","BDTTr","BDTTe","BDTBo",
             "BESTe","BDSTe","BDSBo","AC3Gr","C3Gr","C4Gr","Crp","iCrp","URB","WLD","LKE","GLC"]

        dominant_class = np.argmax(combined_array, axis=1)
        # Convert to xarray.DataArray
        dom_pft = xr.DataArray(
            dominant_class,
            dims=("gridcell"),
            coords={
                "gridcell": np.arange(dominant_class.shape[0]),
            },
            name="dominant_landcover"
        )
        var=dom_pft
        cmap_used = build_colormap_pft()
        
        # check value freq
        values, counts = np.unique(dom_pft.values, return_counts=True)
        freq_dict = {class_labels[v]: c for v, c in zip(values, counts)}
        print(freq_dict)

    elif ( variable == 'soiltexture' ):
       
        pct_clay = ds["PCT_CLAY"]
        pct_sand = ds["PCT_SAND"]
#        frland = ds["PFTDATA_MASK"]
        frland = 1.-ds["PCT_WETLAND"]/100

        lev_soil = 2
        pct_clay = pct_clay.values[lev_soil,:]
        pct_sand = pct_sand.values[lev_soil,:]
        pct_silt = 100.-pct_clay-pct_sand

        #
        usda_scs = np.full(frland.shape, np.nan)
        usda_scs = np.where(( frland >= 0.5) & ( (pct_silt + (1.5*pct_clay)) < 15. ),1, usda_scs)                                                              #  1 - sand
        usda_scs = np.where(( frland >= 0.5) & ( (pct_silt + (1.5*pct_clay)) >= 15. ) & ( (pct_silt + (2.*pct_clay)) < 30. ), 2, usda_scs)                     #  2 - loamy sand
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 7. ) & ( pct_clay < 20. ) & ( pct_sand > 52. ) & ( (pct_silt + (2.*pct_clay)) >= 30. ), 3, usda_scs) #  3 - sandy loam (1)
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay < 7. ) & ( pct_silt < 50. ) & ( (pct_silt + (2.*pct_clay)) >= 30. ), 3, usda_scs)                    #  3 - sandy loam (2)
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 7. ) & ( pct_clay < 27. ) & ( pct_silt >= 28. ) & ( pct_silt < 50. ) & ( pct_sand <= 52. ),4, usda_scs) #  4 - loam
        usda_scs = np.where(( frland >= 0.5) & ( pct_silt >= 50. ) & ( pct_clay >= 12. ) & ( pct_clay < 27. ), 5, usda_scs)                                    #  5 - silt loam (1)
        usda_scs = np.where(( frland >= 0.5) & ( pct_silt >= 50. ) & ( pct_silt < 80. ) & ( pct_clay < 12. ), 5, usda_scs)                                     #  5 - silt loam (2)
        usda_scs = np.where(( frland >= 0.5) & ( pct_silt >= 80. ) & ( pct_clay < 12. ), 6, usda_scs)                                                          #  6 - silt
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 20. ) & ( pct_clay < 35. ) & ( pct_silt < 28. ) & ( pct_sand > 45. ), 7, usda_scs)                #  7 - sandy clay loam
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 27. ) & ( pct_clay < 40. ) & ( pct_sand > 20. ) & ( pct_sand <= 45. ), 8, usda_scs)               #  8 - clay loam
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 27. ) & ( pct_clay < 40. ) & ( pct_sand <= 20. ), 9, usda_scs)                                    #  9 - silty clay loam
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 35. ) & ( pct_sand > 45. ), 10, usda_scs)                                                         # 10 - sandy clay
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 40. ) & ( pct_silt >= 40. ), 11, usda_scs)                                                        # 11 - silty clay
        usda_scs = np.where(( frland >= 0.5) & ( pct_clay >= 40. ) & ( pct_sand <= 45. ) & ( pct_silt < 40. ), 12, usda_scs)                                   # 12 - clay

        var = xr.DataArray(
            usda_scs,
            dims=("gridcell"),
            coords={
                "gridcell": np.arange(usda_scs.shape[0]),
            },
            name="soil_texture"
        )
        cmap_used = build_colormap_soiltexture()

        # check value freq
        #class_labels = ["sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt", "sandy clay loam", "clay loam", "silty clay loam", "sandy clay", "silty clay", "clay"]
        values, counts = np.unique(var.values, return_counts=True)
        freq_dict = {f"{v:.2f}": c for v, c in zip(values, counts)}
        print(freq_dict)

    else:
        var = ds[variable]

        ind_z = 0
        ind_t = 0
        # Reduce variable to first dimension
        dims = var.ndim
        if dims == 3:
            var = var.isel({var.dims[0]: ind_t, var.dims[1]: ind_z})
        elif dims == 2:
            var = var.isel({var.dims[0]: ind_t})
        else:
            var = var
        cmap_used = plt.get_cmap("turbo")

    return var, cmap_used

def build_colormap_pft():

    # create colormap PFT
    pftcol = [  (138/256, 102/256,  66/256, 1),  # 1 BARE
                ( 20/256, 100/256,  40/256, 1),  # 2 NETTe
                (169/256, 169/256, 169/256, 1),  # 3 NETBo
                (110/256, 139/256,  61/256, 1),  # 4 NDTBo
                (169/256, 169/256, 169/256, 1),  # 5 BETTr
                (169/256, 169/256, 169/256, 1),  # 6 BETTe
                (169/256, 169/256, 169/256, 1),  # 7 BDTTr
                ( 80/256, 159/256, 101/256, 1),  # 8 BDTTe
#                ( 40/256, 120/256,  60/256, 1),  # 9 BDTBo
                (169/256, 169/256, 169/256, 1),  # 9 BDTBo
                (188/256, 238/256, 104/256, 1),  #10 BESTe
#                (168/256, 218/256,  84/256, 1),  #11 BDSTe
                (169/256, 169/256, 169/256, 1),  #11 BDSTe
#                (148/256, 188/256,  64/256, 1),  #12 BDSBo
                (169/256, 169/256, 169/256, 1),  #12 BDSBo
#                ( 69/256, 229/256,  69/256, 1),  #13 AC3Gr
                (169/256, 169/256, 169/256, 1),  #13 AC3Gr
                (169/256, 169/256, 169/256, 1),  #14 C3Gr
                ( 20/256, 220/256,  20/256, 1),  #15 C4Gr
                (238/256, 216/256, 174/256, 1),  #16 Crp
                (218/256, 196/256, 154/256, 1),  #17 iCrp
                (205/256,   0/256,   0/256, 1),  #18 URB
                (  0/256,   0/256, 205/256, 1),  #19 WLD
                ( 40/256,  40/256, 245/256, 1),  #20 LKE
                (135/256, 206/256, 235/256, 1)   #21 GLC
            ]

    return ListedColormap(pftcol, name="pftcmp")


def build_colormap_pftnat():

    # create colormap PFT
#    pftcol = np.array( [[0/256, 102/256,  204/256, 1], # 0
#                       [138/256, 102/256,  66/256, 1], # 1
    pftcol = np.array([[138/256, 102/256, 66/256,  1], # 1
                       [110/256, 139/256,  61/256, 1], # 2
                       [188/256, 238/256, 104/256, 1], # 3
                       [  0/256, 205/256,   0/256, 1], # 4
                       [169/256, 169/256, 169/256, 1], # 5
                       [ 91/256,  79/256,  61/256, 1], # 6
                       [169/256, 169/256, 169/256, 1], # 7
                       [ 20/256, 100/256,  40/256, 1], # 8
                       [169/256, 169/256, 169/256, 1], # 9
                       [169/256, 169/256, 169/256, 1], #10
                       [  0/256, 255/256,   0/256, 1], #11
                       [169/256, 169/256, 169/256, 1], #12
                       [169/256, 169/256, 169/256, 1], #13
                       [ 20/256, 220/256,  20/256, 1], #14
                       [169/256, 169/256, 169/256, 1], #15
                       [238/256, 216/256, 174/256, 1], #16
                       [205/256,   0/256,   0/256, 1]])#17

    return ListedColormap(pftcol, name="pftcmp")

def build_colormap_soiltexture():

    # create colormap soiltexture
    sltypcol = np.array([   [1.0, 0.92, 0.64, 1],   # 1 Sand        - pale yellow
                            [0.96, 0.85, 0.50, 1],  # 2 Loamy Sand  - yellowish
                            [0.91, 0.80, 0.55, 1],  # 3 Sandy Loam  - light yellow-brown
                            [0.80, 0.60, 0.40, 1],  # 4 Loam        - medium brown
                            [0.78, 0.75, 0.70, 1],  # 5 Silt Loam   - grayish
                            [0.82, 0.80, 0.75, 1],  # 6 Silt        - light gray
                            [0.85, 0.55, 0.40, 1],  # 7 Sandy Clay Loam - reddish brown
                            [0.80, 0.45, 0.30, 1],  # 8 Clay Loam       - medium reddish
                            [0.70, 0.60, 0.55, 1],  # 9 Silty Clay Loam - grayish brown
                            [0.75, 0.45, 0.25, 1],  #10 Sandy Clay      - reddish
                            [0.60, 0.45, 0.40, 1],  #11 Silty Clay      - dark gray-brown
                            [0.55, 0.35, 0.25, 1]])   #12 Clay            - dark reddish brown

    return ListedColormap(sltypcol, name="sltypcmp")


def build_colormap_terrain(base_cmap="terrain", n_colors=50):
    cmap = plt.get_cmap(base_cmap)
    colors = cmap(np.linspace(0, 1, 500))

    # Detect and filter blue-ish colors
    blue_mask = (colors[:, 2] > colors[:, 0]) & (colors[:, 2] > colors[:, 1])
    non_blue_colors = colors[~blue_mask]

    # Resample to n_colors
    indices = np.linspace(0, len(non_blue_colors) - 1, n_colors).astype(int)
    colors_resampled = non_blue_colors[indices]

    return ListedColormap(colors_resampled, name=f"{base_cmap}_no_blue")


def plot_map(var, triang, mask, vlon_rot_min, vlon_rot_max, vlat_rot_min, vlat_rot_max, cmap, plotFn, model, variable):
    crs_data = ccrs.PlateCarree()
    crs_map = ccrs.RotatedPole(pole_longitude=-162, pole_latitude=39.25)

    fig1 = plt.figure(figsize=(5.0, 5.0))
    ax1 = plt.subplot(111, projection=crs_map)
    ax1.set_aspect('equal')

    ax1.coastlines(resolution='50m', linewidth=0.5)
    ax1.add_feature(cfeature.OCEAN, color='azure')
#    ax1.set_title('ICON external parameters (EUR-12),\nplotting demo with icosahedral grid', fontsize=9)
    
    if (variable=='terrain'):
        levelsVals = (np.arange(51) * 50)
    elif (variable=='pftnat'):
        levelsVals = (np.arange(17+1))
    elif (variable=='pft'):
        levelsVals = (np.arange(21+1))
    elif (variable=='soiltexture'):
        levelsVals = np.arange(13)+1
    else:
        levelsVals = np.linspace(var.min().item(), var.max().item(), 50)

#    print("Shape of triangles array:", triang.triangles.shape)
#    print("Number of points in var:", len(var))
   
    pdo = ax1.tripcolor(
        triang,
        facecolors=var,
        transform=crs_data,
        shading='flat',
#        shading='gouraud',
        edgecolors='none',
        rasterized=False,
        cmap=cmap,
        vmin=levelsVals[0],
        vmax=levelsVals[-1]
        )

    plt.xlim(vlon_rot_min, vlon_rot_max)
    plt.ylim(vlat_rot_min, vlat_rot_max)

    if (variable=='terrain'):
        cb = plt.colorbar(pdo, ax=ax1, extend='both', pad=0.03, shrink=0.8, orientation='horizontal', ticks=levelsVals[::10])
        cb.ax.tick_params(labelsize=8)
        cb.set_label('Surface altitude [m]', fontsize=9)
    elif (variable=='pftnat'):
        tick_labels = ["BARE","NETTe","NETBo","NDTBo","BETTr","BETTe","BDTTr","BDTTe","BDTBo","BESTe","BDSTe","BDSBo","AC3Gr","C3Gr","C4Gr","Crp","URB"]
        tick_positions = (levelsVals[:-1] + levelsVals[1:]) / 2
        cb = plt.colorbar(pdo, ax=ax1, extend='neither', pad=0.03, shrink=0.8, orientation='horizontal', ticks=tick_positions)
        cb.ax.tick_params(labelsize=8)
        cb.ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    elif (variable=='pft'):
        tick_labels = [ "BARE","NETTe","NETBo","NDTBo","BETTr","BETTe","BDTTr","BDTTe","BDTBo",
             "BESTe","BDSTe","BDSBo","AC3Gr","C3Gr","C4Gr","Crp","iCrp","URB","WLD","LKE","GLC"]
        tick_positions = (levelsVals[:-1] + levelsVals[1:]) / 2
        cb = plt.colorbar(pdo, ax=ax1, extend='neither', pad=0.03, shrink=0.8, orientation='horizontal', ticks=tick_positions)
        cb.ax.tick_params(labelsize=8)
        cb.ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    elif (variable=='soiltexture'):
        tick_labels = ["sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt", "sandy clay loam", "clay loam", "silty clay loam", "sandy clay", "silty clay", "clay"]
        tick_positions = (levelsVals[:-1] + levelsVals[1:]) / 2
        cb = plt.colorbar(pdo, ax=ax1, extend='neither', pad=0.03, shrink=0.8, orientation='horizontal', ticks=tick_positions)
        cb.ax.tick_params(labelsize=8)
        cb.ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    else:
        cb = plt.colorbar(pdo, ax=ax1, extend='both', pad=0.03, shrink=0.8, orientation='horizontal', ticks=levelsVals[::10])
        cb.ax.tick_params(labelsize=8)
        labname = var.attrs.get("long_name", "")
        cb.set_label(labname, fontsize=9)

    fig1.savefig(plotFn, bbox_inches='tight', pad_inches=0.1, dpi=2000)
#    fig1.savefig('./map_output_no_blue.png', bbox_inches='tight', pad_inches=0.1, dpi=2000)

    plt.show()


def main():
    t1 = time.time()

    lon_northpole = -162.0
    lat_northpole = 39.25

#    variable = 'terrain'
    variable = 'pft'
    username=getpass.getuser()
    dirname='/p/project1/training2538/'+username+'/simexp_real_CORDEX-EUR-11u_icon-eclm-parflow/dta/geo/'

    if (variable == 'natpft') or (variable == 'soiltexture') or (variable == 'pft'):
       model = 'eclm2'
       dsPnFn = dirname+'/eclm/static/surfdata_ICON-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230302_gcvurb-pfsoil_halo.nc'
#       gridPnFn = dirname+'/eclm/static/domain.lnd.ICON-11_ICON-11.230302_landlake_halo.nc'
       gridPnFn = dirname+'/eclm/static/EUR-R13B05_189976_grid.nc'
#       gridPnFn = '/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/gen_domain_files/domain.lnd.EUR-R13B05_EUR-R13B05.251022.nc'
#       gridPnFn = '/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/mkmapgrids/EUR-R13B05_189976_grid.nc'
    elif (variable == 'terrain'):
       model = 'icon'
       dsPnFn = dirname+'/icon/static/external_parameter_icon_europe011_DOM01_tiles.nc'
       gridPnFn = dirname+'/icon/static/europe011_DOM01.nc'
    else:
       model = 'eclm2'
       dsPnFn = dirname+'/eclm/static/surfdata_ICON-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230302_gcvurb-pfsoil_halo.nc'
       gridPnFn = '/p/project1/cslts/poll1/eclm_coupling/CTSM/eCLM_static-file-generator_regen/mkmapgrids/EUR-R13B05_189976_grid.nc'

    plotFn = f'./map_{model}_{variable}_EUR-12.pdf'

    ds, dsGrid = load_datasets(dsPnFn, gridPnFn)
    clon, clat, vlon, vlat = extract_coordinates(dsGrid, model)
    mask = create_mask(clon, clat, model)

    triang, used_vertices = get_triangulation(dsGrid, mask, vlon, vlat, model)
    vlon_rot, vlat_rot = rotate_coordinates(vlon, vlat, lon_northpole, lat_northpole)
    vlon_rot_min, vlon_rot_max, vlat_rot_min, vlat_rot_max = get_rotated_bounds(vlon_rot, vlat_rot, used_vertices)

    var, cmap_used = select_plotting_var(ds, variable, mask)

    plot_map(var, triang, mask, vlon_rot_min, vlon_rot_max, vlat_rot_min, vlat_rot_max, cmap_used, plotFn, model,variable)

    print('exec wallclock time =  %0.3f s' % (time.time() - t1))


if __name__ == '__main__':
    main()

