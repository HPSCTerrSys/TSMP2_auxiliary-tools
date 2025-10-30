"""
Visualisation script for plotting topography, PFT, and LAI. 

Usage: python InputAnalysis.py case month
e.g. python InputAnalysis.py cordex 5

Author: Stefan POLL s.poll(at)fz-juelich.de
"""

# need to be loaded in env
# module load Python xarray matplotlib Cartopy GEOS/3.9.1

import os
import sys
import copy # only needed in python3.8.5
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def main():

   # check if case is passed
   if len(sys.argv)==1:
       print('Abort: Please add a valid case: cordex, nrw')
       exit()

   # select case
   if sys.argv[1] == 'cordex':
       filename = os.environ['TRAINHOME']+'/TSMP/tsmp_eur11_eraint_eval_v2/input/clm/surfdata_CLM_cordex0.11_436x424_cosmomask_varpft_noice.nc'
       filename_topo = os.environ['TRAINHOME']+'/TSMP/tsmp_eur11_eraint_eval_v2/input/clm/topodata_CLM_cordex0.11_436x424_cosmomask_varpft_noice.nc'
       gridlim = 15
   elif sys.argv[1] == 'nrw':
       filename = os.environ['TRAINHOME']+'/TSMP/tsmp_nrw/input/clm/surfdata_0300x0300.nc'
       filename_topo = os.environ['TRAINHOME']+'/TSMP/tsmp_nrw/input/clm/topodata_0300x0300.nc'
       gridlim = 1
   else:
       print('Abort: Please choose a valid case: cordex, nrw') 
       exit()

   # open netcdf file
   ds = xr.open_dataset(filename)
   dstopo = xr.open_dataset(filename_topo)

   ###########
   ## PFT
   ###########

   # calculate dominant pft
   pct_pft=ds.PCT_PFT
   pct_pft=pct_pft.assign_coords(lsmpft=pct_pft.lsmpft,lsmlat=pct_pft.lsmlat,lsmlon=pct_pft.lsmlon)
   pct_dompftorg=pct_pft.idxmax(dim="lsmpft",skipna=True)
   pftuniq = np.unique(pct_dompftorg)
   if sys.argv[1] == 'cordex':
      landmask=ds.LANDMASK
   elif sys.argv[1] == 'nrw':
      landmask=ds.LANDFRAC_PFT
   pct_dompft=pct_dompftorg*landmask.where(landmask==1)

   # create colormap PFT
   pftcol = np.array( [[138/256, 102/256, 66/256, 1], # 1
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

   pftcmp = ListedColormap(pftcol)
#   print(pftuniq.astype(int))
   
   # create figure and axes object
   fig = plt.figure(figsize=(12,12))
   # choose map projection
   ax = plt.axes(projection=ccrs.PlateCarree())

   ## create outlines
   land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor='none')

   country_50m = cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m',
                                        edgecolor='face',
                                        facecolor='none')

   ax.add_feature(land_50m, edgecolor='gray', linewidth=0.6, zorder=3)
   ax.add_feature(country_50m, edgecolor='gray', linewidth=0.6, zorder=3)
#   ax.add_feature(cfeature.BORDERS, linewidth=0.6, edgecolor='gray', zorder=2)
#   ax.add_feature(cfeature.LAND, color='lightgray',zorder=0,
#   ax.add_feature(land_50m, color='lightgray',zorder=0,
#       linewidth=0.5, edgecolor='black')


   # creat grid in fig 
   ax.gridlines(draw_labels=True,
                linewidth=0.1,
                color='gray',
                xlocs=range(-180,180,gridlim),
                ylocs=range(-90,90,gridlim),
                zorder=4)

   # add title
   ax.set_title('Plant Functional Type (PFT)', fontsize=12, fontweight='bold')

   # create contour line plot
   levels = np.arange(ds.lsmpft.size+1)-0.5
   norm = BoundaryNorm(levels, ncolors=pftcmp.N, clip=True)
#   cnplot = ax.contourf(ds.LONGXY, ds.LATIXY, pct_dompft,
   cnplot = ax.pcolormesh(ds.LONGXY, ds.LATIXY, pct_dompft, 
                          cmap=pftcmp, zorder=1, norm=norm,
                          shading='auto')

   # add colorbar
   cb = plt.colorbar(cnplot, ticks=np.arange(ds.lsmpft.size+1), orientation='horizontal',
      pad=0.05, shrink=0.7)
#   cb.set_label('%')
   cb.ax.set_xticklabels(["BARE","NETT","NEBT","NDBT","","BET","","BDT","","","BDS","","","GRASS","","CROP","URB"],rotation=45, ha='center')

   # save graphic output to file
   plt.savefig('plot_pft_'+sys.argv[1]+'.pdf',
               bbox_inches='tight')


   ###########
   ### LAI
   ###########

   # month to choose
   if len(sys.argv)<3:
       print('Month is not provided for LAI calculation. Take default (5)')
       imon = 5
   else: 
       imon = int(sys.argv[2])
   

   # pct_dompft
   lai = np.empty((ds.lsmlat.size,ds.lsmlon.size))
   for ilat in ds.lsmlat:
      for ilon in ds.lsmlon:
         lai[ilat,ilon] = ds.MONTHLY_LAI[imon,int(pct_dompftorg[ilat,ilon]),ilat,ilon]
   lai=lai*landmask.where(landmask==1)

   # create figure and axes object
   fig = plt.figure(figsize=(12,12))
   # choose map projection
   ax = plt.axes(projection=ccrs.PlateCarree())
   # add feature
   ax.add_feature(land_50m, edgecolor='gray', linewidth=0.6, zorder=3)
   ax.add_feature(country_50m, edgecolor='gray', linewidth=0.6, zorder=3)
   # creat grid in fig 
   ax.gridlines(draw_labels=True, linewidth=0.1, color='gray', xlocs=range(-180,180,gridlim),
                ylocs=range(-90,90,gridlim), zorder=4) 
   # add title
   ax.set_title('Leaf Area Index', fontsize=12, fontweight='bold')
   # create contour line plot
   cnplot = ax.pcolormesh(ds.LONGXY, ds.LATIXY, lai, cmap='jet',shading='auto',zorder=1)

   # add colorbar
   cb = plt.colorbar(cnplot, orientation='horizontal', pad=0.05, shrink=0.7)

   # save graphic output to file
   plt.savefig('plot_lai_'+sys.argv[1]+'.pdf', bbox_inches='tight')

   ###########
   ### TOPO
   ###########

   topo = np.ma.array(dstopo.TOPO,mask=np.isnan(landmask.where(landmask==1)))
   if '3.8.5' in sys.version:
       tcmp = copy.copy(mpl.cm.get_cmap("terrain"))
   else:
       tcmp = plt.cm.get_cmap("terrain").copy() # python >3.9.5
   tcmp = tcmp(np.linspace(0.25, 1, 200))
   tcmp_cut = ListedColormap(tcmp)
   tcmp_cut.set_bad('blue',1.)

   # create figure and axes object
   fig = plt.figure(figsize=(12,12))
   # choose map projection
   ax = plt.axes(projection=ccrs.PlateCarree())
   # add feature
   ax.add_feature(land_50m, edgecolor='gray', linewidth=0.6, zorder=3)
   ax.add_feature(country_50m, edgecolor='gray', linewidth=0.6, zorder=3)
   # creat grid in fig
   ax.gridlines(draw_labels=True, linewidth=0.1, color='gray', xlocs=range(-180,180,gridlim),
                ylocs=range(-90,90,gridlim), zorder=4)
   # add title
   ax.set_title('Topography', fontsize=12, fontweight='bold')
   # create contour line plot
   cnplot = ax.pcolormesh(ds.LONGXY, ds.LATIXY, topo, cmap=tcmp_cut,shading='auto')


# add colorbar
   cb = plt.colorbar(cnplot, orientation='horizontal', pad=0.05, shrink=0.7)
   cb.set_label('m')

   # save graphic output to file
   plt.savefig('plot_topo_'+sys.argv[1]+'.pdf', bbox_inches='tight')

   

if __name__ == '__main__':
    main()
