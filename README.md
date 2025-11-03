# Auxiliary tools for TSMP2

This is the repository of auxiliary software tools for the [TSMP2 model system](https://github.com/HPSCTerrSys/TSMP2), complementing other TSMP2 simulation environment components or existing miscellaneous software tools and libraries. The auxiliary tools for TSMP2 serve primarily the purposes of the TSMP2 users and development team within [HPSC TerrSys](https://github.com/HPSCTerrSys).

## [visualization_icon_2d-map_triangular-grid](visualization_icon_2d-map_triangular-grid/)

Python scripts (based on widely-used libraries, such as matplotlib, cartopy, xarray, numpy) to visualize TSMP2 model output (from ICON or eCLM) on its native icosahedral (triangular) grid on a map, including some examples to demonstrate different settings. We consider a special case, where the model is set up in a limited area mode (LAM) over Europe, covering the EURO-CORDEX pan-European model domain, defined by a rotated longitude-latitude grid.

- [map_2d_icon_icosahedral_vis_template_LAM-setup.py](visualization_icon_2d-map_triangular-grid/map_2d_icon_icosahedral_vis_template_LAM-setup.py): Visualisation template for plotting ICON data on its icosahedral native grid on a map using matplotlib and cartopy. This is a template and demonstrator on how TSMP2 ICON (or eCLM) data on its native simulation grid can be plotted with a rotated map projection as a basis, including proper colorbar, etc. Data is plotted without distortion, either vectorized or rasterized, the map edges align with the data. The tool is currently configured to plot 3km or 12km horizontal resolution ICON data (R13B07, R13B05) over a EURO-CORDEX EUR-12 model domain. A recommendation for various use cases is given at the end.

- [map_2d_icon_icosahedral_grid_plotting_testing_and_demo.py](visualization_icon_2d-map_triangular-grid/map_2d_icon_icosahedral_grid_plotting_testing_and_demo.py): Visualisation testing and demo tool on how to plot data on a triangular grid using various plotting and triangulation options; use this as a basis for further visualisation tool developments, when, e.g. ICON data is to be plotted on its native grid. Deliberately lots of hardcoding with many comments. Many explanations added during exploration of the functionalities.

## [visualization_parflow_2d-map_global](visualization_parflow_2d-map_global/)

Python script(s) to visualize global ParFlow IHM simulation results at 1km grid spacing on a map.  

- [map_2D_ParFlow-global-highres.py](visualization_parflow_2d-map_global/map_2D_ParFlow-global-highres.py): Read all necessary data into a dataclass object and provide visualisation for different variables and spatial subsets. Deliberately not a generic script, but tailored to the plots as used in Kollet al al. (under review, JHX).

## [visualization_and_handling_of_TSMP2-sim-outputs](visualization_and_handling_of_TSMP2-sim-outputs/)

These are scriptsi and the accompagnying software environment under `env/` (`source <inifile>` on JSC HPC systems), to work with TSMP2 simulation results, primarily in the context of the TSMP2 introduction by Poll et al. (to be submitted to GMDD).

- [map_2d_trigrid_vis_modularized_LAM.py](real-cordex/map_2d_trigrid_vis_modularized_LAM.py): Python script to visualize TSMP2 model output (from ICON or eCLM) on its native icosahedral (triangular) grid on a map for generic variables as well special cases for terrain, pft and soiltexture.

- [plot_ico_2dfield_fast.ipynb](real-cordex/plot_ico_2dfield_fast.ipynb): Python Juypter notebook to visualize ICON data on its native icosahedral grid. Different approach and projection compared to `map_2d_trigrid_vis_modularized_LAM.py`.

- [InputAnalysis.py](real-cordex/InputAnalysis.py): Python script to visualize CLM input data on a curvilinear grid. 

- [TSMP2_fs-idealnwp_visualization.ipynb](ideal-fs/TSMP2_fs-idealnwp_visualization.ipynb): Python Juypter notebook to investigate fs-idealnwp case.

- [convert_phead2volsoilmoist.py](ideal-fs/convert_phead2volsoilmoist.py): Python script convert pressure from ParFlow to volumetric soil moisture and relative saturation.

