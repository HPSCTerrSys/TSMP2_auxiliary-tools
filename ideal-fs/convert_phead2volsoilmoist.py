# Program to convert pressure from ParFlow to volumetric soil moisture and relative saturation

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
import os

# -------------------------
# van Genuchten Parameters from ParFlow namelist
# -------------------------
alpha = 2.69              # 
n = 1.41
m = 1 - 1 / n
porosity = 0.389          # Geom.domain.Porosity.Value
S_res = 0.08              # Geom.domain.Saturation.SRes

theta_s = porosity                    # Saturated water content
theta_r = porosity * S_res            # Residual water content

# -------------------------
# Vertical dzScale (cell thicknesses in cm)
# -------------------------
dz_scale = np.array([
    10, 7.5, 7.5, 5, 5, 5, 2, 2, 2, 2, 1.14, 1.04, 0.94, 0.84, 0.74,
    0.64, 0.54, 0.44, 0.4, 0.36, 0.32, 0.28, 0.24, 0.2, 0.16, 0.12,
    0.08, 0.06, 0.04, 0.02
])

# -------------------------
# van Genuchten Function
# -------------------------
def theta_vg(h, alpha, n, m, theta_s, theta_r):
    h = np.abs(h)
    theta = theta_r + (theta_s - theta_r) / ((1 + (alpha * h)**n)**m)
    return np.clip(theta, theta_r, theta_s)

# -------------------------
# Process Single File
# -------------------------
def process_netcdf(file_path, pressure_head_var='pressure'):
    ds = xr.open_dataset(file_path)

    # Extract pressure head
    h = np.squeeze(ds[pressure_head_var])

    # Compute volumetric water content θ
    theta = theta_vg(h, alpha, n, m, theta_s, theta_r)
    ds['theta'] = (h.dims, theta.data)
    ds['theta'].attrs['units'] = 'cm³/cm³'
    ds['theta'].attrs['long_name'] = 'Volumetric Water Content'

    # Compute Saturation: S = θ / φ
    saturation = theta / porosity
    ds['saturation'] = (h.dims, saturation.data)
    ds['saturation'].attrs['units'] = 'fraction'
    ds['saturation'].attrs['long_name'] = 'Saturation (theta / porosity)'

    # Add dzScale if it matches vertical dimension
    vertical_dim = None
    for dim in h.dims:
        if ds[dim].size == dz_scale.size:
            vertical_dim = dim
            break

    if vertical_dim:
        ds['dz'] = (vertical_dim, dz_scale)
        ds['dz'].attrs['units'] = 'cm'
        ds['dz'].attrs['long_name'] = 'Layer Thickness'

    return ds

# -------------------------
# Batch Process All Files
# -------------------------
if __name__ == "__main__":
    dirname = "/p/scratch/cslts/poll1/sim/ideal/simexp_ideal_fs-eclmcurv_ICON-eCLM-ParFlow/run/sim_pft00-sid02-sv06_eclmparflow_20150701/"
    input_files = sorted(glob.glob(dirname+"fs-idealnwp.out.*.nc"))
    output_dir = "processed_output"
    os.makedirs(output_dir, exist_ok=True)

    for file_path in input_files:
        print(f"Processing: {file_path}")
        try:
            ds_out = process_netcdf(file_path, pressure_head_var='pressure')

            # Create output filename
            base_name = os.path.basename(file_path)
            output_path = os.path.join(output_dir, f"processed_{base_name}")

            # Save processed dataset
            ds_out.to_netcdf(output_path)
            print(f"Saved: {output_path}")
        except Exception as e:
            print(f"Failed to process {file_path}: {e}")

