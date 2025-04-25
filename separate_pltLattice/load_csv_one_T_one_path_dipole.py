import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import math
from matplotlib.ticker import FormatStrFormatter  # Import formatter

from pathlib import Path
# This script plots dipoles for one T, one path

if len(sys.argv) != 5:
    print("wrong number of arguments")
    exit(1)
N0=int(sys.argv[1])
N1=int(sys.argv[2])
TStr=sys.argv[3]
init_path=sys.argv[4]
csvDataFolderRoot = f"../dataAll/N0_{N0}_N1_{N1}/T{TStr}/csvOutAll/"
dipole_each_site_dir=csvDataFolderRoot+"/dipole_each_site/"
Path(dipole_each_site_dir).mkdir(exist_ok=True,parents=True)
one_path_folder=csvDataFolderRoot+f"/init_path{init_path}/"
dipole_csv_file_name=one_path_folder+"avg_dipole_combined.csv"
if not os.path.exists(dipole_csv_file_name):
    print(f"avg_dipole_combined.csv does not exist for {TStr}")
    exit(1)
dipole_arr = np.array(pd.read_csv(dipole_csv_file_name, header=None))
# The four rows: [Px, Py, Qx, Qy]
Px = dipole_arr[0, :]
Py = dipole_arr[1, :]

# Compute the mean of the combined arrays
avg_polarization_x = np.mean(Px)  # Mean of Px and Qx combined
avg_polarization_y = np.mean(Py)  # Mean of Py and Qy combined

# Print the results
print(f"Average polarization along x (Px and Qx combined): {avg_polarization_x}")
print(f"Average polarization along y (Py and Qy combined): {avg_polarization_y}")


# Reshape the dipole components
Px_arr = Px.reshape((N0, N1))
Py_arr = Py.reshape((N0, N1))

# Define the lattice constant
a = 2
# Instead of a flat list followed by meshgrid, generate index grids first.
# Let n0 and n1 be the integer indices corresponding to the two lattice directions.
n0 = np.arange(N0)
n1 = np.arange(N1)

# Use index ordering consistent with how the CSV was written
# For a triangular (non-square) lattice:
i_grid, j_grid = np.meshgrid(n0, n1, indexing="ij")

X_O = a * i_grid
Y_O=a*j_grid

mag_A = np.sqrt(Px_arr**2 + Py_arr**2)

mag_min = mag_A.min()
mag_max = mag_A.max()
# Plot using quiver; the 5th argument is the color array.
plt.figure(figsize=(90, 60))
scale=2

# Plot dipoles for sublattice A with a colormap for the magnitude
qA = plt.quiver(
    X_O, Y_O,
    Px_arr, Py_arr,
    mag_A,
    cmap='viridis',
    scale=scale,
    scale_units='xy',
    angles='xy'
)

plt.xlabel("x", fontsize=100)
plt.ylabel("y", fontsize=100)
avg_polarization_x_str=np.round(avg_polarization_x,3)
avg_polarization_y_str=np.round(avg_polarization_y,3)
plt.title(f"Dipole on each site for T = {TStr}, init_path{init_path}, p={avg_polarization_x_str,avg_polarization_y_str}", fontsize=120)
plt.axis("equal")

# Add colorbar from one of the quiver plots and increase number size on the colorbar.
cbar = plt.colorbar(qA)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))

cbar.ax.tick_params(labelsize=120)

plt.savefig(one_path_folder+f"/dipole_each_site_T{TStr}_init_path{init_path}.png")
plt.savefig(dipole_each_site_dir+f"/dipole_each_site_T{TStr}_init_path{init_path}.png")
plt.close()