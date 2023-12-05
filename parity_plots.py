from ase.io import read
from ase.db import connect

from deepmd.infer import DeepPot as DP

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from dptools.utils import colors
import glob

def get_mse(data):
        return np.mean((data[:, 0] - data[:, 1])**2)

def plot_yx(dft, ax):
    xrng = max(dft) - min(dft)
    xmin = min(dft) - 0.5 * xrng
    xmax = max(dft) + 0.5 * xrng
    ax.plot([xmin, xmax], [xmin, xmax], "--k", zorder=1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([xmin, xmax])



#dir = "/home/sausiva/Gas_adsorbates_surface_ML/1_retraining_tests/for_saurabh" #CHANGE DIR

systems = read(f"/global/cfs/cdirs/m3548/tdprice/\
10_ML_training_loop/iter1/rxn_1_iter1.traj",":") #CHANGE

type_map = {"Ag": 0,"C": 1,"H": 2,"O": 3}

dp=DP("/global/cfs/cdirs/m3548/tdprice/10_ML_training_loop/\
iter1/train/00/graph.pb")

energies = []
forces = []

pos = np.array([a.get_positions().flatten() for a in systems])
cell = np.array([a.cell.array.flatten() for a in systems])
types = np.array([type_map[a.symbol] for a in systems[0]])

e_dp, f_dp, v_dp = dp.eval(pos, cell, types)
e_dp = e_dp.flatten()
f_dp = f_dp.flatten()
v_dp = v_dp.flatten()

e_dft = np.array([a.get_potential_energy() for a in systems])
f_dft = np.array([a.get_forces() for a in systems]).flatten()

energies = np.append(e_dft[:, np.newaxis], e_dp[:, np.newaxis], axis=1)
forces = np.append(f_dft[:, np.newaxis], f_dp[:, np.newaxis], axis=1)

width = 0.5 + 4 * 2
fig, axs = plt.subplots(1, 2, figsize=(width, 3.5))

e_data = np.vstack(energies)
f_data = np.vstack(forces)

print(e_data[:, 0])
print(e_data[:, 1])


err1 = get_mse(e_data)
err2 = get_mse(f_data)

print(err1,err2)

axs[0].plot(e_data[:, 0], e_data[:, 1], "ro", ms=3)
plot_yx(e_data[:, 0], axs[0])
axs[0].annotate(f"MSE = {err1:.3e}", xy=(0.1, 0.85),xycoords="axes fraction", fontsize=12)
axs[0].set_ylabel(f"DP Energy (eV)", fontsize=14)
axs[0].set_xlabel(f"DFT Energy (eV)", fontsize=14)
axs[0].tick_params(labelsize=12)

axs[1].plot(f_data[:, 0], f_data[:, 1], "bo", ms=3)
plot_yx(f_data[:, 0], axs[1])
axs[1].annotate(f"MSE = {err2:.3e}", xy=(0.1, 0.85),xycoords="axes fraction", fontsize=12)
axs[1].set_ylabel(f"DP Forces (eV/A)", fontsize=14)
axs[1].set_xlabel(f"DFT Forces (eV/A)", fontsize=14)
axs[1].tick_params(labelsize=12)


plt.tight_layout()
plt.subplots_adjust(wspace=0.5)

plt.savefig('HFSP_TS_parity.png')
plt.close()