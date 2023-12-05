#!/usr/bin/env python

# SBATCH --exclude=agate-12,agate-13,agate-14,agate-15,agate-16,agate-18,agate-19,agate-28,agate-41,agate-42
# SBATCH --nodes=1 --partition=high
# SBATCH --ntasks-per-node=32
# SBATCH --ntasks-per-core=1
# SBATCH --threads-per-core=1
# SBATCH --output=job.out
# SBATCH --error=job.err
# SBATCH --time=24:00:00
# SBATCH --verbose


from ase.io import write, read
import copy
import shutil
from glob import glob
import numpy as np
import pickle
import math
from ase.db import connect
import glob, os, sys
import subprocess

working_dir = '/home/sausiva/Diffusion_Paper/iter_3/Diffusion_MD_0.1'  # CHANGE ITER

files = [dirs for dirs in os.listdir(working_dir)]
files.sort()

for i, each_file in enumerate(files):
    # smiles = each_file.split("_")[-1]
    # file_number = each_file.split("_")[0]

    folder_path = f'{working_dir}/{each_file}'
    os.chdir(folder_path)
    bashCommand = f"dptools sample -m diffusion_iter3 -n 100 --lo 0.25 --hi 2.00 -o Uncertain_configs.traj {working_dir}/{each_file}/diffusion.traj"  # CHANGE ITER

    # print('b')

    subprocess.run(bashCommand.split())

    # print('a')
    # bashCommand1 = "rm -rf dev.npy"
    # subprocess.run(bashCommand1.split())
    # print('b')