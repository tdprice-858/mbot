#!/usr/bin/env bash


source /home/sausiva/.bashrc
cd /home/sausiva/Diffusion_Paper
conda activate dpmd
sbatch -J sampling /home/sausiva/Diffusion_Paper/generate_uncertain_configs.py