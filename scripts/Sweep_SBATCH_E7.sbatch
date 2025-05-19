#!/bin/bash
#SBATCH --job-name=TissueR_ParamSweep_1000_v2
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --partition=ncpu
#SBATCH --mail-user=harry.booth@crick.ac.uk
#SBATCH --mail-type=ALL

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia ParameterSweep_E7.jl