#!/bin/bash
#SBATCH --job-name=TissueR_ParamSweep
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --partition=ncpu
#SBATCH --mail-user=harry.booth@crick.ac.uk
#SBATCH --mail-type=ALL

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia ParameterSweep.jl