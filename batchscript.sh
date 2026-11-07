#!/bin/bash
#SBATCH -A f202500010hpcvlabuminhoa
#SBATCH -p normal-arm
#SBATCH -t 10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=out.o%j
#SBATCH --error=err.e%j

ml GCC/15.2.0
make clean
make
export OMP_NUM_THREADS=1
perf record -g ./md < inputdata.txt
perf report > perfil_relatorio.txt
