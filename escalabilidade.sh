#!/bin/bash
#SBATCH -A f202500010hpcvlabuminhoa
#SBATCH -p normal-arm
#SBATCH -t 30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=resultados_tempos.txt
#SBATCH --error=err_escala.e%j

# Preparar ambiente
ml purge
ml GCC/15.2.0
make clean
make

# Testar com diferentes números de threads e medir o tempo (wall-clock time)
for t in 1 2 4 8 16 32 48; do
    export OMP_NUM_THREADS=$t
    echo "=== Teste com $t threads ==="
    /usr/bin/time -p ./md < inputdata.txt
    echo ""
done
