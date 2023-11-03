#!/bin/bash
#SBATCH --job-name=TR100
#SBATCH --partition=normal-32g
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16      # Use 12 processors
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=axf200007@utdallas.edu

# load environment-modules itself and compiler-specific modules
. /etc/profile.d/modules.sh

module load intel


./FN.exe $SLURM_ARRAY_TASK_ID


# ./FFMQ.exe
# ./FN.exe > c2_Gshock_lp2.out
