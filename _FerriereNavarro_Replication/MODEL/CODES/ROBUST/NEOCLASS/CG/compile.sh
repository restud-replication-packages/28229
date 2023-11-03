#!/bin/bash
#SBATCH --job-name=TR67
#SBATCH --partition=normal-64g
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1      # Use 12 processors
#SBATCH --mem=4g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=axf200007@utdallas.edu

# load environment-modules itself and compiler-specific modules
. /etc/profile.d/modules.sh

module load intel

rm FN.exe 

ifort -o FN.exe normcdf.f Toolbox.f90 GLOBALS.f90 ModuleSAVE.f90 ModuleINIT.f90 FUNCTIONS.f90 ModuleCOMPUTATIONS.f90 ModuleMPC.f90 ModuleSTEADY.f90 FUNCTIONS_TR.f90 ModuleCOMPUTATIONS_TR.f90 ModuleQERRORS_TR.f90 ModuleTRANSITION.f90 MAIN.f90 -mcmodel=large -shared-intel -qopenmp -mkl


#srun ./FN.exe


# ./FFMQ.exe
# ./FN.exe > c2_Gshock_lp2.out
