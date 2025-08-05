#!/bin/bash
#SBATCH -J gamlss
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH -p q_cn
export R_LIBS_USER=/GPFS/cuizaixu_lab_temp/xuhaoshu/R/packages
module load R/4.2.2
n=$1
Rscript /ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/figure1/Step2_2_bootstrap_onebackACC.R $n
