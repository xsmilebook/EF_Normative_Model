#!/bin/bash

#SBATCH --job-name=oneback_array_analysis 
#SBATCH --output=/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/logs/fig2/oneback_analysis_%A_%a.out 
#SBATCH --error=/ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/logs/fig2/oneback_analysis_%A_%a.err  
#SBATCH --cpus-per-task=56
#SBATCH --mem=60G           
#SBATCH --array=1-5
#SBATCH --partition=q_cn_2


# 1. psy variable list
PSYC_VARIABLES=(
  "SDQ_PB_sum_z"
  "SDQ_H_sum_z"
  "SDQ_CP_sum_z"
  "SDQ_PP_sum_z"
  "SDQ_ES_sum_z"
)

# 2. run R script
R_SCRIPT="/ibmgpfs/cuizaixu_lab/xuhaoshu/EF_Normative_Model/figure2/figure2_oneback.R"

TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

CURRENT_VAR=${PSYC_VARIABLES[$TASK_ID]}

if [ -z "$CURRENT_VAR" ]; then
  echo "error: failed to get $SLURM_ARRAY_TASK_ID."
  exit 1
fi

echo "----------------------------------------------------"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "variable: $CURRENT_VAR"
echo "----------------------------------------------------"

export R_LIBS_USER=/GPFS/cuizaixu_lab_temp/xuhaoshu/R/packages
module load R/4.2.2

Rscript "$R_SCRIPT" "$CURRENT_VAR"

echo "Task of $CURRENT_VAR completed successfully."```

