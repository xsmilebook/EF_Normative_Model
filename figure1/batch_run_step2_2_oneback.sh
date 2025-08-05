#! /bin/bash

for num in {1..1000}
do
sbatch -o /ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/logs/fig1/S2_bootstrap/Oneback/OnebackAcc_boot${num}.out Run_step2_2_oneback.sh $num
done

