#! /bin/bash

for num in {1..1000}
do
sbatch -o /ibmgpfs/cuizaixu_lab/xuhaoshu/datasets/yunfu/logs/fig1/S2_bootstrap/GNGd/GNGd_prime_boot${num}.out Run_step2_2_GNGd.sh $num
done

