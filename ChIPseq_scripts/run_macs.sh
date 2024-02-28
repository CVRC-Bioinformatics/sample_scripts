#!/bin/bash -l

#SBATCH -J run_macs2
#SBATCH --mem=50000
#SBATCH -n 4
#SBATCH --array=1-15
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=700

#date
d1=$(date +%s)

mytreat=`sed -n "${SLURM_ARRAY_TASK_ID}p" macs_params.txt | cut -f 1`
myName=`sed -n "${SLURM_ARRAY_TASK_ID}p" macs_params.txt | cut -f 2`

echo $HOSTNAME
echo $mytreat
echo $myName

module load macs2/2.1.1

macs2 callpeak -t $mytreat -n $myName --outdir /gpfs/data/fisherlab/Bianca_ATAC/macs/ -g 2.7e9 -f BAMPE

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour =$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
