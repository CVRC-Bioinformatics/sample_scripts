#!/bin/bash -l

#SBATCH -J findMotifs
#SBATCH --mem=50000
#SBATCH -n 4
#SBATCH --array=1-4
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=700

#date
d1=$(date +%s)

myBed=`sed -n "${SLURM_ARRAY_TASK_ID}p" homer_params.txt | cut -f 1`
myOutDir=`sed -n "${SLURM_ARRAY_TASK_ID}p" homer_params.txt | cut -f 2`

echo $HOSTNAME
echo $myBed
echo $myOutDir

module load homer/4.10

findMotifsGenome.pl $myBed genomes/mm10/mm10.fa $myOutDir -size 200 -p 4

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
