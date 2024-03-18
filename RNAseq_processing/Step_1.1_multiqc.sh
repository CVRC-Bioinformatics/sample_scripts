#!/bin/bash -l
 
#SBATCH -J MultiQC
#SBATCH --mem=10000
#SBATCH -n 2
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=1:00:00
#SBATCH --output="MultiQC-raw-data-%A.out"

#date
d1=$(date +%s)

echo $HOSTNAME

module load python/cpu/3.6.5

multiqc /gpfs/data/labfolder/projectfolder/* -n multiqc_report.html

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
