#!/bin/bash -l

#SBATCH -J merge
#SBATCH --mem=250000
#SBATCH -n 1
#SBATCH --array=1-18
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=600
#SBATCH --output="merge_data-%A_%a.out"

#date
d1=$(date +%s)

##### merge files as a job
path='/gpfs/data/labfolder/projectfolder'

lane1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 1`
lane2=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 2`
merged=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 3`

cat $path/raw_data/$lane1.fastq.gz $path/raw_data/$lane2.fastq.gz > $path/merged_data/$merged.fastq.gz

echo $merged

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)