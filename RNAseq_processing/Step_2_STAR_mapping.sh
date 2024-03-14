#!/bin/bash -l
 
#SBATCH -J STAR_map
#SBATCH --mem=250000
#SBATCH -n 8
#SBATCH --array=1-9
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=600
#SBATCH --output="STAR_map-%A_%a.out"
 
#date
d1=$(date +%s)

path='/gpfs/data/labfolder/projectfolder'

read1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 1`
read2=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 2`
output_name=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 3`

echo $HOSTNAME
echo $read1
echo $read2
echo $output_name
 
module load star/2.6.1d
 
STAR --genomeDir /genomes/mm10/STAR --runThreadN 8 --readFilesCommand zcat --readFilesIn $path/raw_data/$read1.fastq.gz $path/raw_data/$read2.fastq.gz --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $path/STAR_mapped_data/$output_name. --quantMode GeneCounts

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
