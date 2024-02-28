#!/bin/bash -l
 
#SBATCH -J picard
#SBATCH --mem=250000
#SBATCH -n 8
#SBATCH --array=1-9
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=600
#SBATCH --output="picard-%A_%a.out"
 
#date
d1=$(date +%s)

##### remove duplicates from STAR mapped data using picard #### ONLY USE IF NECESSARY ###

path='/gpfs/data/labfolder/projectfolder'

sample=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 3`
input_file=$path/data/STAR_mapped_data/$sample.Aligned.sortedByCoord.out.bam
output_file=$path/data/mapped_data_no_duplicates/$sample.Aligned.sortedByCoord.out.filtered.bam
metrics_file=$path/data/mapped_data_no_duplicates/$sample.dup_metrics.txt

echo $HOSTNAME
echo $sample
 
module load picard-tools/2.18.20
 
java -jar $PICARD_ROOT/libs/picard.jar MarkDuplicates INPUT=$input_file OUTPUT=$output_file METRICS_FILE=$metrics_file REMOVE_DUPLICATES=true 
 
#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
