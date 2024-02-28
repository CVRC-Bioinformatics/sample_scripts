#!/bin/bash -l
 
#SBATCH -J FeatCount
#SBATCH --mem=50000
#SBATCH -n 12
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=12:00:00
#SBATCH --output="FeatCounts-%A.out"

#date
d1=$(date +%s)
 
echo $HOSTNAME
 
myoutfile='projectfolder/FeatureCounts_data/FeatureCounts_table'
mybams=`cat projectfolder/key_to_FeatureCounts.txt`

module load subread/1.6.3
 
featureCounts -T 12 -p -g gene_id --extraAttributes gene_type,gene_name -M -a genomes/hg38/gencode.v28.annotation.gtf -o $myoutfile $mybams
 
#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)