#!/bin/bash -l

#SBATCH -J run_bwa
#SBATCH --mem=50000
#SBATCH -n 8
#SBATCH --array=1-15
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=700

#date
d1=$(date +%s)

myread1=`sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_files_paired.txt | cut -f 1`
myread2=`sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_files_paired.txt | cut -f 2`
mybam=`sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_files_paired.txt | cut -f 3`

echo $HOSTNAME
echo $myread1
echo $myread2
echo $mybam

module load bwa/0.7.17
module load samtools/1.9

bwa mem -t 8 genomes/mm10/mm10.fa $myread1 $myread2 | samtools view -bS - | samtools sort - > $mybam

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour =$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
