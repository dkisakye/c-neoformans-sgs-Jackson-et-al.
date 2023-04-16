#!/bin/bash
#BATCH -N 1
#SBATCH -n 16
#SBATCH --partition=amd2tb
#SBATCH --mem=40000mb
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=sgs_variant_calling_pipeline.sh
#SBATCH -A nielsenk
#SBATCH -o pipeline%j.output
#SBATCH -e pipeline_%j.error

##### Variables ####

OUTDIR="../results"
INDIR="../samples"
REF="../H99/FungiDB-57_CneoformansH99_Genome.fasta" #reference genome

module load fastqc

mkdir -p $OUTDIR/{fastqc,trimmed_reads,aligned_data,vcfs}

####--Perform quality assessment-- #####

for SAMPLE in $(cat ./sample_list.txt)

do

R1=${INDIR}/${SAMPLE}/*R1*.gz
R2=${INDIR}/${SAMPLE}/*R2*.gz

fastqc $R1 $R2 -o ${OUTDIR}/fastqc

rm -rf ${OUTDIR}/fastqc/*.zip 

####--Trim the fastq reads--####


mkdir -p $OUTDIR/trimmed_reads/${SAMPLE}

trim_galore -j 8 --paired $R1 $R2 -q 28 --length 50 --fastqc -o ${OUTDIR}/trimmed_reads/${SAMPLE}  

###--Perform Alignment--####

module load bwa

#bwa index $REF

mkdir -p $OUTDIR/aligned_data/${SAMPLE}

bwa mem -t 4 -v 1 ${REF} ${OUTDIR}/trimmed_reads/${SAMPLE}/*R1*.gz ${OUTDIR}/trimmed_reads/${SAMPLE}/*R2*.gz > ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}.sam

####--Convert aligned SAM file to BAM--###

samtools view -h -b ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}.sam > ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}.bam

done

