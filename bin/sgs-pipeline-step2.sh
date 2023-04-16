#!/bin/bash
#BATCH -N 1
#SBATCH -n 16
#SBATCH --partition=amd2tb
#SBATCH --mem=40000mb
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=sgs_variant_calling_pipeline-step2.sh
#SBATCH -A nielsenk
#SBATCH -o pipeline_2%j.output
#SBATCH -e pipeline_2%j.error

#Variables
OUTDIR="../results" 
INDIR="../samples"
REF="../H99/FungiDB-57_CneoformansH99_Genome.fasta" #reference genome

# I didn't assign read groups at mapping so I used other script to assign RGs for the bam files generated. Each pool of samples(1-4,A) has different read group information. Assigning read groups created *RG.bam and *RG.bai files.

#### Sort the bam files ####

module load samtools

for SAMPLE in `cat ./sample_list.txt`
do

samtools sort ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_RG.bam > ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_RG_sorted.bam

#### Mark--duplicates--####
module load picard-tools

java -Xmx4g -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar MarkDuplicates I=${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_RG_sorted.bam O=${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_marked_sorted_RG.bam M=${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_marked_dup_metrics.txt CREATE_INDEX=True


####--Perform variant calling--####

#Index the reference genome using 

samtools faidx $REF

module load gatk/4.1.2

mkdir -p  ${OUTDIR}/vcfs/${SAMPLE}

gatk --java-options "-Xmx4g" HaplotypeCaller -R /home/nielsenk/shared/MinION_project/RIS_analysis_2/H99/fungidb-genomes-2022-05-16/FungiDB-57_CneoformansH99_Genome.fasta -I ${OUTDIR}/aligned_data/${SAMPLE}/${SAMPLE}_marked_sorted_RG.bam  -O ${OUTDIR}/vcfs/${SAMPLE}/${SAMPLE}.vcf.gz

#Filter to retain SNPs

gatk SelectVariants \
  -V ${OUTDIR}/vcfs/${SAMPLE}/${SAMPLE}.vcf.gz \
    -select-type SNP \
    -O ${OUTDIR}/vcfs/${SAMPLE}/${SAMPLE}_snps.vcf.gz 

#Hard filter SNPs
gatk VariantFiltration \
  -V ${OUTDIR}/vcfs/${SAMPLE}/${SAMPLE}_snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${OUTDIR}/vcfs/${SAMPLE}/${SAMPLE}_filtered_snps.vcf.gz

done






