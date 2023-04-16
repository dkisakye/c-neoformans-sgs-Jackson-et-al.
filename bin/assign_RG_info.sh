#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --partition=amd2tb
#SBATCH --mem=40000mb
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=assign_RG1
#SBATCH -o assign_RG1%j.output
#SBATCH -e assign_RG1%j.error


#Assign read group information to reads. RGs are needed for GATK to run. The samples had 5 different read groups (RG1, RG2,RG3, RG4, RGA)

module load picard/2.25.6

for SAMPLE in $(cat ../RG_info/RG1.txt)
	do

java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar AddOrReplaceReadGroups -I ../results/aligned_data/${SAMPLE}/${SAMPLE}.bam -O ../results/aligned_data/${SAMPLE}/${SAMPLE}_RG.bam -SO coordinate -ID M00262:35:000000000 -LB 1 -PU AGB80 -PL Illumina -SM ${SAMPLE} -CREATE_INDEX True

done

#for SAMPLE in `cat ../RG_info/RG2.txt`
#	do

#java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar AddOrReplaceReadGroups -I ../results/aligned_data/${SAMPLE}/${SAMPLE}.bam -O ../results/aligned_data/${SAMPLE}/${SAMPLE}_RG.bam -SO coordinate -ID M00784:236:000000000 -LB 2 -PU AGBCE -PL Illumina -SM ${SAMPLE} -CREATE_INDEX True

	#done

#for SAMPLE in `cat ../RG_info/RG3.txt`
#        do

#java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar AddOrReplaceReadGroups -I ../results/aligned_data/${SAMPLE}/${SAMPLE}.bam -O ../results/aligned_data/${SAMPLE}/${SAMPLE}_RG.bam -SO coordinate -ID M00262:37:000000000 -LB 3 -PU AH5AB -PL Illumina -SM ${SAMPLE} -CREATE_INDEX True
	#done

#for SAMPLE in `cat ../RG_info/RG4.txt`
#	do

#java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar AddOrReplaceReadGroups -I ../aligned_data/${SAMPLE}/${SAMPLE}.bam -O ../aligned_data/${SAMPLE}/${SAMPLE}_RG.bam -SO coordinate -ID M00784:239:000000000 -LB 4 -PU AH2MU -PL Illumina -SM ${SAMPLE} -CREATE_INDEX True

#done

#for SAMPLE in `cat ../RG_info/RGA.txt`
#	do

#java -jar /panfs/roc/msisoft/picard/2.25.6/picard.jar AddOrReplaceReadGroups -I ../aligned_data/${SAMPLE}/${SAMPLE}.bam -O ../aligned_data/${SAMPLE}/${SAMPLE}_RG.bam -SO coordinate -ID R0209721 -LB A -PU 0216 -PL Illumina -SM ${SAMPLE} -CREATE_INDEX True

	#done
