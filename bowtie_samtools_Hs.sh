#!/bin/bash

## map the reads to H.subflexa reference genome (H.sub_2000.fasta)
## convert sam file to bam file

## Hs_tomatillo
cd /media/megan/Elements3/Rong/WGS/Hs_tom/
for sample in *_R1_trimmed.fastq.gz
do
        echo $sample
        describer=$(echo ${sample} | sed 's/_R1_trimmed.fastq.gz//')
        echo $describer
        bowtie2 -x ~/Rong/H_vir/H_sub -1 ${describer}_R1_trimmed.fastq.gz -2 ${describer}_R2_trimmed.fastq.gz --end-to-end --very-sensitive --no-unal -p 10 > ${describer}.sam

     #convert file from SAM to BAM format
        samtools view -bS ${describer}.sam -o ${describer}.bam

        #Sort BAM file
        samtools sort ${describer}.bam > ${describer}.sorted.bam

        #Index BAM file
        samtools index ${describer}.sorted.bam

	rm ${describer}.sam

done

for sample in *_R1.fastq.gz
do
        echo $sample
        describer=$(echo ${sample} | sed 's/_R1.fastq.gz//')
        echo $describer
        bowtie2 -x ~/Rong/H_vir/H_sub -1 ${describer}_R1.fastq.gz -2 ${describer}_R2.fastq.gz --end-to-end --very-sensitive --no-unal -p 10 > ${describer}.sam

     #convert file from SAM to BAM format
        samtools view -bS ${describer}.sam -o ${describer}.bam

        #Sort BAM file
        samtools sort ${describer}.bam > ${describer}.sorted.bam

        #Index BAM file
        samtools index ${describer}.sorted.bam
	
	rm ${describer}.sam

done


## wild_Hs, repeat
cd /media/megan/Elements3/Rong/WGS/wild_Hs/
for sample in *_R1_trimmed.fastq.gz
do
        echo $sample
        describer=$(echo ${sample} | sed 's/_R1_trimmed.fastq.gz//')
        echo $describer
        bowtie2 -x ~/Rong/H_vir/H_sub -1 ${describer}_R1_trimmed.fastq.gz -2 ${describer}_R2_trimmed.fastq.gz --end-to-end --very-sensitive --no-unal -p 10 > ${describer}.sam

     #convert file from SAM to BAM format
        samtools view -bS ${describer}.sam -o ${describer}.bam

        #Sort BAM file
        samtools sort ${describer}.bam > ${describer}.sorted.bam

        #Index BAM file
        samtools index ${describer}.sorted.bam

	rm ${describer}.sam

done

for sample in *_R1.fastq.gz
do
        echo $sample
        describer=$(echo ${sample} | sed 's/_R1.fastq.gz//')
        echo $describer
        bowtie2 -x ~/Rong/H_vir/H_sub -1 ${describer}_R1.fastq.gz -2 ${describer}_R2.fastq.gz --end-to-end --very-sensitive --no-unal -p 10 > ${describer}.sam

     #convert file from SAM to BAM format
        samtools view -bS ${describer}.sam -o ${describer}.bam

        #Sort BAM file
        samtools sort ${describer}.bam > ${describer}.sorted.bam

        #Index BAM file
        samtools index ${describer}.sorted.bam

	rm ${describer}.sam

done
