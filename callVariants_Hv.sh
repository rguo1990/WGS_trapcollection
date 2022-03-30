## the script was used to call variants from the Hv population

#!/bin/bash
## getting a list of Hv bam files
ls Hv_tom/*.sorted.bam > Hvtom_bamfiles.txt
ls Hv_tob/*.sorted.bam > Hvtob_bamfiles.txt
cat Hvtom_bamfiles.txt Hvtob_bamfiles.txt > Hv_bamfiles.txt

## samtools 1.12
## bcftools 1.7

## collect summary information in the bam files, compute the likelihoods of the data and store the likelihoods in a BCF file
nohup bcftools mpileup -b Hv_bamfiles.txt -d 150 -f ~/Rong/H_vir/GCA_002382865.1_K63_refined_pacbio_genomic.fna -o Hvpopulation.bcf -O b
## -b: list of input BAM filenames, -d: max per-file depth; avoids excessive memory usage, -O: output format, b:compressed bcf

## call SNPs/convert to VCF file
## -O v: output uncompressed VCF file, -v: output variant sites only, -m: alternative model for multiallelic and rare-variant calling (conflicts with -c)
bcftools call -vmO v -o Hvpopulation.vcf Hvpopulation.bcf

## filtering by vcftools, version 0.1.15
## --minDP include only genotypes greater than or equal to the "--minDP" value 
## --maf include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value, to exclude sequencing errors
## --max-missing exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed)
vcftools --vcf Hvpopulation.vcf --out Filtered_Hvpopulation --recode --minDP 3 --maf 0.1 --max-missing 0.5
