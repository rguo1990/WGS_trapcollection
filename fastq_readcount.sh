## the command was used to count the number of reads of the fastq files

for i in `ls *.fastq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done
