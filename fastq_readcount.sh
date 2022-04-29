## the command was used to count the number of reads of the fastq files

for i in `ls *.fastq.gz`; do echo $(zcat ${i} | wc -l)/4|bc; done

## the command was used to measure GC% of the fastq file
awk '(NR%4==2) {N1+=length($0);gsub(/[AT]/,"");N2+=length($0);}END{print N2/N1;}' in.fastq
