## the script was used to measure the mean depth of the whole genes vs. their CDS

## OR6 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:26260-30528 > depth/OR6_gene.depth
## OR6 CDS
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:26260-26300 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:26725-26983 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:27443-27553 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:28260-28492 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:28790-29026 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:29550-29652 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:29719-29814 | samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:30373-30528 > depth/OR6_CDS.depth 

## OR14 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:78506-83663 > depth/OR14_gene.depth
## OR14 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:83472-83515 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:82882-83158 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:82202-82312 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:81610-81842 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:80866-81120 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:80475-80577 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:79723-79818 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:79566-79571 > depth/OR14_CDS.depth

## OR16 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58307-64441 > depth/OR16_gene.depth
## OR16 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58367-58407 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58691-58949 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:59928-60038 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:61380-61612 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62019-62240 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62720-62821 | samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62894-62894 > depth/OR16_CDS.depth

## OR13 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:97-4097 > depth/OR13_gene.depth
## OR13 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:97-709 | samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:2387-2582 | samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:2951-3165 | samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:3370-4097 > depth/OR13_CDS.depth

## orco gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:13533-34538 > depth/orco_gene.depth
## orco CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:29259-29358 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:28112-28322 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:21276-21412 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:19216-19376 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:18769-18955 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:18130-18340 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:16132-16231 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:14799-14903 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:14241-14396 | samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:13665-13715 > > depth/orco_CDS.depth
