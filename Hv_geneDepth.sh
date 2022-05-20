## the script was used to measure the mean depth of the whole genes vs. their CDS

## OR6 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:26260-30918 > depth/OR6_gene.depth
## OR6 CDS
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:26260-26300 > depth/OR6_CDS_1.depth 
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:26725-26983 > depth/OR6_CDS_2.depth 
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:27443-27553 > depth/OR6_CDS_3.depth 
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:28260-28492 > depth/OR6_CDS_4.depth
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:28790-29026 > depth/OR6_CDS_5.depth
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:29550-29652 > depth/OR6_CDS_6.depth
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:29719-29814 > depth/OR6_CDS_7.depth 
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:30373-30528 > depth/OR6_CDS_8.depth 
samtools depth -f Hv_bamfiles.txt  -r NWSH01000007.1:30783-30842 > depth/OR6_CDS_9.depth
cd depth/
cat OR6_CDS_*.depth > OR6_CDS.depth

## OR14 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:78506-83663 > depth/OR14_gene.depth
## OR14 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:83472-83515 > depth/OR14_CDS_1.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:82882-83158 > depth/OR14_CDS_2.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:82202-82312 > depth/OR14_CDS_3.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:81610-81842 > depth/OR14_CDS_4.depth 
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:80866-81120 > depth/OR14_CDS_5.depth 
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:80475-80577 > depth/OR14_CDS_6.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:79723-79818 > depth/OR14_CDS_7.depth 
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:79566-79571 > depth/OR14_CDS_8.depth
cd depth/
cat OR14_CDS_*.depth > OR14_CDS.depth

## OR16 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58307-64441 > depth/OR16_gene.depth
## OR16 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58367-58407 > depth/OR16_CDS_1.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:58691-58949 > depth/OR16_CDS_2.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:59928-60038 > depth/OR16_CDS_3.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:61380-61612 > depth/OR16_CDS_4.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62019-62240 > depth/OR16_CDS_5.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62720-62821 > depth/OR16_CDS_6.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000007.1:62894-62894 > depth/OR16_CDS_7.depth
cd depth/
cat OR16_CDS_*.depth > OR16_CDS.depth

## OR13 gene
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:97-4097 > depth/OR13_gene.depth
## OR13 CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:97-709 > depth/OR13_CDS_1.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:2387-2582 > depth/OR13_CDS_2.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:2951-3165 > depth/OR13_CDS_3.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01002174.1:3370-4097 > depth/OR13_CDS_4.depth
cd depth/
cat OR13_CDS_*.depth > OR13_CDS.depth

## orco gene
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:13533-34538 > depth/orco_gene.depth
## orco CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:29259-29358 > depth/orco_CDS_1.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:28112-28322 > depth/orco_CDS_2.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:21276-21412 > depth/orco_CDS_3.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:19216-19376 > depth/orco_CDS_4.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:18769-18955 > depth/orco_CDS_5.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:18130-18340 > depth/orco_CDS_6.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:16132-16231 > depth/orco_CDS_7.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:14799-14903 > depth/orco_CDS_8.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:14241-14396 > depth/orco_CDS_9.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01000897.1:13665-13715 > depth/orco_CDS_10.depth
cd depth/
cat orco_CDS_*.depth > orco_CDS.depth

## OR11 
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:45142-49299 > depth/OR11_gene.depth

## CDS
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:49235-49266 > depth/OR11_CDS_1.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:48784-49021 > depth/OR11_CDS_2.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:48251-48358 > depth/OR11_CDS_3.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:47610-48094 > depth/OR11_CDS_4.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:47256-47358 > depth/OR11_CDS_5.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:46147-46242 > depth/OR11_CDS_6.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:45726-45881 > depth/OR11_CDS_7.depth
samtools depth -f Hv_bamfiles.txt -r NWSH01001866.1:45351-45410 > depth/OR11_CDS_8.depth
cd depth/
cat OR11_CDS_*.depth > OR11_CDS.depth

