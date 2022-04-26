## the script was used to measure the mean depth of the whole genes vs. their CDS (in Hs population)

## OR6 gene
cd /media/megan/Elements3/Rong/WGS
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:10583-21176	 > Hs_tom/depth/OR6_gene.depth
## OR6 CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:21136-21176 > Hs_tom/depth/OR6_CDS_1.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19954-20042 > Hs_tom/depth/OR6_CDS_2.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19633-19802 > Hs_tom/depth/OR6_CDS_3.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19056-19166 > Hs_tom/depth/OR6_CDS_4.depth
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:18496-18728 > Hs_tom/depth/OR6_CDS_5.depth
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:16941-17156 > Hs_tom/depth/OR6_CDS_6.depth
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:16455-16557 > Hs_tom/depth/OR6_CDS_7.depth 
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:16313-16375 > Hs_tom/depth/OR6_CDS_8.depth 
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:16145-16177 > Hs_tom/depth/OR6_CDS_9.depth 
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:15655-15728 > Hs_tom/depth/OR6_CDS_10.depth 
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:15449-15530 > Hs_tom/depth/OR6_CDS_11.depth 
samtools depth -f Hs_bamfiles.txt  -r Scaffold_4031_rc.Scaffold_489:10583-10642 > Hs_tom/depth/OR6_CDS_12.depth 

cd Hs_tom/depth
cat OR6_CDS_*.depth > OR6_CDS.depth

## OR14 gene
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:173309-182119 > Hs_tom/depth/OR14_gene.depth
## OR14 CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:173309-173352 > Hs_tom/depth/OR14_CDS_1.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:173687-173963 > Hs_tom/depth/OR14_CDS_2.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:174468-174578 > Hs_tom/depth/OR14_CDS_3.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:174735-174967 > Hs_tom/depth/OR14_CDS_4.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:175197-175451 > Hs_tom/depth/OR14_CDS_5.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:175758-175860 > Hs_tom/depth/OR14_CDS_6.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:179832-179853 > Hs_tom/depth/OR14_CDS_7.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:179991-180064 > Hs_tom/depth/OR14_CDS_8.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:180218-180373 > Hs_tom/depth/OR14_CDS_9.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:182072-182119 > Hs_tom/depth/OR14_CDS_10.depth

cd Hs_tom/depth/
cat OR14_CDS_*.depth > OR14_CDS.depth

