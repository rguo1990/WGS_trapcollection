## the script was used to measure the mean depth of the whole genes vs. their CDS (in Hs population)

## OR6 gene
cd /media/megan/Elements3/Rong/WGS
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:10583-21176	 > Hs_tom/depth/OR6_gene.depth
## OR6 CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:21136-21176 > Hs_tom/depth/OR6_CDS_1.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19954-20042 > Hs_tom/depth/OR6_CDS_2.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19633-19802 > Hs_tom/depth/OR6_CDS_3.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:19056-19166 > Hs_tom/depth/OR6_CDS_4.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:18496-18728 > Hs_tom/depth/OR6_CDS_5.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:16941-17156 > Hs_tom/depth/OR6_CDS_6.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:16455-16557 > Hs_tom/depth/OR6_CDS_7.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:16313-16375 > Hs_tom/depth/OR6_CDS_8.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:16145-16177 > Hs_tom/depth/OR6_CDS_9.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:15655-15728 > Hs_tom/depth/OR6_CDS_10.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:15449-15530 > Hs_tom/depth/OR6_CDS_11.depth 
samtools depth -f Hs_bamfiles.txt -r Scaffold_4031_rc.Scaffold_489:10583-10642 > Hs_tom/depth/OR6_CDS_12.depth 

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

## OR16 gene
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:199741-206825 > Hs_tom/depth/OR16_gene.depth
## OR16 CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:206203-206243 > Hs_tom/depth/OR16_CDS_1.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:205652-205910 > Hs_tom/depth/OR16_CDS_2.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:205127-205237 > Hs_tom/depth/OR16_CDS_3.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:204102-204334 > Hs_tom/depth/OR16_CDS_4.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:203280-203501 > Hs_tom/depth/OR16_CDS_5.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:202381-202483 > Hs_tom/depth/OR16_CDS_6.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:202210-202305 > Hs_tom/depth/OR16_CDS_7.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:201615-201770 > Hs_tom/depth/OR16_CDS_8.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_656_new:199926-199973 > Hs_tom/depth/OR16_CDS_9.depth

cd Hs_tom/depth/
cat OR16_CDS_*.depth > OR16_CDS.depth

## OR13 gene
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:135542-140490 > Hs_tom/depth/OR13_gene.depth
## OR13 CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:135623-135654 > Hs_tom/depth/OR13_CDS_1.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:135862-136099 > Hs_tom/depth/OR13_CDS_2.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:136920-137027 > Hs_tom/depth/OR13_CDS_3.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:137199-137680 > Hs_tom/depth/OR13_CDS_4.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:137904-138006 > Hs_tom/depth/OR13_CDS_5.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:138158-138253 > Hs_tom/depth/OR13_CDS_6.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:138520-138675 > Hs_tom/depth/OR13_CDS_7.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_2038:140086-140148 > Hs_tom/depth/OR13_CDS_8.depth
cd Hs_tom/depth/
cat OR13_CDS_*.depth > OR13_CDS.depth

## orco gene
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:38994-54701 > Hs_tom/depth/orco_gene.depth
## orco CDS
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:54602-54701 > Hs_tom/depth/orco_CDS_1.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:53477-53687 > Hs_tom/depth/orco_CDS_2.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:49768-49904 > Hs_tom/depth/orco_CDS_3.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:46998-47158 > Hs_tom/depth/orco_CDS_4.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:46320-46506 > Hs_tom/depth/orco_CDS_5.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:45656-45866 > Hs_tom/depth/orco_CDS_6.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:41461-41560 > Hs_tom/depth/orco_CDS_7.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:40055-40159 > Hs_tom/depth/orco_CDS_8.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:39205-39360 > Hs_tom/depth/orco_CDS_9.depth
samtools depth -f Hs_bamfiles.txt -r Scaffold_68:38994-39044 > Hs_tom/depth/orco_CDS_10.depth
cd Hs_tom/depth/
cat orco_CDS_*.depth > orco_CDS.depth

