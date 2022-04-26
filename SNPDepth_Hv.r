## the script was used to plot the read alignment depth of SNPs in the coding regions

library(ggplot2)
library(reshape2)
## OR6_CDS.depth: site depth of coding regions in OR6, measured by samtools depth
## OR6.recode.vcf: vcf file of OR6 gene

or6_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR6_CDS.depth",sep = "\t",header = FALSE)
or6 <- read.loci("~/Downloads/OR6.recode.vcf",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or6_depth) <- names
nrow(or6_depth)
## 1236

## only for SNPs in the coding regions
or6_depth_snp <- or6_depth[or6_depth$pos %in% or6$V2,] 
## nrow(or6_depth_snp) = 32
a=as.numeric(or6_depth_snp[1,3:22])
mean(a)

x <- rep(or6_depth_snp$pos,20)
y <- c(or6_depth_snp$tom10,or6_depth_snp$tom1,or6_depth_snp$tom2,or6_depth_snp$tom3,or6_depth_snp$tom4,
       or6_depth_snp$tom5,or6_depth_snp$tom6,or6_depth_snp$tom7,or6_depth_snp$tom8,or6_depth_snp$tom9,
       or6_depth_snp$tob10,or6_depth_snp$tob1,or6_depth_snp$tob2,or6_depth_snp$tob3,or6_depth_snp$tob4,
       or6_depth_snp$tob5,or6_depth_snp$tob6,or6_depth_snp$tob7,or6_depth_snp$tob8,or6_depth_snp$tob9)
plants <- c(rep("tom10",32),rep("tom1",32),rep("tom2",32),rep("tom3",32),rep("tom4",32),rep("tom5",32),rep("tom6",32),rep("tom7",32),
            rep("tom8",32),rep("tom9",32),rep("tob10",32),rep("tob1",32),rep("tob2",32),rep("tob3",32),rep("tob4",32),
            rep("tob5",32),rep("tob6",32),rep("tob7",32),rep("tob8",32),rep("tob9",32))
df6 <- data.frame(cbind(x,y,plants))
colnames(df6) <- c("pos","depth","plants")

png("OR6_CDS_SNP.png",units="in", width=13, height=6,res=300)
## make boxplot for each SNP site
ggplot(df6, aes(x=factor(x), y=y,group=x)) + geom_boxplot()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="SNPs in OR6 coding regions")
dev.off()

or14_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR14_CDS.depth",sep = "\t", header=FALSE)
or14 <- read.loci("~/Downloads/OR14.recode.vcf",header = FALSE)
colnames(or14_depth) <- names
nrow(or14_depth)
## 1125
or14_depth_snp <- or14_depth[or14_depth$pos %in% or14$V2,] 
## nrow(or14_depth_snp) = 45

x <- rep(or14_depth_snp$pos,20)
y <- c(or14_depth_snp$tom10,or14_depth_snp$tom1,or14_depth_snp$tom2,or14_depth_snp$tom3,or14_depth_snp$tom4,
       or14_depth_snp$tom5,or14_depth_snp$tom6,or14_depth_snp$tom7,or14_depth_snp$tom8,or14_depth_snp$tom9,
       or14_depth_snp$tob10,or14_depth_snp$tob1,or14_depth_snp$tob2,or14_depth_snp$tob3,or14_depth_snp$tob4,
       or14_depth_snp$tob5,or14_depth_snp$tob6,or14_depth_snp$tob7,or14_depth_snp$tob8,or14_depth_snp$tob9)
plants <- c(rep("tom10",45),rep("tom1",45),rep("tom2",45),rep("tom3",45),rep("tom4",45),rep("tom5",45),rep("tom6",45),rep("tom7",45),
            rep("tom8",45),rep("tom9",45),rep("tob10",45),rep("tob1",45),rep("tob2",45),rep("tob3",45),rep("tob4",45),
            rep("tob5",45),rep("tob6",45),rep("tob7",45),rep("tob8",45),rep("tob9",45))
df14 <- data.frame(cbind(x,y,plants))
colnames(df14) <- c("pos","depth","plants")

png("OR14_CDS_SNP.png",units="in", width=17, height=6,res=300)
ggplot(df14, aes(x=factor(x), y=y,group=x)) + geom_boxplot()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="SNPs in OR14 coding regions")
dev.off()

or16_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR16_CDS.depth",sep = "\t", header=FALSE)
or16 <- read.loci("~/Downloads/OR16.recode.vcf",header = FALSE)
colnames(or16_depth) <- names
nrow(or16_depth)
## 969
or16_depth_snp <- or16_depth[or16_depth$pos %in% or16$V2,] 
## nrow(or16_depth_snp) = 34

x <- rep(or16_depth_snp$pos,20)
y <- c(or16_depth_snp$tom10,or16_depth_snp$tom1,or16_depth_snp$tom2,or16_depth_snp$tom3,or16_depth_snp$tom4,
       or16_depth_snp$tom5,or16_depth_snp$tom6,or16_depth_snp$tom7,or16_depth_snp$tom8,or16_depth_snp$tom9,
       or16_depth_snp$tob10,or16_depth_snp$tob1,or16_depth_snp$tob2,or16_depth_snp$tob3,or16_depth_snp$tob4,
       or16_depth_snp$tob5,or16_depth_snp$tob6,or16_depth_snp$tob7,or16_depth_snp$tob8,or16_depth_snp$tob9)
plants <- c(rep("tom10",34),rep("tom1",34),rep("tom2",34),rep("tom3",34),rep("tom4",34),rep("tom5",34),rep("tom6",34),rep("tom7",34),
            rep("tom8",34),rep("tom9",34),rep("tob10",34),rep("tob1",34),rep("tob2",34),rep("tob3",34),rep("tob4",34),
            rep("tob5",34),rep("tob6",34),rep("tob7",34),rep("tob8",34),rep("tob9",34))
df16 <- data.frame(cbind(x,y,plants))
colnames(df16) <- c("pos","depth","plants")

png("OR16_CDS_SNP.png",units="in", width=13, height=6,res=300)
ggplot(df16, aes(x=factor(x), y=y,group=x)) + geom_boxplot()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="SNPs in OR16 coding regions")
dev.off()

or13_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR13_CDS.depth",sep = "\t", header=FALSE)
or13 <- read.loci("~/Downloads/OR13.recode.vcf",header = FALSE)
colnames(or13_depth) <- names
nrow(or13_depth)
## 1752
or13_depth_snp <- or13_depth[or13_depth$pos %in% or13$V2,] 
## nrow(or13_depth_snp) = 83

x <- rep(or13_depth_snp$pos,20)
y <- c(or13_depth_snp$tom10,or13_depth_snp$tom1,or13_depth_snp$tom2,or13_depth_snp$tom3,or13_depth_snp$tom4,
       or13_depth_snp$tom5,or13_depth_snp$tom6,or13_depth_snp$tom7,or13_depth_snp$tom8,or13_depth_snp$tom9,
       or13_depth_snp$tob10,or13_depth_snp$tob1,or13_depth_snp$tob2,or13_depth_snp$tob3,or13_depth_snp$tob4,
       or13_depth_snp$tob5,or13_depth_snp$tob6,or13_depth_snp$tob7,or13_depth_snp$tob8,or13_depth_snp$tob9)
plants <- c(rep("tom10",83),rep("tom1",83),rep("tom2",83),rep("tom3",83),rep("tom4",83),rep("tom5",83),rep("tom6",83),rep("tom7",83),
            rep("tom8",83),rep("tom9",83),rep("tob10",83),rep("tob1",83),rep("tob2",83),rep("tob3",83),rep("tob4",83),
            rep("tob5",83),rep("tob6",83),rep("tob7",83),rep("tob8",83),rep("tob9",83))
df13 <- data.frame(cbind(x,y,plants))
colnames(df13) <- c("pos","depth","plants")

png("OR13_CDS_SNP.png",units="in", width=20, height=8,res=300)
ggplot(df13, aes(x=factor(x), y=y,group=x)) + geom_boxplot()+xlab("Loci: Scaffold NWSH01002174.1")+ylab("Depth")+
  labs(title="SNPs in OR13 coding regions") + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

orco_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/orco_CDS.depth",sep = "\t", header=FALSE)
orco <- read.loci("~/Downloads/orco.recode.vcf",header = FALSE)
colnames(orco_depth) <- names
nrow(orco_depth)
## 1419
orco_depth_snp <- orco_depth[orco_depth$pos %in% orco$V2,] 
## nrow(orco_depth_snp) = 53

x <- rep(orco_depth_snp$pos,20)
y <- c(orco_depth_snp$tom10,orco_depth_snp$tom1,orco_depth_snp$tom2,orco_depth_snp$tom3,orco_depth_snp$tom4,
       orco_depth_snp$tom5,orco_depth_snp$tom6,orco_depth_snp$tom7,orco_depth_snp$tom8,orco_depth_snp$tom9,
       orco_depth_snp$tob10,orco_depth_snp$tob1,orco_depth_snp$tob2,orco_depth_snp$tob3,orco_depth_snp$tob4,
       orco_depth_snp$tob5,orco_depth_snp$tob6,orco_depth_snp$tob7,orco_depth_snp$tob8,orco_depth_snp$tob9)
plants <- c(rep("tom10",53),rep("tom1",53),rep("tom2",53),rep("tom3",53),rep("tom4",53),rep("tom5",53),rep("tom6",53),rep("tom7",53),
            rep("tom8",53),rep("tom9",53),rep("tob10",53),rep("tob1",53),rep("tob2",53),rep("tob3",53),rep("tob4",53),
            rep("tob5",53),rep("tob6",53),rep("tob7",53),rep("tob8",53),rep("tob9",53))
dforco <- data.frame(cbind(x,y,plants))
colnames(dforco) <- c("pos","depth","plants")

png("orco_CDS_SNP.png",units="in", width=20, height=8,res=300)
ggplot(dforco, aes(x=factor(x), y=y,group=x)) + geom_boxplot()+xlab("Loci: Scaffold NWSH01000897.1")+ylab("Depth")+
  labs(title="SNPs in orco coding regions") + theme(axis.text.x = element_text(angle=45,hjust=1))
dev.off()

