## the script was used to plot read alignment depth of the coding regions of OR genes
## where x axis is positions, y axis is the read alignment depth
## .depth files were generated by samtools.

library(ggplot2)
library(reshape2)
or6_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR6_CDS.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or6_depth) <- names
nrow(or6_depth)
## 1236 
ggplot(or6_depth,aes(pos,tob4))+geom_line(color="red")

loci <- rep(or6_depth$pos,20)
y <- c(or6_depth$tom10,or6_depth$tom1,or6_depth$tom2,or6_depth$tom3,or6_depth$tom4,or6_depth$tom5,or6_depth$tom6,or6_depth$tom7,
       or6_depth$tom8,or6_depth$tom9,or6_depth$tob10,or6_depth$tob1,or6_depth$tob2,or6_depth$tob3,or6_depth$tob4,or6_depth$tob5,
       or6_depth$tob6,or6_depth$tob7,or6_depth$tob8,or6_depth$tob9)
plants <- c(rep("tom10",1236),rep("tom1",1236),rep("tom2",1236),rep("tom3",1236),rep("tom4",1236),rep("tom5",1236),rep("tom6",1236),rep("tom7",1236),
            rep("tom8",1236),rep("tom9",1236),rep("tob10",1236),rep("tob1",1236),rep("tob2",1236),rep("tob3",1236),rep("tob4",1236),
            rep("tob5",1236),rep("tob6",1236),rep("tob7",1236),rep("tob8",1236),rep("tob9",1236))
df6 <- data.frame(cbind(loci,y,plants))
colnames(df6) = c("Loci","Depth","Group")
          
png("OR6_CDS.png",units="in", width=8, height=6,res=300)
ggplot(df6,aes(x=loci,y=y,colour = plants))+geom_point()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="Read alignment depth of OR6 coding regions")
dev.off()

or14_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR14_CDS.depth",sep = "\t", header=FALSE)
colnames(or14_depth) <- names
nrow(or14_depth)
## 1125 

loci <- rep(or14_depth$pos,20)
y <- c(or14_depth$tom10,or14_depth$tom1,or14_depth$tom2,or14_depth$tom3,or14_depth$tom4,or14_depth$tom5,or14_depth$tom6,or14_depth$tom7,
       or14_depth$tom8,or14_depth$tom9,or14_depth$tob10,or14_depth$tob1,or14_depth$tob2,or14_depth$tob3,or14_depth$tob4,or14_depth$tob5,
       or14_depth$tob6,or14_depth$tob7,or14_depth$tob8,or14_depth$tob9)
plants <- c(rep("tom10",1125),rep("tom1",1125),rep("tom2",1125),rep("tom3",1125),rep("tom4",1125),rep("tom5",1125),rep("tom6",1125),rep("tom7",1125),
            rep("tom8",1125),rep("tom9",1125),rep("tob10",1125),rep("tob1",1125),rep("tob2",1125),rep("tob3",1125),rep("tob4",1125),
            rep("tob5",1125),rep("tob6",1125),rep("tob7",1125),rep("tob8",1125),rep("tob9",1125))
df14 <- data.frame(cbind(loci,y,plants))
colnames(df14) = c("Loci","Depth","Group")

png("OR14_CDS.png",units="in", width=8, height=6,res=300)
ggplot(df14,aes(x=loci,y=y,colour = plants))+geom_point()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="Read alignment depth of OR14 coding regions")
dev.off()

or16_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR16_CDS.depth",sep = "\t",header = FALSE)
colnames(or16_depth) <- names
nrow(or16_depth)
## 969
loci <- rep(or16_depth$pos,20)
y <- c(or16_depth$tom10,or16_depth$tom1,or16_depth$tom2,or16_depth$tom3,or16_depth$tom4,or16_depth$tom5,or16_depth$tom6,or16_depth$tom7,
       or16_depth$tom8,or16_depth$tom9,or16_depth$tob10,or16_depth$tob1,or16_depth$tob2,or16_depth$tob3,or16_depth$tob4,or16_depth$tob5,
       or16_depth$tob6,or16_depth$tob7,or16_depth$tob8,or16_depth$tob9)
plants <- c(rep("tom10",969),rep("tom1",969),rep("tom2",969),rep("tom3",969),rep("tom4",969),rep("tom5",969),rep("tom6",969),rep("tom7",969),
            rep("tom8",969),rep("tom9",969),rep("tob10",969),rep("tob1",969),rep("tob2",969),rep("tob3",969),rep("tob4",969),
            rep("tob5",969),rep("tob6",969),rep("tob7",969),rep("tob8",969),rep("tob9",969))
df16 <- data.frame(cbind(loci,y,plants))
colnames(df16) = c("Loci","Depth","Group")

png("OR16_CDS.png",units="in", width=8, height=6,res=300)
ggplot(df16,aes(x=loci,y=y,colour = plants))+geom_point()+xlab("Loci: Scaffold NWSH01000007.1")+ylab("Depth")+
  labs(title="Read alignment depth of OR16 coding regions")
dev.off()

or13_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/OR13_CDS.depth",sep = "\t",header = FALSE)
colnames(or13_depth) <- names
nrow(or13_depth)
## 1752
loci <- rep(or13_depth$pos,20)
y <- c(or13_depth$tom10,or13_depth$tom1,or13_depth$tom2,or13_depth$tom3,or13_depth$tom4,or13_depth$tom5,or13_depth$tom6,or13_depth$tom7,
       or13_depth$tom8,or13_depth$tom9,or13_depth$tob10,or13_depth$tob1,or13_depth$tob2,or13_depth$tob3,or13_depth$tob4,or13_depth$tob5,
       or13_depth$tob6,or13_depth$tob7,or13_depth$tob8,or13_depth$tob9)
plants <- c(rep("tom10",1752),rep("tom1",1752),rep("tom2",1752),rep("tom3",1752),rep("tom4",1752),rep("tom5",1752),rep("tom6",1752),rep("tom7",1752),
            rep("tom8",1752),rep("tom9",1752),rep("tob10",1752),rep("tob1",1752),rep("tob2",1752),rep("tob3",1752),rep("tob4",1752),
            rep("tob5",1752),rep("tob6",1752),rep("tob7",1752),rep("tob8",1752),rep("tob9",1752))
df13 <- data.frame(cbind(loci,y,plants))
colnames(df13) = c("Loci","Depth","Group")

png("OR13_CDS.png",units="in", width=8, height=6,res=300)
ggplot(df13,aes(x=loci,y=y,colour = plants))+geom_point()+xlab("Loci: Scaffold NWSH01002174.1")+ylab("Depth")+
  labs(title="Read alignment depth of OR13 coding regions")
dev.off()

orco_depth <- read.csv("/Users/rongguo/UMD/WGS/Depth/orco_CDS.depth",sep = "\t",header = FALSE)
colnames(orco_depth) <- names
nrow(orco_depth)
## 1419
loci <- rep(orco_depth$pos,20)
y <- c(orco_depth$tom10,orco_depth$tom1,orco_depth$tom2,orco_depth$tom3,orco_depth$tom4,orco_depth$tom5,orco_depth$tom6,orco_depth$tom7,
       orco_depth$tom8,orco_depth$tom9,orco_depth$tob10,orco_depth$tob1,orco_depth$tob2,orco_depth$tob3,orco_depth$tob4,orco_depth$tob5,
       orco_depth$tob6,orco_depth$tob7,orco_depth$tob8,orco_depth$tob9)
plants <- c(rep("tom10",1419),rep("tom1",1419),rep("tom2",1419),rep("tom3",1419),rep("tom4",1419),rep("tom5",1419),rep("tom6",1419),rep("tom7",1419),
            rep("tom8",1419),rep("tom9",1419),rep("tob10",1419),rep("tob1",1419),rep("tob2",1419),rep("tob3",1419),rep("tob4",1419),
            rep("tob5",1419),rep("tob6",1419),rep("tob7",1419),rep("tob8",1419),rep("tob9",1419))
dforco <- data.frame(cbind(loci,y,plants))
colnames(dforco) = c("Loci","Depth","Group")

png("orco_CDS.png",units="in", width=8, height=6,res=300)
ggplot(dforco,aes(x=loci,y=y,colour = plants))+geom_point()+xlab("Loci: Scaffold NWSH01000897.1")+ylab("Depth")+
  labs(title="Read alignment depth of orco coding regions")
dev.off()