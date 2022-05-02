## the script was used to barplot the average read alignment depth of OR genes (Hs population)
## x axis: gene loci
## y axis: average read alignment depth among 20 Hs individuals

library(ggplot2)
## OR6_gene.depth was made by samtools depth
or6_depth <- read.csv("/Users/rongguo/Downloads/Hs/OR6_gene.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or6_depth) <- names
l=nrow(or6_depth)
## 7405
a=NULL
for (i in 1:l){
  y <- c(or6_depth[i,]$tom10,or6_depth[i,]$tom1,or6_depth[i,]$tom2,or6_depth[i,]$tom3,or6_depth[i,]$tom4,
         or6_depth[i,]$tom5,or6_depth[i,]$tom6,or6_depth[i,]$tom7,or6_depth[i,]$tom8,or6_depth[i,]$tom9,
         or6_depth[i,]$tob10,or6_depth[i,]$tob1,or6_depth[i,]$tob2,or6_depth[i,]$tob3,or6_depth[i,]$tob4,
         or6_depth[i,]$tob5,or6_depth[i,]$tob6,or6_depth[i,]$tob7,or6_depth[i,]$tob8,or6_depth[i,]$tob9)
  a <- c(a, mean(y))
}
loci <- or6_depth$pos
cds <- c(21136:21176, 19954:20042, 19633:19802, 19056:19166, 18496:18728, 16941:17156, 16455:16557, 16313:16375,
         16145:16177, 15655:15728, 15449:15530, 10583:10642) 
IfCDS <- or6_depth$pos %in% cds
df6 <- data.frame(cbind(loci,a,IfCDS))
colnames(df6) <- c("position","mean_depth","CDS")

png("OR6_Hs_Depth.png",units="in", width=15, height=8,res=300)
ggplot(df6,aes(loci,a,fill=IfCDS))+geom_bar(stat="identity")+xlab("positions in Scaffold_4031_rc.Scaffold_489")+
  ylab("Mean Depth") + labs(title="Average read alignment depth of gene OR6")
dev.off()

or13_depth <- read.csv("/Users/rongguo/Downloads/Hs/OR13_gene.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or13_depth) <- names
l=nrow(or13_depth)
## 3742
a=NULL
for (i in 1:l){
  y <- c(or13_depth[i,]$tom10,or13_depth[i,]$tom1,or13_depth[i,]$tom2,or13_depth[i,]$tom3,or13_depth[i,]$tom4,
         or13_depth[i,]$tom5,or13_depth[i,]$tom6,or13_depth[i,]$tom7,or13_depth[i,]$tom8,or13_depth[i,]$tom9,
         or13_depth[i,]$tob10,or13_depth[i,]$tob1,or13_depth[i,]$tob2,or13_depth[i,]$tob3,or13_depth[i,]$tob4,
         or13_depth[i,]$tob5,or13_depth[i,]$tob6,or13_depth[i,]$tob7,or13_depth[i,]$tob8,or13_depth[i,]$tob9)
  a <- c(a, mean(y))
}
loci <- or13_depth$pos
cds <- c(135623:135654, 135862:136099, 136920:137027, 137199:137680, 137904:138006, 138158:138253, 
         138520:138675,140086:140148) 
IfCDS <- or13_depth$pos %in% cds
df13 <- data.frame(cbind(loci,a,IfCDS))
colnames(df13) <- c("position","mean_depth","CDS")

png("OR13_Hs_Depth.png",units="in", width=15, height=8,res=300)
ggplot(df13,aes(loci,a,fill=IfCDS))+geom_bar(stat="identity")+xlab("positions in Scaffold_2038")+
  ylab("Mean Depth") + labs(title="Average read alignment depth of gene OR13")
dev.off()

or14_depth <- read.csv("/Users/rongguo/Downloads/Hs/OR14_gene.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or14_depth) <- names
l=nrow(or14_depth)
## 4219
a=NULL
for (i in 1:l){
  y <- c(or14_depth[i,]$tom10,or14_depth[i,]$tom1,or14_depth[i,]$tom2,or14_depth[i,]$tom3,or14_depth[i,]$tom4,
         or14_depth[i,]$tom5,or14_depth[i,]$tom6,or14_depth[i,]$tom7,or14_depth[i,]$tom8,or14_depth[i,]$tom9,
         or14_depth[i,]$tob10,or14_depth[i,]$tob1,or14_depth[i,]$tob2,or14_depth[i,]$tob3,or14_depth[i,]$tob4,
         or14_depth[i,]$tob5,or14_depth[i,]$tob6,or14_depth[i,]$tob7,or14_depth[i,]$tob8,or14_depth[i,]$tob9)
  a <- c(a, mean(y))
}
loci <- or14_depth$pos
cds <- c(173309:173352, 173687:173963, 174468:174578, 174735:174967, 175197:175451, 175758:175860, 
         179832:179853, 179991:180064, 180218:180373, 182072:182119) 
IfCDS <- or14_depth$pos %in% cds
df14 <- data.frame(cbind(loci,a,IfCDS))
colnames(df14) <- c("position","mean_depth","CDS")

png("OR14_Hs_Depth.png",units="in", width=15, height=8,res=300)
ggplot(df14,aes(loci,a,fill=IfCDS))+geom_bar(stat="identity")+xlab("positions in Scaffold_656")+
  ylab("Mean Depth") + labs(title="Average read alignment depth of gene OR14")
dev.off()

or16_depth <- read.csv("/Users/rongguo/Downloads/Hs/OR16_gene.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(or16_depth) <- names
l=nrow(or16_depth)
## 5253
a=NULL
for (i in 1:l){
  y <- c(or16_depth[i,]$tom10,or16_depth[i,]$tom1,or16_depth[i,]$tom2,or16_depth[i,]$tom3,or16_depth[i,]$tom4,
         or16_depth[i,]$tom5,or16_depth[i,]$tom6,or16_depth[i,]$tom7,or16_depth[i,]$tom8,or16_depth[i,]$tom9,
         or16_depth[i,]$tob10,or16_depth[i,]$tob1,or16_depth[i,]$tob2,or16_depth[i,]$tob3,or16_depth[i,]$tob4,
         or16_depth[i,]$tob5,or16_depth[i,]$tob6,or16_depth[i,]$tob7,or16_depth[i,]$tob8,or16_depth[i,]$tob9)
  a <- c(a, mean(y))
}
loci <- or16_depth$pos
cds <- c(206203:206243, 205652:205910, 205127:205237, 204102:204334, 203280:203501, 202381:202483,
         202210:202305, 201615:201770, 199926:199973) 
IfCDS <- or16_depth$pos %in% cds
df16 <- data.frame(cbind(loci,a,IfCDS))
colnames(df16) <- c("position","mean_depth","CDS")

png("OR16_Hs_Depth.png",units="in", width=15, height=8,res=300)
ggplot(df16,aes(loci,a,fill=IfCDS))+geom_bar(stat="identity")+xlab("positions in Scaffold_656")+
  ylab("Mean Depth") + labs(title="Average read alignment depth of gene OR16")
dev.off()

orco_depth <- read.csv("/Users/rongguo/Downloads/Hs/orco_gene.depth",sep = "\t",header = FALSE)
names <- c("Scaffold","pos","tom10","tom1","tom2","tom3","tom4","tom5","tom6","tom7","tom8","tom9",
           "tob10","tob1","tob2","tob3","tob4","tob5","tob6","tob7","tob8","tob9")
colnames(orco_depth) <- names
l=nrow(orco_depth)
## 6950
a=NULL
for (i in 1:l){
  y <- c(orco_depth[i,]$tom10,orco_depth[i,]$tom1,orco_depth[i,]$tom2,orco_depth[i,]$tom3,orco_depth[i,]$tom4,
         orco_depth[i,]$tom5,orco_depth[i,]$tom6,orco_depth[i,]$tom7,orco_depth[i,]$tom8,orco_depth[i,]$tom9,
         orco_depth[i,]$tob10,orco_depth[i,]$tob1,orco_depth[i,]$tob2,orco_depth[i,]$tob3,orco_depth[i,]$tob4,
         orco_depth[i,]$tob5,orco_depth[i,]$tob6,orco_depth[i,]$tob7,orco_depth[i,]$tob8,orco_depth[i,]$tob9)
  a <- c(a, mean(y))
}
loci <- orco_depth$pos
cds <- c(54602:54701, 53477:53687, 49768:49904, 46998:47158, 46320:46506, 45656:45866, 41461:41560, 
         40055:40159, 39205:39360, 38994:39044) 
IfCDS <- orco_depth$pos %in% cds
dforco <- data.frame(cbind(loci,a,IfCDS))
colnames(dforco) <- c("position","mean_depth","CDS")

png("orco_Hs_Depth.png",units="in", width=15, height=8,res=300)
ggplot(dforco,aes(loci,a,fill=IfCDS))+geom_bar(stat="identity")+xlab("positions in Scaffold_8")+
  ylab("Mean Depth") + labs(title="Average read alignment depth of gene orco")
dev.off()
