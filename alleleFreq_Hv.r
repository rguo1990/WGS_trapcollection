## the script was used plot allele frequency of the nonsynonymous SNPs in the coding regions of ORs in Hv populations

library(foreign)
library("Biostrings")
library("pegas")
library("seqinr")
library("stringr")
library(ggplot2)

## OR6 Hv tobacco
f <- file.path("/Users/rongguo/Downloads/HvOR6_updated", c("Hvtob1.OR6_cDNA_1.fa","Hvtob1.OR6_cDNA_2.fa", 
                                                   "Hvtob3.OR6_cDNA_1.fa", "Hvtob3.OR6_cDNA_2.fa",
                                                   "Hvtob4.OR6_cDNA_1.fa", "Hvtob4.OR6_cDNA_2.fa",
                                                   "Hvtob5.OR6_cDNA_1.fa", "Hvtob5.OR6_cDNA_2.fa",
                                                   "Hvtob6.OR6_cDNA_1.fa", "Hvtob6.OR6_cDNA_2.fa",
                                                   "Hvtob7.OR6_cDNA_1.fa", "Hvtob7.OR6_cDNA_2.fa",
                                                   "Hvtob8.OR6_cDNA_1.fa", "Hvtob8.OR6_cDNA_2.fa",
                                                   "Hvtob9.OR6_cDNA_1.fa", "Hvtob9.OR6_cDNA_2.fa",
                                                   "Hvtob10.OR6_cDNA_1.fa", "Hvtob10.OR6_cDNA_2.fa"))
or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}
a1=a2=a3=NULL
for (i in 1:18){
  a1 = c(a1,substr(or[i],275*3+1,275*3+1))
  a2 = c(a2,substr(or[i],275*3+2,275*3+2))
  a3 = c(a3,substr(or[i],275*3+3,275*3+3))
}
## only 1 position: 267
positions <- c("267_1","267_1","267_1","267_1","267_2","267_2","267_2","267_2","267_3","267_3","267_3","267_3")
agct <- c(table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
          table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
          table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
d=cbind(rep(c("A","G","C","T"),3),agct,positions)
row.names(d) <- 1:12
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")

png("AF_OR6_Hvtob.png",units="in", width=7, height=8,res=300)
ggplot(df,aes(positions,agct,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR6 in Hv tobacco population")
dev.off()

## OR6 Hv tomatillo
f <- file.path("/Users/rongguo/Downloads/HvOR6_updated/", c("Hvtom1.OR6_cDNA_1.fa","Hvtom1.OR6_cDNA_2.fa", 
                                                   "Hvtom2.OR6_cDNA_1.fa","Hvtom2.OR6_cDNA_2.fa",
                                                   "Hvtom4.OR6_cDNA_1.fa","Hvtom4.OR6_cDNA_2.fa",
                                                   "Hvtom5.OR6_cDNA_1.fa","Hvtom5.OR6_cDNA_2.fa",
                                                   "Hvtom6.OR6_cDNA_1.fa","Hvtom6.OR6_cDNA_2.fa",
                                                   "Hvtom7.OR6_cDNA_1.fa","Hvtom7.OR6_cDNA_2.fa",
                                                   "Hvtom8.OR6_cDNA_1.fa","Hvtom8.OR6_cDNA_2.fa",
                                                   "Hvtom9.OR6_cDNA_1.fa","Hvtom9.OR6_cDNA_2.fa",
                                                   "Hvtom10.OR6_cDNA_1.fa","Hvtom10.OR6_cDNA_2.fa"))
or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}
a1=a2=a3=NULL
for (i in 1:18){
  a1 = c(a1,substr(or[i],275*3+1,275*3+1))
  a2 = c(a2,substr(or[i],275*3+2,275*3+2))
  a3 = c(a3,substr(or[i],275*3+3,275*3+3))
}
positions <- c("267_1","267_1","267_1","267_1","267_2","267_2","267_2","267_2","267_3","267_3","267_3","267_3")
agct <- c(table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
          table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
          table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
d=cbind(rep(c("A","G","C","T"),3),agct,positions)
row.names(d) <- 1:12
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")

png("AF_OR6_Hvtom.png",units="in", width=7, height=8,res=300)
ggplot(df,aes(positions,agct,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency in Hv tomatillo population")
dev.off()

## Fisher-exact test
Group1 <- c(11,7,3,15)
Group2 <- c(12,6,2,16)
f <- cbind(Group1,Group2)
fisher.test(f[1:2,])$p.value
## p-value=1
fisher.test(f[3:4,])$p.value
## p-value=1

## OR14 Hv tobacco
f <- file.path("/Users/rongguo/Downloads/HvOR14_updated", c("Hvtob1.OR14_cDNA_1.fa","Hvtob1.OR14_cDNA_2.fa", 
                                                           "Hvtob3.OR14_cDNA_1.fa", "Hvtob3.OR14_cDNA_2.fa",
                                                           "Hvtob4.OR14_cDNA_1.fa", "Hvtob4.OR14_cDNA_2.fa",
                                                           "Hvtob5.OR14_cDNA_1.fa", "Hvtob5.OR14_cDNA_2.fa",
                                                           "Hvtob6.OR14_cDNA_1.fa", "Hvtob6.OR14_cDNA_2.fa",
                                                           "Hvtob7.OR14_cDNA_1.fa", "Hvtob7.OR14_cDNA_2.fa",
                                                           "Hvtob8.OR14_cDNA_1.fa", "Hvtob8.OR14_cDNA_2.fa",
                                                           "Hvtob9.OR14_cDNA_1.fa", "Hvtob9.OR14_cDNA_2.fa",
                                                           "Hvtob10.OR14_cDNA_1.fa", "Hvtob10.OR14_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

## 20 diversity positions in total
var=c(49,56,57,74,82,123,135,167,180,187,198,207,212,229,264,267,311,314,348,355)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
    }
    agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
              table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
              table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])

}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*20),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(1,17,12,6,12,2,4,9,9,10,8,9,9,16,2,10,8,3,15,17,1,9,9,6,12,10,8,7,11,4,14,4,14,1,17,4,14,15,3,11,7,13,5,6,12)
Group1 <- Frequency

png("AF_OR14_Hvtob.png",units="in", width=16, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR14 in Hv tobacco population")
dev.off()

## OR14 Hv tomatillo
f <- file.path("/Users/rongguo/Downloads/HvOR14_updated", c("Hvtom1.OR14_cDNA_1.fa","Hvtom1.OR14_cDNA_2.fa",
                                                            "Hvtom2.OR14_cDNA_1.fa","Hvtom2.OR14_cDNA_2.fa",
                                                            "Hvtom4.OR14_cDNA_1.fa","Hvtom4.OR14_cDNA_2.fa",
                                                            "Hvtom5.OR14_cDNA_1.fa","Hvtom5.OR14_cDNA_2.fa",
                                                            "Hvtom6.OR14_cDNA_1.fa","Hvtom6.OR14_cDNA_2.fa",
                                                            "Hvtom7.OR14_cDNA_1.fa","Hvtom7.OR14_cDNA_2.fa",
                                                            "Hvtom8.OR14_cDNA_1.fa","Hvtom8.OR14_cDNA_2.fa",
                                                            "Hvtom9.OR14_cDNA_1.fa","Hvtom9.OR14_cDNA_2.fa",
                                                            "Hvtom10.OR14_cDNA_1.fa","Hvtom10.OR14_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

var=c(49,56,57,74,82,123,135,167,180,187,198,207,212,229,264,267,311,314,348,355)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
  }
  agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
            table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
            table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
  
}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*20),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(6,12,9,9,8,2,8,8,10,10,8,10,8,16,2,10,8,4,14,13,5,10,8,7,11,11,7,8,10,1,17,5,13,3,15,2,16,10,8,15,3,12,6,4,14)
Group2 <- Frequency

png("AF_OR14_Hvtom.png",units="in", width=16, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR14 in Hv tomatillo population")
dev.off()

## Fisher-exact test
f <- cbind(Group1,Group2)
fisher.test(f[1:2,])$p.value
## p-value=0.08768328
fisher.test(f[3:4,])$p.value
## p-value=0.499792
fisher.test(f[5:7,])$p.value
## p-value=0.3574862
fisher.test(f[8:9,])$p.value
## p-value=1
fisher.test(f[10:11,])$p.value
## p-value=1
fisher.test(f[12:13,])$p.value
## p-value=1
fisher.test(f[14:15,])$p.value
## p-value=1
fisher.test(f[16:17,])$p.value
## p-value=1
fisher.test(f[18:19,])$p.value
## p-value=1
fisher.test(f[20:21,])$p.value
## p-value=0.1774194
fisher.test(f[22:23,])$p.value
## p-value=1
fisher.test(f[24:25,])$p.value
## p-value=1
fisher.test(f[26:27,])$p.value
## p-value=1
fisher.test(f[28:29,])$p.value
## p-value=1
fisher.test(f[30:31,])$p.value
## p-value=0.3376623
fisher.test(f[32:33,])$p.value
## p-value=1
fisher.test(f[34:35,])$p.value
## p-value=0.6025974
fisher.test(f[36:37,])$p.value
## p-value=0.6581483
fisher.test(f[38:39,])$p.value
## p-value=0.1463531
fisher.test(f[40:41,])$p.value
## p-value=0.264233
fisher.test(f[42:43,])$p.value
## p-value=1
fisher.test(f[44:45,])$p.value
## p-value=0.7111943


## OR16 Hv tobacco
f <- file.path("/Users/rongguo/Downloads/HvOR16_updated", c("Hvtob1.OR16_cDNA_1.fa","Hvtob1.OR16_cDNA_2.fa", 
                                                            "Hvtob3.OR16_cDNA_1.fa", "Hvtob3.OR16_cDNA_2.fa",
                                                            "Hvtob4.OR16_cDNA_1.fa", "Hvtob4.OR16_cDNA_2.fa",
                                                            "Hvtob5.OR16_cDNA_1.fa", "Hvtob5.OR16_cDNA_2.fa",
                                                            "Hvtob6.OR16_cDNA_1.fa", "Hvtob6.OR16_cDNA_2.fa",
                                                            "Hvtob7.OR16_cDNA_1.fa", "Hvtob7.OR16_cDNA_2.fa",
                                                            "Hvtob8.OR16_cDNA_1.fa", "Hvtob8.OR16_cDNA_2.fa",
                                                            "Hvtob9.OR16_cDNA_1.fa", "Hvtob9.OR16_cDNA_2.fa",
                                                            "Hvtob10.OR16_cDNA_1.fa", "Hvtob10.OR16_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

## 9 diversity positions in total
var=c(53,259,260,264,282,290,292,300,307)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
  }
  agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
            table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
            table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
  
}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*l),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(13,5,16,2,5,13,5,13,12,6,1,17,16,2,10,8,10,8,16,2,3,15)
Group1 <- Frequency

png("AF_OR16_Hvtob.png",units="in", width=14, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR16 in Hv tobacco population")
dev.off()

## OR16 Hv tomatillo
f <- file.path("/Users/rongguo/Downloads/HvOR16_updated", c("Hvtom1.OR16_cDNA_1.fa","Hvtom1.OR16_cDNA_2.fa",
                                                            "Hvtom2.OR16_cDNA_1.fa","Hvtom2.OR16_cDNA_2.fa",
                                                            "Hvtom4.OR16_cDNA_1.fa","Hvtom4.OR16_cDNA_2.fa",
                                                            "Hvtom5.OR16_cDNA_1.fa","Hvtom5.OR16_cDNA_2.fa",
                                                            "Hvtom6.OR16_cDNA_1.fa","Hvtom6.OR16_cDNA_2.fa",
                                                            "Hvtom7.OR16_cDNA_1.fa","Hvtom7.OR16_cDNA_2.fa",
                                                            "Hvtom8.OR16_cDNA_1.fa","Hvtom8.OR16_cDNA_2.fa",
                                                            "Hvtom9.OR16_cDNA_1.fa","Hvtom9.OR16_cDNA_2.fa",
                                                            "Hvtom10.OR16_cDNA_1.fa","Hvtom10.OR16_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

var=c(53,259,260,264,282,290,292,300,307)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
  }
  agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
            table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
            table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
  
}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*l),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(10,8,15,3,3,15,3,15,10,8,4,14,13,5,12,6,12,6,13,5,1,17)
Group2 <- Frequency

png("AF_OR16_Hvtom.png",units="in", width=14, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR16 in Hv tomatillo population")
dev.off()

## Fisher-exact test
f <- cbind(Group1,Group2)
fisher.test(f[1:2,])$p.value
## p-value=0.4886763
fisher.test(f[3:4,])$p.value
## p-value=1
fisher.test(f[5:6,])$p.value
## p-value=0.6905653
fisher.test(f[7:8,])$p.value
## p-value=0.6905653
fisher.test(f[9:10,])$p.value
## p-value=0.7332224
fisher.test(f[11:12,])$p.value
## p-value=0.3376623
fisher.test(f[13:14,])$p.value
## p-value=0.4017595
fisher.test(f[15:16,])$p.value
## p-value=0.7332224
fisher.test(f[17:18,])$p.value
## p-value=0.7332224
fisher.test(f[19:20,])$p.value
## p-value=0.4017595
fisher.test(f[21:22,])$p.value
## p-value=0.6025974

## OR13 Hv tobacco
f <- file.path("/Users/rongguo/Downloads/HvOR13_updated", c("Hvtob1.OR13_cDNA_1.fa","Hvtob1.OR13_cDNA_2.fa", 
                                                            "Hvtob3.OR13_cDNA_1.fa", "Hvtob3.OR13_cDNA_2.fa",
                                                            "Hvtob4.OR13_cDNA_1.fa", "Hvtob4.OR13_cDNA_2.fa",
                                                            "Hvtob5.OR13_cDNA_1.fa", "Hvtob5.OR13_cDNA_2.fa",
                                                            "Hvtob6.OR13_cDNA_1.fa", "Hvtob6.OR13_cDNA_2.fa",
                                                            "Hvtob7.OR13_cDNA_1.fa", "Hvtob7.OR13_cDNA_2.fa",
                                                            "Hvtob8.OR13_cDNA_1.fa", "Hvtob8.OR13_cDNA_2.fa",
                                                            "Hvtob9.OR13_cDNA_1.fa", "Hvtob9.OR13_cDNA_2.fa",
                                                            "Hvtob10.OR13_cDNA_1.fa", "Hvtob10.OR13_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

## 13 diversity positions in total
var=c(60,64,74,81,334,343,411,416,423,429,432,463,546)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
  }
  agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
            table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
            table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
  
}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*l),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(9,9,6,12,6,12,11,7,10,8,11,7,1,17,16,2,15,3,7,11,11,7,16,2,17,1,14,4,3,15)
Group1 <- Frequency

png("AF_OR13_Hvtob.png",units="in", width=8, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR13 in Hv tobacco population")
dev.off()

## OR13 Hv tomatillo
f <- file.path("/Users/rongguo/Downloads/HvOR13_updated", c("Hvtom1.OR13_cDNA_1.fa","Hvtom1.OR13_cDNA_2.fa",
                                                            "Hvtom2.OR13_cDNA_1.fa","Hvtom2.OR13_cDNA_2.fa",
                                                            "Hvtom4.OR13_cDNA_1.fa","Hvtom4.OR13_cDNA_2.fa",
                                                            "Hvtom5.OR13_cDNA_1.fa","Hvtom5.OR13_cDNA_2.fa",
                                                            "Hvtom6.OR13_cDNA_1.fa","Hvtom6.OR13_cDNA_2.fa",
                                                            "Hvtom7.OR13_cDNA_1.fa","Hvtom7.OR13_cDNA_2.fa",
                                                            "Hvtom8.OR13_cDNA_1.fa","Hvtom8.OR13_cDNA_2.fa",
                                                            "Hvtom9.OR13_cDNA_1.fa","Hvtom9.OR13_cDNA_2.fa",
                                                            "Hvtom10.OR13_cDNA_1.fa","Hvtom10.OR13_cDNA_2.fa"))

or=NULL
for (i in 1:18){
  fa <- readDNAStringSet(f[i])
  or[i] = paste(fa)
}

var=c(60,64,74,81,334,343,411,416,423,429,432,463,546)
agct=NULL
l=length(var)
for (j in 1:l){
  a1=a2=a3=NULL
  for (i in 1:18){
    a1 = c(a1,substr(or[i],(var[j]-1)*3+1,(var[j]-1)*3+1))
    a2 = c(a2,substr(or[i],(var[j]-1)*3+2,(var[j]-1)*3+2))
    a3 = c(a3,substr(or[i],(var[j]-1)*3+3,(var[j]-1)*3+3))
  }
  agct <- c(agct,table(a1)["A"],table(a1)["G"],table(a1)["C"],table(a1)["T"],
            table(a2)["A"],table(a2)["G"],table(a2)["C"],table(a2)["T"],
            table(a3)["A"],table(a3)["G"],table(a3)["C"],table(a3)["T"])
  
}
positions=NULL
for (i in 1:l){
  positions <- c(positions,rep(paste(as.character(var[i]),"_1",sep=""),4),rep(paste(as.character(var[i]),"_2",sep=""),4),
                 rep(paste(as.character(var[i]),"_3",sep=""),4))
}
d=cbind(rep(c("A","G","C","T"),3*l),agct,positions)
row.names(d) <- 1:nrow(d)
df<-data.frame(d)
colnames(df) <- c("AGCT","Freq","Pos")
df_new <- subset(df,Freq!=18)
## df_new$Freq
Frequency <- c(5,13,4,14,4,14,15,3,4,14,10,8,4,14,15,3,16,2,14,4,4,14,13,5,12,6,12,6,5,13)
Group2 <- Frequency

png("AF_OR13_Hvtom.png",units="in", width=8, height=8,res=300)
ggplot(df_new,aes(Pos,Frequency,fill=AGCT))+geom_bar(stat="identity",width = 0.5)+xlab("positions")+
  scale_y_continuous(name="Frequency",breaks = 1:20)+labs(title = "Nonsynonymous SNP allele frequency of OR13 in Hv tomatillo population")
dev.off()

## Fisher-exact test
f <- cbind(Group1,Group2)
fisher.test(f[1:2,])$p.value
## p-value=0.3052667
fisher.test(f[3:4,])$p.value
## p-value=0.7111943
fisher.test(f[5:6,])$p.value
## p-value=0.7111943
fisher.test(f[7:8,])$p.value
## p-value=0.264233
fisher.test(f[9:10,])$p.value
## p-value=0.08580226
fisher.test(f[11:12,])$p.value
## p-value=1
fisher.test(f[13:14,])$p.value
## p-value=0.3376623
fisher.test(f[15:16,])$p.value
## p-value=1
fisher.test(f[17:18,])$p.value
## p-value=1
fisher.test(f[19:20,])$p.value
## p-value=0.0409118
fisher.test(f[21:22,])$p.value
## p-value=0.0409118
fisher.test(f[23:24,])$p.value
## p-value=0.4017595
fisher.test(f[25:26,])$p.value
## p-value=0.08768328
fisher.test(f[27:28,])$p.value
## p-value=0.7111943
fisher.test(f[29:30,])$p.value
## p-value=0.6905653

