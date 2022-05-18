## the script was used to extract the coding sequences of OR6, OR13, OR14, OR16 in the corresponding scaffolds from their vcf files, 
## and translate them into protein sequences

library("Biostrings")
library("pegas")
library("seqinr")
library("stringr")
fa <- readDNAStringSet("~/Downloads/Scaffold_4031_rc.Scaffold_489.fa")
sequence = paste(fa)

or6_vcf <- read.vcf("~/Downloads/Hs/OR6-Hs.recode.vcf")
or13_vcf <- read.vcf("~/Downloads/Hs/OR13-Hs.recode.vcf")
or14_vcf <- read.vcf("~/Downloads/Hs/OR14-Hs.recode.vcf")
or16_vcf <- read.vcf("~/Downloads/Hs/OR16-Hs.recode.vcf")

or6 <- read.loci("~/Downloads/Hs/OR6-Hs.recode.vcf",header = FALSE)
or13 <- read.loci("~/Downloads/Hs/OR13-Hs.recode.vcf",header = FALSE)
or14 <- read.loci("~/Downloads/Hs/OR14-Hs.recode.vcf",header = FALSE)
or16 <- read.loci("~/Downloads/Hs/OR16-Hs.recode.vcf",header = FALSE)

## OR6
or6_pos <- or6$V2
colnames(or6_vcf) <- or6_pos
or6_cds <- or6_vcf[,or6_pos %in% c(21136:21176, 19954:20042, 19633:19802, 19056:19166, 
                                   18496:18728, 16941:17156, 16455:16557, 16313:16375, 
                                   16145:16177, 15655:15728, 15449:15530, 10583:10642)]
positions <- colnames(or6_cds)
df6_cds <- data.frame(or6_cds)
colnames(df6_cds)=positions

individuals <- row.names(df6_cds)

## cDNA of OR6
Hs00=str_sub(sequence,1,25000)
a=as.character(reverseComplement(DNAString(subseq(Hs00,21136,21176))))
b=as.character(reverseComplement(DNAString(subseq(Hs00,19954,20042))))
c=as.character(reverseComplement(DNAString(subseq(Hs00,19633,19802))))
d=as.character(reverseComplement(DNAString(subseq(Hs00,19056,19166))))
e=as.character(reverseComplement(DNAString(subseq(Hs00,18496,18728))))
f=as.character(reverseComplement(DNAString(subseq(Hs00,16941,17156))))
g=as.character(reverseComplement(DNAString(subseq(Hs00,16455,16557))))
h=as.character(reverseComplement(DNAString(subseq(Hs00,16313,16375))))
k=as.character(reverseComplement(DNAString(subseq(Hs00,16145,16177))))
l=as.character(reverseComplement(DNAString(subseq(Hs00,15655,15728))))
m=as.character(reverseComplement(DNAString(subseq(Hs00,15449,15530))))
n=as.character(reverseComplement(DNAString(subseq(Hs00,10583,10642))))

ori_or6_cdna <- paste(a,b,c,d,e,f,g,h,k,l,m,n,sep="")
write.fasta(ori_or6_cdna, "OR6_cDNA", "~/Downloads/Hs/OR6_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR6
Hs_1=Hs00
Hs_2=Hs00
for (j in 1:20){
  for (i in positions){
    l=as.integer(i)
    Hs_1 <- paste(subseq(Hs_1,1,l-1),str_split(df6_cds[,i][j],"/")[[1]][1],subseq(Hs_1,l+1,25000),sep = "")
    Hs_2 <- paste(subseq(Hs_2,1,l-1),str_split(df6_cds[,i][j],"/")[[1]][2],subseq(Hs_2,l+1,25000),sep = "")
    }
    
    a=as.character(reverseComplement(DNAString(subseq(Hs_1,21136,21176))))
    b=as.character(reverseComplement(DNAString(subseq(Hs_1,19954,20042))))
    c=as.character(reverseComplement(DNAString(subseq(Hs_1,19633,19802))))
    d=as.character(reverseComplement(DNAString(subseq(Hs_1,19056,19166))))
    e=as.character(reverseComplement(DNAString(subseq(Hs_1,18496,18728))))
    f=as.character(reverseComplement(DNAString(subseq(Hs_1,16941,17156))))
    g=as.character(reverseComplement(DNAString(subseq(Hs_1,16455,16557))))
    h=as.character(reverseComplement(DNAString(subseq(Hs_1,16313,16375))))
    k=as.character(reverseComplement(DNAString(subseq(Hs_1,16145,16177))))
    o=as.character(reverseComplement(DNAString(subseq(Hs_1,15655,15728))))
    m=as.character(reverseComplement(DNAString(subseq(Hs_1,15449,15530))))
    n=as.character(reverseComplement(DNAString(subseq(Hs_1,10583,10642))))
    
    Hs_or6_cdna_1 <- paste(a,b,c,d,e,f,g,h,k,o,m,n,sep="")
    
    a=str_split(individuals[j],".sorted.")[[1]][1]
    name = paste(a,"_OR6_cDNA_1",sep="")
    path = paste("~/Downloads/Hs/",a,"_OR6_cDNA_1.fa",sep="")
  
    write.fasta(Hs_or6_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
    detach("package:seqinr", unload=TRUE)
    or6_protein_1 <- as.character(translate(DNAString(Hs_or6_cdna_1),if.fuzzy.codon="X"))
    name = paste(a,"_OR6_pep_1",sep="")
    path = paste("~/Downloads/Hs/",a,"_OR6_pep_1.fa",sep="")
    library("seqinr")
    write.fasta(or6_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    a=as.character(reverseComplement(DNAString(subseq(Hs_2,21136,21176))))
    b=as.character(reverseComplement(DNAString(subseq(Hs_2,19954,20042))))
    c=as.character(reverseComplement(DNAString(subseq(Hs_2,19633,19802))))
    d=as.character(reverseComplement(DNAString(subseq(Hs_2,19056,19166))))
    e=as.character(reverseComplement(DNAString(subseq(Hs_2,18496,18728))))
    f=as.character(reverseComplement(DNAString(subseq(Hs_2,16941,17156))))
    g=as.character(reverseComplement(DNAString(subseq(Hs_2,16455,16557))))
    h=as.character(reverseComplement(DNAString(subseq(Hs_2,16313,16375))))
    k=as.character(reverseComplement(DNAString(subseq(Hs_2,16145,16177))))
    o=as.character(reverseComplement(DNAString(subseq(Hs_2,15655,15728))))
    m=as.character(reverseComplement(DNAString(subseq(Hs_2,15449,15530))))
    n=as.character(reverseComplement(DNAString(subseq(Hs_2,10583,10642))))
  
    Hs_or6_cdna_2 <- paste(a,b,c,d,e,f,g,h,k,o,m,n,sep="")
    
    a=str_split(individuals[j],".sorted.")[[1]][1]
    name = paste(a,"_OR6_cDNA_2",sep="")
    path = paste("~/Downloads/Hs/",a,"_OR6_cDNA_2.fa",sep="")
    
    write.fasta(Hs_or6_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or6_protein_2 <- as.character(translate(DNAString(Hs_or6_cdna_2),if.fuzzy.codon="X"))
    name = paste(a,"_OR6_pep_2",sep="")
    path = paste("~/Downloads/Hs/",a,"_OR6_pep_2.fa",sep="")
    library("seqinr")
    write.fasta(or6_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}

## OR14
or14_pos <- or14$V2
colnames(or14_vcf) <- or14_pos
or14_cds <- or14_vcf[,or14_pos %in% c(173309:173352, 173687:173963, 174468:174578, 
                                      174735:174967, 175197:175451, 175758:175860, 
                                      179832:179853, 179991:180064, 180218:180373, 182072:182119)]
positions_2 <- colnames(or14_cds)
df14_cds <- data.frame(or14_cds)
colnames(df14_cds)=positions_2

fa <- readDNAStringSet("~/Downloads/Scaffold_656_new.fa")
sequence = paste(fa)

## cDNA of OR14
Hs00=str_sub(sequence,1,185000)
ori_or14_cdna <- paste(subseq(Hs00,173309,173352), subseq(Hs00,173687,173963), subseq(Hs00,174468,174578), 
                       subseq(Hs00,174735,174967), subseq(Hs00,175197,175451), subseq(Hs00,175758,175860), 
                       subseq(Hs00,179832,179853), subseq(Hs00,179991,180064), subseq(Hs00,180218,180373), 
                       subseq(Hs00,182072,182119),sep="")
write.fasta(ori_or14_cdna, "OR14_cDNA", "~/Downloads/Hs/OR14_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR14
Hs_1=Hs00
Hs_2=Hs00
positions_2=positions_2[! positions_2 %in% c("173899")]
for (j in 1:20){
  for (i in positions_2){
    l=as.integer(i)
    Hs_1 <- paste(subseq(Hs_1,1,l-1),str_split(df14_cds[,i][j],"/")[[1]][1],subseq(Hs_1,l+1,185000),sep = "")
    Hs_2 <- paste(subseq(Hs_2,1,l-1),str_split(df14_cds[,i][j],"/")[[1]][2],subseq(Hs_2,l+1,185000),sep = "")
  }
  
  Hs_or14_cdna_1 <- paste(subseq(Hs_1,173309,173352), subseq(Hs_1,173687,173963), subseq(Hs_1,174468,174578), 
                         subseq(Hs_1,174735,174967), subseq(Hs_1,175197,175451), subseq(Hs_1,175758,175860), 
                         subseq(Hs_1,179832,179853), subseq(Hs_1,179991,180064), subseq(Hs_1,180218,180373), 
                         subseq(Hs_1,182072,182119),sep="")

  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR14_cDNA_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR14_cDNA_1.fa",sep="")
  
  write.fasta(Hs_or14_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or14_protein_1 <- as.character(translate(DNAString(Hs_or14_cdna_1),if.fuzzy.codon="X"))
  name = paste(a,"_OR14_pep_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR14_pep_1.fa",sep="")
  library("seqinr")
  write.fasta(or14_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  Hs_or14_cdna_2 <- paste(subseq(Hs_2,173309,173352), subseq(Hs_2,173687,173963), subseq(Hs_2,174468,174578), 
                           subseq(Hs_2,174735,174967), subseq(Hs_2,175197,175451), subseq(Hs_2,175758,175860), 
                           subseq(Hs_2,179832,179853), subseq(Hs_2,179991,180064), subseq(Hs_2,180218,180373), 
                           subseq(Hs_2,182072,182119),sep="")
  
  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR14_cDNA_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR14_cDNA_2.fa",sep="")
  
  write.fasta(Hs_or14_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or14_protein_2 <- as.character(translate(DNAString(Hs_or14_cdna_2),if.fuzzy.codon="X"))
  name = paste(a,"_OR14_pep_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR14_pep_2.fa",sep="")
  library("seqinr")
  write.fasta(or14_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}


## OR16
or16_pos <- or16$V2
colnames(or16_vcf) <- or16_pos
or16_cds <- or16_vcf[,or16_pos %in% c(206203:206243, 205652:205910, 205127:205237, 
                                      204102:204334, 203280:203501, 202381:202483,
                                      202210:202305, 201615:201770, 199926:199973)]
positions_3 <- colnames(or16_cds)
df16_cds <- data.frame(or16_cds)
colnames(df16_cds)=positions_3

## cDNA of OR16
Hs00=str_sub(sequence,1,210000)
a=as.character(reverseComplement(DNAString(subseq(Hs00,206203,206243))))
b=as.character(reverseComplement(DNAString(subseq(Hs00,205652,205910))))
c=as.character(reverseComplement(DNAString(subseq(Hs00,205127,205237))))
d=as.character(reverseComplement(DNAString(subseq(Hs00,204102,204334))))
e=as.character(reverseComplement(DNAString(subseq(Hs00,203280,203501))))
f=as.character(reverseComplement(DNAString(subseq(Hs00,202381,202483))))
g=as.character(reverseComplement(DNAString(subseq(Hs00,202210,202305))))
h=as.character(reverseComplement(DNAString(subseq(Hs00,201615,201770))))
k=as.character(reverseComplement(DNAString(subseq(Hs00,199926,199973))))

ori_or16_cdna <- paste(a,b,c,d,e,f,g,h,k,sep="")
write.fasta(ori_or6_cdna, "OR16_cDNA", "~/Downloads/Hs/OR16_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR16
Hs_1=Hs00
Hs_2=Hs00
for (j in 1:20){
  for (i in positions_3){
    l=as.integer(i)
    if (str_split(df16_cds[,i][j],"/")[[1]][1]!="." && str_split(df16_cds[,i][j],"/")[[1]][2] !="."){
      Hs_1 <- paste(subseq(Hs_1,1,l-1),str_split(df16_cds[,i][j],"/")[[1]][1],subseq(Hs_1,l+1,210000),sep = "")
      Hs_2 <- paste(subseq(Hs_2,1,l-1),str_split(df16_cds[,i][j],"/")[[1]][2],subseq(Hs_2,l+1,210000),sep = "")
    }}
  
  a=as.character(reverseComplement(DNAString(subseq(Hs_1,206203,206243))))
  b=as.character(reverseComplement(DNAString(subseq(Hs_1,205652,205910))))
  c=as.character(reverseComplement(DNAString(subseq(Hs_1,205127,205237))))
  d=as.character(reverseComplement(DNAString(subseq(Hs_1,204102,204334))))
  e=as.character(reverseComplement(DNAString(subseq(Hs_1,203280,203501))))
  f=as.character(reverseComplement(DNAString(subseq(Hs_1,202381,202483))))
  g=as.character(reverseComplement(DNAString(subseq(Hs_1,202210,202305))))
  h=as.character(reverseComplement(DNAString(subseq(Hs_1,201615,201770))))
  k=as.character(reverseComplement(DNAString(subseq(Hs_1,199926,199973))))
  
  Hs_or16_cdna_1 <- paste(a,b,c,d,e,f,g,h,k,sep="")
  
  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR16_cDNA_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR16_cDNA_1.fa",sep="")
  
  write.fasta(Hs_or16_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or16_protein_1 <- as.character(translate(DNAString(Hs_or16_cdna_1),if.fuzzy.codon="X"))
  name = paste(a,"_OR16_pep_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR16_pep_1.fa",sep="")
  library("seqinr")
  write.fasta(or16_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  a=as.character(reverseComplement(DNAString(subseq(Hs_2,206203,206243))))
  b=as.character(reverseComplement(DNAString(subseq(Hs_2,205652,205910))))
  c=as.character(reverseComplement(DNAString(subseq(Hs_2,205127,205237))))
  d=as.character(reverseComplement(DNAString(subseq(Hs_2,204102,204334))))
  e=as.character(reverseComplement(DNAString(subseq(Hs_2,203280,203501))))
  f=as.character(reverseComplement(DNAString(subseq(Hs_2,202381,202483))))
  g=as.character(reverseComplement(DNAString(subseq(Hs_2,202210,202305))))
  h=as.character(reverseComplement(DNAString(subseq(Hs_2,201615,201770))))
  k=as.character(reverseComplement(DNAString(subseq(Hs_2,199926,199973))))
  
  Hs_or16_cdna_2 <- paste(a,b,c,d,e,f,g,h,k,sep="")
  
  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR16_cDNA_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR16_cDNA_2.fa",sep="")
  
  write.fasta(Hs_or16_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or16_protein_2 <- as.character(translate(DNAString(Hs_or16_cdna_2),if.fuzzy.codon="X"))
  name = paste(a,"_OR16_pep_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR16_pep_2.fa",sep="")
  library("seqinr")
  write.fasta(or16_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}

## OR13
or13_pos <- or13$V2
colnames(or13_vcf) <- or13_pos
or13_cds <- or13_vcf[,or13_pos %in% c(135623:135654, 135862:136099, 136920:137027, 137199:137680, 
                                      137904:138006, 138158:138253, 138520:138675,140086:140148)]
positions_4 <- colnames(or13_cds)
df13_cds <- data.frame(or13_cds)
colnames(df13_cds)=positions_4

fa <- readDNAStringSet("~/Downloads/Scaffold_2038.fa")
sequence = paste(fa)

## cDNA of OR13
Hs00=str_sub(sequence,1,150000)
ori_or13_cdna <- paste(subseq(Hs00,135623,135654), subseq(Hs00,135862,136099), subseq(Hs00,136920,137027), 
                       subseq(Hs00,137199,137680), subseq(Hs00,137904,138006), subseq(Hs00,138158,138253), 
                       subseq(Hs00,138520,138675), subseq(Hs00,140086,140148),sep="")
write.fasta(ori_or13_cdna, "OR13_cDNA", "~/Downloads/Hs/OR13_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR13
Hs_1=Hs00
Hs_2=Hs00
for (j in 1:20){
  for (i in positions_4){
    l=as.integer(i)
    Hs_1 <- paste(subseq(Hs_1,1,l-1),str_split(df13_cds[,i][j],"/")[[1]][1],subseq(Hs_1,l+1,str_length(Hs_1)),sep = "")
    Hs_2 <- paste(subseq(Hs_2,1,l-1),str_split(df13_cds[,i][j],"/")[[1]][2],subseq(Hs_2,l+1,str_length(Hs_2)),sep = "")
  }
  
  Hs_or13_cdna_1 <- paste(subseq(Hs_1,135623,135654), subseq(Hs_1,135862,136099), subseq(Hs_1,136920,137027), 
                          subseq(Hs_1,137199,137680), subseq(Hs_1,137904,138006), subseq(Hs_1,138158,138253), 
                          subseq(Hs_1,138520,138675), subseq(Hs_1,140086,140148),sep="")
  
  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR13_cDNA_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR13_cDNA_1.fa",sep="")
  
  write.fasta(Hs_or13_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or13_protein_1 <- as.character(translate(DNAString(Hs_or13_cdna_1),if.fuzzy.codon="X"))
  name = paste(a,"_OR13_pep_1",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR13_pep_1.fa",sep="")
  library("seqinr")
  write.fasta(or13_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  Hs_or13_cdna_2 <- paste(subseq(Hs_2,135623,135654), subseq(Hs_2,135862,136099), subseq(Hs_2,136920,137027), 
                          subseq(Hs_2,137199,137680), subseq(Hs_2,137904,138006), subseq(Hs_2,138158,138253), 
                          subseq(Hs_2,138520,138675), subseq(Hs_2,140086,140148),sep="")
  
  a=str_split(individuals[j],".sorted.")[[1]][1]
  name = paste(a,"_OR13_cDNA_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR13_cDNA_2.fa",sep="")
  
  write.fasta(Hs_or13_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or13_protein_2 <- as.character(translate(DNAString(Hs_or13_cdna_2),if.fuzzy.codon="X"))
  name = paste(a,"_OR13_pep_2",sep="")
  path = paste("~/Downloads/Hs/",a,"_OR13_pep_2.fa",sep="")
  library("seqinr")
  write.fasta(or13_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}

