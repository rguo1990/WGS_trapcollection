## the script was used to extract the coding sequences of OR13 and Orco genes from their vcf files, and translate them into protein sequences

library("Biostrings")
library("pegas")
library("seqinr")
library("stringr")
fa <- readDNAStringSet("~/Downloads/NWSH01000897.1.fa")
sequence_1 = paste(fa)
fa <- readDNAStringSet("~/Downloads/NWSH01002174.1.fa")
sequence_2 = paste(fa)

orco_vcf <- read.vcf("~/Downloads/Orco.recode.vcf")
or13_vcf <- read.vcf("~/Downloads/OR13.recode.vcf")

orco <- read.loci("~/Downloads/Orco.recode.vcf",header = FALSE)
or13 <- read.loci("~/Downloads/OR13.recode.vcf",header = FALSE)

## Orco
orco_pos <- orco$V2
colnames(orco_vcf) <- orco_pos
orco_cds <- orco_vcf[,orco_pos %in% c(29259:29358, 28112:28322, 21276:21412, 19216:19376, 18769:18955, 
                                      18130:18340, 16132:16231, 14799:14903, 14241:14396, 13665:13715)]
positions_1 <- colnames(orco_cds)
dforco_cds <- data.frame(orco_cds)
colnames(dforco_cds)=positions_1

## OR13
or13_pos <- or13$V2
colnames(or13_vcf) <- or13_pos
or13_cds <- or13_vcf[,or13_pos %in% c(97:709, 2387:2582, 2951:3165, 3370:4097)]
positions_2 <- colnames(or13_cds)
df13_cds <- data.frame(or13_cds)
colnames(df13_cds)=positions_2

individuals <- row.names(df13_cds)

## the reference cDNA of Orco
Hvtom1=str_sub(sequence_1,1,35000)
a=as.character(reverseComplement(DNAString(subseq(Hvtom1,29259,29358))))
b=as.character(reverseComplement(DNAString(subseq(Hvtom1,28112,28322))))
c=as.character(reverseComplement(DNAString(subseq(Hvtom1,21276,21412))))
d=as.character(reverseComplement(DNAString(subseq(Hvtom1,19216,19376))))
e=as.character(reverseComplement(DNAString(subseq(Hvtom1,18769,18955))))
f=as.character(reverseComplement(DNAString(subseq(Hvtom1,18130,18340))))
g=as.character(reverseComplement(DNAString(subseq(Hvtom1,16132,16231))))
h=as.character(reverseComplement(DNAString(subseq(Hvtom1,14799,14903))))
k=as.character(reverseComplement(DNAString(subseq(Hvtom1,14241,14396))))
m=as.character(reverseComplement(DNAString(subseq(Hvtom1,13665,13715))))

ori_orco_cdna <- paste(a,b,c,d,e,f,g,h,k,m,sep="")
write.fasta(ori_orco_cdna, "Orco_cDNA", "~/Downloads/Orco_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## the reference cDNA of OR13
Hvtom2=str_sub(sequence_2,1,5000)
ori_or13_cdna <- paste(subseq(Hvtom2,97,709),subseq(Hvtom2,2387,2582),subseq(Hvtom2,2951,3165),
                      subseq(Hvtom2,3370,4097),sep="")
write.fasta(ori_or13_cdna, "OR13_cDNA", "~/Downloads/OR13_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## Orco
Hvtom_1=Hvtom1
Hvtom_2=Hvtom1
for (j in 1:18){
  for (i in positions_1){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(dforco_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,35000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(dforco_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,35000),sep = "")
  }
  
  a=as.character(reverseComplement(DNAString(subseq(Hvtom_1,29259,29358))))
  b=as.character(reverseComplement(DNAString(subseq(Hvtom_1,28112,28322))))
  c=as.character(reverseComplement(DNAString(subseq(Hvtom_1,21276,21412))))
  d=as.character(reverseComplement(DNAString(subseq(Hvtom_1,19216,19376))))
  e=as.character(reverseComplement(DNAString(subseq(Hvtom_1,18769,18955))))
  f=as.character(reverseComplement(DNAString(subseq(Hvtom_1,18130,18340))))
  g=as.character(reverseComplement(DNAString(subseq(Hvtom_1,16132,16231))))
  h=as.character(reverseComplement(DNAString(subseq(Hvtom_1,14799,14903))))
  k=as.character(reverseComplement(DNAString(subseq(Hvtom_1,14241,14396))))
  m=as.character(reverseComplement(DNAString(subseq(Hvtom_1,13665,13715))))
  
  Hvtom_orco_cdna_1 <- paste(a,b,c,d,e,f,g,h,k,m,sep="")
  Hvtom_orco_cdna_1 <- str_replace_all(Hvtom_orco_cdna_1, "[^[:alnum:]]", "")
  
  a=str_split(individuals[j],"/")[[1]][2]
  b=str_split(a,"sorted")[[1]][1]
  name = paste(b,"Orco_cDNA_1",sep="")
  path = paste("~/Downloads/",b,"Orco_cDNA_1.fa",sep="")
  
  write.fasta(Hvtom_orco_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  orco_protein_1 <- as.character(translate(DNAString(Hvtom_orco_cdna_1)))
  name = paste(b,"Orco_pep_1",sep="")
  path = paste("~/Downloads/",b,"Orco_pep_1.fa",sep="")
  library("seqinr")
  write.fasta(orco_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  a=as.character(reverseComplement(DNAString(subseq(Hvtom_1,29259,29358))))
  b=as.character(reverseComplement(DNAString(subseq(Hvtom_1,28112,28322))))
  c=as.character(reverseComplement(DNAString(subseq(Hvtom_1,21276,21412))))
  d=as.character(reverseComplement(DNAString(subseq(Hvtom_1,19216,19376))))
  e=as.character(reverseComplement(DNAString(subseq(Hvtom_1,18769,18955))))
  f=as.character(reverseComplement(DNAString(subseq(Hvtom_1,18130,18340))))
  g=as.character(reverseComplement(DNAString(subseq(Hvtom_1,16132,16231))))
  h=as.character(reverseComplement(DNAString(subseq(Hvtom_1,14799,14903))))
  k=as.character(reverseComplement(DNAString(subseq(Hvtom_1,14241,14396))))
  m=as.character(reverseComplement(DNAString(subseq(Hvtom_1,13665,13715))))
  
  Hvtom_orco_cdna_2 <- paste(a,b,c,d,e,f,g,h,k,m,sep="")
  Hvtom_orco_cdna_2 <- str_replace_all(Hvtom_orco_cdna_2, "[^[:alnum:]]", "")
  
  a=str_split(individuals[j],"/")[[1]][2]
  b=str_split(a,"sorted")[[1]][1]
  name = paste(b,"Orco_cDNA_2",sep="")
  path = paste("~/Downloads/",b,"Orco_cDNA_2.fa",sep="")
  write.fasta(Hvtom_orco_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  orco_protein_2 <- as.character(translate(DNAString(Hvtom_orco_cdna_2)))
  name = paste(b,"Orco_pep_2",sep="")
  path = paste("~/Downloads/",b,"Orco_pep_2.fa",sep="")
  library("seqinr")
  write.fasta(orco_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}

## OR13
Hvtom_1=Hvtom2
Hvtom_2=Hvtom2
positions_2=positions_2[! positions_2 %in% c("324")]
## position 324 is an insertion
for (j in 1:18){
  for (i in positions_2){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(df13_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,5000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(df13_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,5000),sep = "")
  }

  Hvtom_or13_cdna_1 <- paste(subseq(Hvtom_1,97,709),subseq(Hvtom_1,2387,2582),subseq(Hvtom_1,2951,3165),subseq(Hvtom_1,3370,4097),sep="")
  Hvtom_or13_cdna_1 <- str_replace_all(Hvtom_or13_cdna_1, "[^[:alnum:]]", "")
  
  a=str_split(individuals[j],"/")[[1]][2]
  b=str_split(a,"sorted")[[1]][1]
  name = paste(b,"OR13_cDNA_1",sep="")
  path = paste("~/Downloads/",b,"OR13_cDNA_1.fa",sep="")
  
  write.fasta(Hvtom_or13_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or13_protein_1 <- as.character(translate(DNAString(Hvtom_or13_cdna_1),if.fuzzy.codon="X"))
  name = paste(b,"OR13_pep_1",sep="")
  path = paste("~/Downloads/",b,"OR13_pep_1.fa",sep="")
  library("seqinr")
  write.fasta(or13_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  Hvtom_or13_cdna_2 <- paste(subseq(Hvtom_2,97,709),subseq(Hvtom_2,2387,2582),subseq(Hvtom_2,2951,3165),subseq(Hvtom_2,3370,4097),sep="")
  Hvtom_or13_cdna_2 <- str_replace_all(Hvtom_or13_cdna_2, "[^[:alnum:]]", "")
  
  name = paste(b,"OR13_cDNA_2",sep="")
  path = paste("~/Downloads/",b,"OR13_cDNA_2.fa",sep="")
  write.fasta(Hvtom_or13_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
  detach("package:seqinr", unload=TRUE)
  or13_protein_2 <- as.character(translate(DNAString(Hvtom_or13_cdna_2),if.fuzzy.codon="X"))
  name = paste(b,"OR13_pep_2",sep="")
  path = paste("~/Downloads/",b,"OR13_pep_2.fa",sep="")
  library("seqinr")
  write.fasta(or13_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}

