## the script was used to extract the coding sequences of OR6, OR14-16 in scaffold NWSH01000007 from their vcf files, and translate them into protein sequences

library("Biostrings")
library("pegas")
library("seqinr")
library("stringr")
fa <- readDNAStringSet("~/Downloads/NWSH01000007.fa")
sequence = paste(fa)

or6_vcf <- read.vcf("~/Downloads/OR6.recode.vcf")
or14_vcf <- read.vcf("~/Downloads/OR14.recode.vcf")
or15_vcf <- read.vcf("~/Downloads/OR15.recode.vcf")
or16_vcf <- read.vcf("~/Downloads/OR16.recode.vcf")

or6 <- read.loci("~/Downloads/OR6.recode.vcf",header = FALSE)
or14 <- read.loci("~/Downloads/OR14.recode.vcf",header = FALSE)
or15 <- read.loci("~/Downloads/OR15.recode.vcf",header = FALSE)
or16 <- read.loci("~/Downloads/OR16.recode.vcf",header = FALSE)

## OR6
or6_pos <- or6$V2
colnames(or6_vcf) <- or6_pos
or6_cds <- or6_vcf[,or6_pos %in% c(26260:26300,26725:26983,27443:27553,28260:28492,
                                              28790:29026,29550:29652,29719:29814,30373:30528)]
positions <- colnames(or6_cds)
df6_cds <- data.frame(or6_cds)
colnames(df6_cds)=positions

## OR14
or14_pos <- or14$V2
colnames(or14_vcf) <- or14_pos
or14_cds <- or14_vcf[,or14_pos %in% c(83472:83515,82882:83158,82202:82312,81610:81842,80866:81120,
                                      80475:80577,79723:79818,79566:79571)]
positions_2 <- colnames(or14_cds)
df14_cds <- data.frame(or14_cds)
colnames(df14_cds)=positions_2

## OR15
or15_pos <- or15$V2
colnames(or15_vcf) <- or15_pos
or15_cds <- or15_vcf[,or15_pos %in% c(68317:68360,69576:69852,70225:70335,71006:71238,74624:74875,
                                      75389:75491,76003:76098,76316:76471,77092:77139)]
positions_3 <- colnames(or15_cds)
df15_cds <- data.frame(or15_cds)
colnames(df15_cds)=positions_3

## OR16
or16_pos <- or16$V2
colnames(or16_vcf) <- or16_pos
or16_cds <- or16_vcf[,or16_pos %in% c(58367:58407,58691:58949,59928:60038,61380:61612,62019:62240,
                                      62720:62821,62894:62894)]
positions_4 <- colnames(or16_cds)
df16_cds <- data.frame(or16_cds)
colnames(df16_cds)=positions_4

individuals <- row.names(df6_cds)

Hvtom10=str_sub(sequence,1,100000)
ori_or6_cdna <- paste(subseq(Hvtom10,26260,26300),subseq(Hvtom10,26725,26983),subseq(Hvtom10,27443,27553),
      subseq(Hvtom10,28260,28492),subseq(Hvtom10,28790,29026),subseq(Hvtom10,29550,29652),
      subseq(Hvtom10,29719,29814),subseq(Hvtom10,30373,30528),sep="")
write.fasta(ori_or6_cdna, "OR6_cDNA", "~/Downloads/OR6_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR14 (-)
a=as.character(reverseComplement(DNAString(subseq(Hvtom10,83472,83515))))
b=as.character(reverseComplement(DNAString(subseq(Hvtom10,82882,83158))))
c=as.character(reverseComplement(DNAString(subseq(Hvtom10,82202,82312))))
d=as.character(reverseComplement(DNAString(subseq(Hvtom10,81610,81842))))
e=as.character(reverseComplement(DNAString(subseq(Hvtom10,80866,81120))))
f=as.character(reverseComplement(DNAString(subseq(Hvtom10,80475,80577))))
g=as.character(reverseComplement(DNAString(subseq(Hvtom10,79723,79818))))
h=as.character(reverseComplement(DNAString(subseq(Hvtom10,79566,79571))))

ori_or14_cdna <- paste(a,b,c,d,e,f,g,h,sep="")
write.fasta(ori_or14_cdna, "OR14_cDNA", "~/Downloads/OR14_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

ori_or15_cdna <- paste(subseq(Hvtom10,68317,68360),subseq(Hvtom10,69576,69852),subseq(Hvtom10,70225,70335),
                       subseq(Hvtom10,71006,71238),subseq(Hvtom10,74624,74875),subseq(Hvtom10,75389,75491),
                       subseq(Hvtom10,76003,76098),subseq(Hvtom10,76316,76471),subseq(Hvtom10,77092,77139),sep="")
write.fasta(ori_or15_cdna, "OR15_cDNA", "~/Downloads/OR15_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

ori_or16_cdna <- paste(subseq(Hvtom10,58367,58407),subseq(Hvtom10,58691,58949),subseq(Hvtom10,59928,60038),
                       subseq(Hvtom10,61380,61612),subseq(Hvtom10,62019,62240),subseq(Hvtom10,62720,62821),
                       subseq(Hvtom10,62894,62894),sep="")
write.fasta(ori_or16_cdna, "OR16_cDNA", "~/Downloads/OR16_cDNA.fa", open = "w", nbchar = 60, as.string = FALSE)

## OR6
Hvtom_1=Hvtom10
Hvtom_2=Hvtom10
for (j in 1:18){
  for (i in positions){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(df6_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,31000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(df6_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,31000),sep = "")
    }
    
    Hvtom_or6_cdna_1 <- paste(subseq(Hvtom_1,26260,26300),subseq(Hvtom_1,26725,26983),subseq(Hvtom_1,27443,27553),
                              subseq(Hvtom_1,28260,28492),subseq(Hvtom_1,28790,29026),subseq(Hvtom_1,29550,29652),
                              subseq(Hvtom_1,29719,29814),subseq(Hvtom_1,30373,30528),sep="")
    a=str_split(individuals[j],"/")[[1]][2]
    b=str_split(a,"sorted")[[1]][1]
    name = paste(b,"OR6_cDNA_1",sep="")
    path = paste("~/Downloads/",b,"OR6_cDNA_1.fa",sep="")
    
    write.fasta(Hvtom_or6_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or6_protein_1 <- as.character(translate(DNAString(Hvtom_or6_cdna_1),if.fuzzy.codon="X"))
    name = paste(b,"OR6_pep_1",sep="")
    path = paste("~/Downloads/",b,"OR6_pep_1.fa",sep="")
    library("seqinr")
    write.fasta(or6_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)

    Hvtom_or6_cdna_2 <- paste(subseq(Hvtom_2,26260,26300),subseq(Hvtom_2,26725,26983),subseq(Hvtom_2,27443,27553),
                              subseq(Hvtom_2,28260,28492),subseq(Hvtom_2,28790,29026),subseq(Hvtom_2,29550,29652),
                              subseq(Hvtom_2,29719,29814),subseq(Hvtom_2,30373,30528),sep="")
    name = paste(b,"OR6_cDNA_2",sep="")
    path = paste("~/Downloads/",b,"OR6_cDNA_2.fa",sep="")
    write.fasta(Hvtom_or6_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or6_protein_2 <- as.character(translate(DNAString(Hvtom_or6_cdna_2),if.fuzzy.codon="X"))
    name = paste(b,"OR6_pep_2",sep="")
    path = paste("~/Downloads/",b,"OR6_pep_2.fa",sep="")
    library("seqinr")
    write.fasta(or6_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
}
  
## OR14
Hvtom_1=Hvtom10
Hvtom_2=Hvtom10
for (j in 1:18){
  for (i in positions_2){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(df14_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,90000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(df14_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,90000),sep = "")
    }
    
    a=as.character(reverseComplement(DNAString(subseq(Hvtom_1,83472,83515))))
    b=as.character(reverseComplement(DNAString(subseq(Hvtom_1,82882,83158))))
    c=as.character(reverseComplement(DNAString(subseq(Hvtom_1,82202,82312))))
    d=as.character(reverseComplement(DNAString(subseq(Hvtom_1,81610,81842))))
    e=as.character(reverseComplement(DNAString(subseq(Hvtom_1,80866,81120))))
    f=as.character(reverseComplement(DNAString(subseq(Hvtom_1,80475,80577))))
    g=as.character(reverseComplement(DNAString(subseq(Hvtom_1,79723,79818))))
    h=as.character(reverseComplement(DNAString(subseq(Hvtom_1,79566,79571))))
    
    Hvtom_or14_cdna_1 <- paste(a,b,c,d,e,f,g,h,sep="")
    
    a=str_split(individuals[j],"/")[[1]][2]
    b=str_split(a,"sorted")[[1]][1]
    name = paste(b,"OR14_cDNA_1",sep="")
    path = paste("~/Downloads/",b,"OR14_cDNA_1.fa",sep="")
    
    write.fasta(Hvtom_or14_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or14_protein_1 <- as.character(translate(DNAString(Hvtom_or14_cdna_1),if.fuzzy.codon="X"))
    name = paste(b,"OR14_pep_1",sep="")
    path = paste("~/Downloads/",b,"OR14_pep_1.fa",sep="")
    library("seqinr")
    write.fasta(or14_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    a=as.character(reverseComplement(DNAString(subseq(Hvtom_2,83472,83515))))
    b=as.character(reverseComplement(DNAString(subseq(Hvtom_2,82882,83158))))
    c=as.character(reverseComplement(DNAString(subseq(Hvtom_2,82202,82312))))
    d=as.character(reverseComplement(DNAString(subseq(Hvtom_2,81610,81842))))
    e=as.character(reverseComplement(DNAString(subseq(Hvtom_2,80866,81120))))
    f=as.character(reverseComplement(DNAString(subseq(Hvtom_2,80475,80577))))
    g=as.character(reverseComplement(DNAString(subseq(Hvtom_2,79723,79818))))
    h=as.character(reverseComplement(DNAString(subseq(Hvtom_2,79566,79571))))
    
    Hvtom_or14_cdna_2 <- paste(a,b,c,d,e,f,g,h,sep="")
    
    a=str_split(individuals[j],"/")[[1]][2]
    b=str_split(a,"sorted")[[1]][1]
    name = paste(b,"OR14_cDNA_2",sep="")
    path = paste("~/Downloads/",b,"OR14_cDNA_2.fa",sep="")
    write.fasta(Hvtom_or14_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or14_protein_2 <- as.character(translate(DNAString(Hvtom_or14_cdna_2),if.fuzzy.codon="X"))
    name = paste(b,"OR14_pep_2",sep="")
    path = paste("~/Downloads/",b,"OR14_pep_2.fa",sep="")
    library("seqinr")
    write.fasta(or14_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)

}

## OR15
Hvtom_1=Hvtom10
Hvtom_2=Hvtom10
for (j in 1:18){
  for (i in positions_3){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(df15_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,80000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(df15_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,80000),sep = "")
    }
    
    Hvtom_or15_cdna_1 <- paste(subseq(Hvtom_1,68317,68360),subseq(Hvtom_1,69576,69852),subseq(Hvtom_1,70225,70335),
                           subseq(Hvtom_1,71006,71238),subseq(Hvtom_1,74624,74875),subseq(Hvtom_1,75389,75491),
                           subseq(Hvtom_1,76003,76098),subseq(Hvtom_1,76316,76471),subseq(Hvtom_1,77092,77139),sep="")
    Hvtom_or15_cdna_1 <- str_replace_all(Hvtom_or15_cdna_1, "[^[:alnum:]]", "")
    
    a=str_split(individuals[j],"/")[[1]][2]
    b=str_split(a,"sorted")[[1]][1]
    name = paste(b,"OR15_cDNA_1",sep="")
    path = paste("~/Downloads/",b,"OR15_cDNA_1.fa",sep="")
    
    write.fasta(Hvtom_or15_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or15_protein_1 <- as.character(translate(DNAString(Hvtom_or15_cdna_1)))
    name = paste(b,"OR15_pep_1",sep="")
    path = paste("~/Downloads/",b,"OR15_pep_1.fa",sep="")
    library("seqinr")
    write.fasta(or15_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    Hvtom_or15_cdna_2 <- paste(subseq(Hvtom_2,68317,68360),subseq(Hvtom_2,69576,69852),subseq(Hvtom_2,70225,70335),
                               subseq(Hvtom_2,71006,71238),subseq(Hvtom_2,74624,74875),subseq(Hvtom_2,75389,75491),
                               subseq(Hvtom_2,76003,76098),subseq(Hvtom_2,76316,76471),subseq(Hvtom_2,77092,77139),sep="")
    Hvtom_or15_cdna_2 <- str_replace_all(Hvtom_or15_cdna_1, "[^[:alnum:]]", "")
    
    name = paste(b,"OR15_cDNA_2",sep="")
    path = paste("~/Downloads/",b,"OR15_cDNA_2.fa",sep="")
    write.fasta(Hvtom_or15_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or15_protein_2 <- as.character(translate(DNAString(Hvtom_or15_cdna_2)))
    name = paste(b,"OR15_pep_2",sep="")
    path = paste("~/Downloads/",b,"OR15_pep_2.fa",sep="")
    library("seqinr")
    write.fasta(or15_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
  
}


## OR16
Hvtom_1=Hvtom10
Hvtom_2=Hvtom10
for (j in 1:18){
  for (i in positions_4){
    l=as.integer(i)
    Hvtom_1 <- paste(subseq(Hvtom_1,1,l-1),str_split(df16_cds[,i][j],"/")[[1]][1],subseq(Hvtom_1,l+1,65000),sep = "")
    Hvtom_2 <- paste(subseq(Hvtom_2,1,l-1),str_split(df16_cds[,i][j],"/")[[1]][2],subseq(Hvtom_2,l+1,65000),sep = "")
  }
    
    Hvtom_or16_cdna_1 <- paste(subseq(Hvtom_1,58367,58407),subseq(Hvtom_1,58691,58949),subseq(Hvtom_1,59928,60038),
                           subseq(Hvtom_1,61380,61612),subseq(Hvtom_1,62019,62240),subseq(Hvtom_1,62720,62821),
                           subseq(Hvtom_1,62894,62894),sep="")
    Hvtom_or16_cdna_1 <- str_replace_all(Hvtom_or16_cdna_1, "[^[:alnum:]]", "")
    a=str_split(individuals[j],"/")[[1]][2]
    b=str_split(a,"sorted")[[1]][1]
    name = paste(b,"OR16_cDNA_1",sep="")
    path = paste("~/Downloads/",b,"OR16_cDNA_1.fa",sep="")
    
    write.fasta(Hvtom_or16_cdna_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or16_protein_1 <- as.character(translate(DNAString(Hvtom_or16_cdna_1),if.fuzzy.codon="X"))
    name = paste(b,"OR16_pep_1",sep="")
    path = paste("~/Downloads/",b,"OR16_pep_1.fa",sep="")
    library("seqinr")
    write.fasta(or16_protein_1, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    Hvtom_or16_cdna_2 <- paste(subseq(Hvtom_2,58367,58407),subseq(Hvtom_2,58691,58949),subseq(Hvtom_2,59928,60038),
                               subseq(Hvtom_2,61380,61612),subseq(Hvtom_2,62019,62240),subseq(Hvtom_2,62720,62821),
                               subseq(Hvtom_2,62894,62894),sep="")
    Hvtom_or16_cdna_2 <- str_replace_all(Hvtom_or16_cdna_2, "[^[:alnum:]]", "")
    
    name = paste(b,"OR16_cDNA_2",sep="")
    path = paste("~/Downloads/",b,"OR16_cDNA_2.fa",sep="")
    write.fasta(Hvtom_or16_cdna_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    
    detach("package:seqinr", unload=TRUE)
    or16_protein_2 <- as.character(translate(DNAString(Hvtom_or16_cdna_2),if.fuzzy.codon="X"))
    name = paste(b,"OR16_pep_2",sep="")
    path = paste("~/Downloads/",b,"OR16_pep_2.fa",sep="")
    library("seqinr")
    write.fasta(or16_protein_2, name, path, open = "w", nbchar = 60, as.string = FALSE)
    

}

