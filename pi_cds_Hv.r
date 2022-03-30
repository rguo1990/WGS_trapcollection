## the script was used to calculate pi values for the whole OR genes and their cds

## read the site-pi data made by vcftools
or6 <- read.csv("/Users/rongguo/Downloads/Pi-OR6.sites.pi",sep = "\t")
or14 <- read.csv("/Users/rongguo/Downloads/Pi-OR14.sites.pi",sep = "\t")
or15 <- read.csv("/Users/rongguo/Downloads/Pi-OR15.sites.pi",sep = "\t")
or16 <- read.csv("/Users/rongguo/Downloads/Pi-OR16.sites.pi",sep = "\t")
or11 <- read.csv("/Users/rongguo/Downloads/Pi-OR11.sites.pi",sep = "\t")
or13 <- read.csv("/Users/rongguo/Downloads/Pi-OR13.sites.pi",sep = "\t")
orco <- read.csv("/Users/rongguo/Downloads/Pi-orco.sites.pi",sep = "\t")

## calculate pi-values of the whole genes
pi_or6 <- round(mean(or6$PI),3)
pi_or14 <- round(mean(or14$PI),3)
pi_or15 <- round(mean(or15$PI),3)
pi_or16 <- round(mean(or16$PI),3)
pi_or11 <- round(mean(or11$PI),3)
pi_or13 <- round(mean(or13$PI),3)
pi_orco <- round(mean(orco$PI),3)
sd_or6 <- round(sd(or6$PI),3)
sd_or14 <- round(sd(or14$PI),3)
sd_or15 <- round(sd(or15$PI),3)
sd_or16 <- round(sd(or16$PI),3)
sd_or11 <- round(sd(or11$PI),3)
sd_or13 <- round(sd(or13$PI),3)
sd_orco <- round(sd(orco$PI),3)

## calculate pi-values for CDS
or6_cds <- or6[or6$POS %in% 26260:26300 | or6$POS %in% 26725:26983 | or6$POS %in% 27443:27553 | 
      or6$POS %in% 28260:28492 | or6$POS %in% 28790:29026 | or6$POS %in% 29550:29652 | 
      or6$POS %in% 29719:29814 | or6$POS %in% 30373:30528,]

pi_or6_cds <- round(mean(or6_cds$PI),3)
sd_or6_cds <- round(sd(or6_cds$PI),3)

or14_cds <- or14[or14$POS %in% 83472:83515 | or14$POS %in% 82882:83158 | or14$POS %in% 82202:82312 |
       or14$POS %in% 81610:81842 |or14$POS %in% 80866:81120 | or14$POS %in% 80475:80577 | 
       or14$POS %in% 79723:79818 | or14$POS %in% 79566:79571,]
pi_or14_cds <- round(mean(or14_cds$PI),3)
sd_or14_cds <- round(sd(or14_cds$PI),3)

or15_cds <- or15[or15$POS %in% 68317:68360 | or15$POS %in% 69576:69852 | or15$POS %in% 70225:70335 | or15$POS %in% 71006:71238 |
  or15$POS %in% 74624:74875 | or15$POS %in% 75389:75491 | or15$POS %in% 76003:76098 | or15$POS %in% 76316:76471 |
  or15$POS %in% 77092:77139,]
pi_or15_cds <- round(mean(or15_cds$PI),3)
sd_or15_cds <- round(sd(or15_cds$PI),3)

or16_cds <- or16[or16$POS %in% 58367:58407 | or16$POS %in% 58691:58949 | or16$POS %in% 59928:60038 | or16$POS %in% 61380:61612 |
  or16$POS %in% 62019:62240 | or16$POS %in% 62720:62821 | or16$POS %in% 62894:62894,]
pi_or16_cds <- round(mean(or16_cds$PI),3)
sd_or16_cds <- round(sd(or16_cds$PI),3)

orco_cds <- orco[orco$POS %in% 29259:29358 | orco$POS %in% 28112:28322 | orco$POS %in% 21276:21412 | orco$POS %in% 19216:19376 |
  orco$POS %in% 18769:18955 | orco$POS %in% 18130:18340 | orco$POS %in% 16132:16231 | orco$POS %in% 14799:14903 |
  orco$POS %in% 14241:14396 | orco$POS %in% 13665:13715,]
pi_orco_cds <- round(mean(orco_cds$PI),3)
sd_orco_cds <- round(sd(orco_cds$PI),3)

or11_cds <- or11[or11$POS %in% 49235:49266 | or11$POS %in% 48784:49021 | or11$POS %in% 48251:48358 | or11$POS %in% 47610:48094 |
  or11$POS %in% 47257:47358 | or11$POS %in% 46203:46242,]
pi_or11_cds <- round(mean(or11_cds$PI),3)
sd_or11_cds <- round(sd(or11_cds$PI),3)

or13_cds <- or13[or13$POS %in% 97:709 | or13$POS %in% 2387:2582 | or13$POS %in% 2951:3165 | or13$POS %in% 3370:4097,]
pi_or13_cds <- round(mean(or13_cds$PI),3)
sd_or13_cds <- round(sd(or13_cds$PI),3)
