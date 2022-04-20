## the script was used to calculate pi values for the whole OR genes and their cds, and plot the figures.

## read the site-pi data made by vcftools
or6 <- read.csv("/Users/rongguo/Downloads/Pi-OR6-Hs.sites.pi",sep = "\t")
or14 <- read.csv("/Users/rongguo/Downloads/Pi-OR14-Hs.sites.pi",sep = "\t")
or16 <- read.csv("/Users/rongguo/Downloads/Pi-OR16-Hs.sites.pi",sep = "\t")
or13 <- read.csv("/Users/rongguo/Downloads/Pi-OR13-Hs.sites.pi",sep = "\t")
orco <- read.csv("/Users/rongguo/Downloads/Pi-orco-Hs.sites.pi",sep = "\t")

## calculate pi-values of the whole genes
pi_or6 <- round(mean(or6$PI),3)
pi_or14 <- round(mean(or14$PI),3)
pi_or16 <- round(mean(or16$PI),3)
pi_or13 <- round(mean(or13$PI),3)
pi_orco <- round(mean(orco$PI),3)
sd_or6 <- round(sd(or6$PI),3)
sd_or14 <- round(sd(or14$PI),3)
sd_or16 <- round(sd(or16$PI),3)
sd_or13 <- round(sd(or13$PI),3)
sd_orco <- round(sd(orco$PI),3)

## calculate pi-values for CDS
or6_cds <- or6[or6$POS %in% 21136:21176 | or6$POS %in% 19954:20042 | or6$POS %in% 19633:19802 | or6$POS %in% 19056:19166 
               | or6$POS %in% 18496:18728 | or6$POS %in% 16941:17156 | or6$POS %in% 16455:16557 | or6$POS %in% 16313:16375 
               | or6$POS %in% 16145:16177 | or6$POS %in% 15655:15728 | or6$POS %in% 15449:15530 | or6$POS %in% 10583:10642,]

pi_or6_cds <- round(mean(or6_cds$PI),3)
sd_or6_cds <- round(sd(or6_cds$PI),3)

or14_cds <- or14[or14$POS %in% 173309:173352 | or14$POS %in% 173687:173963 | or14$POS %in% 174468:174578 | or14$POS %in% 174735:174967 |
                   or14$POS %in% 175197:175451 | or14$POS %in% 175758:175860 | or14$POS %in% 179832:179853 | or14$POS %in% 179991:180064 |
                   or14$POS %in% 180218:180373 | or14$POS %in% 182072:182119,]
pi_or14_cds <- round(mean(or14_cds$PI),3)
sd_or14_cds <- round(sd(or14_cds$PI),3)

or13_cds <- or13[or13$POS %in% 135623:135654 | or13$POS %in% 135862:136099 | or13$POS %in% 136920:137027 | or13$POS %in% 137199:137680 |
                   or13$POS %in% 137904:138006 | or13$POS %in% 138158:138253 | or13$POS %in% 138520:138675 | or13$POS %in% 140086:140148,]
pi_or13_cds <- round(mean(or13_cds$PI),3)
sd_or13_cds <- round(sd(or13_cds$PI),3)

or16_cds <- or16[or16$POS %in% 206203:206243 | or16$POS %in% 205652:205910 | or16$POS %in% 205127:205237 | or16$POS %in% 204102:204334 |
                   or16$POS %in% 203280:203501 | or16$POS %in% 202381:202483 | or16$POS %in% 202210:202305 | or16$POS %in% 201615:201770 |
                   or16$POS %in% 199926:199973,]
pi_or16_cds <- round(mean(or16_cds$PI),3)
sd_or16_cds <- round(sd(or16_cds$PI),3)

orco_cds <- orco[orco$POS %in% 54602:54701 | orco$POS %in% 53477:53687 | orco$POS %in% 49768:49904 | orco$POS %in% 46998:47158 |
                   orco$POS %in% 46320:46506 | orco$POS %in% 45656:45866 | orco$POS %in% 41461:41560 | orco$POS %in% 40055:40159 |
                   orco$POS %in% 39205:39360 | orco$POS %in% 38994:39044,]
pi_orco_cds <- round(mean(orco_cds$PI),3)
sd_orco_cds <- round(sd(orco_cds$PI),3)

library(ggplot2)
pi_gene <- c(pi_orco,pi_or6,pi_or13,pi_or14,pi_or16)
pi_cds <- c(pi_orco_cds,pi_or6_cds,pi_or13_cds,pi_or14_cds,pi_or16_cds)
sd_gene <- c(sd_orco,sd_or6,sd_or13,sd_or14,sd_or16)
sd_cds <- c(sd_orco_cds,sd_or6_cds,sd_or13_cds,sd_or14_cds,sd_or16_cds)
gene <- c("orco","OR6","OR13","OR14","OR16")
d <- data.frame(cbind(gene,pi_gene,pi_cds,sd_gene,sd_cds))

d$pi_gene=as.numeric(as.character(d$pi_gene))
d$pi_cds=as.numeric(as.character(d$pi_cds))

ymax=pi_cds+sd_cds
ymin=pi_cds-sd_cds  
p2 <- ggplot(d)+geom_point(aes(x=c(1:5),y=pi_cds),size=3)+
  geom_errorbar(aes( x=c(1:5),ymax = ymax, ymin = ymin),width=.1)+
  labs(x="Genes",y="PI(CDS)",title="PI values of CDS")
p2+expand_limits(y=c(0,1))+
  scale_y_continuous(breaks = seq(0,1,0.2))+scale_x_discrete(limits = gene)


