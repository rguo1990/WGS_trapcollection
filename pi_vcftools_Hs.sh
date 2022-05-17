## the script was used to calculate site-pi values for OR genes in Hs population

## calculate pi values for OR genes
## OR6
vcftools --chr Scaffold_4031_rc.Scaffold_489 --from-bp 10583 --to-bp 21176 --out Pi-OR6-Hs --site-pi --vcf Filtered_Hspopulation.recode.vcf
vcftools --chr Scaffold_4031_rc.Scaffold_489 --from-bp 10583 --to-bp 21176 --out OR6-Hs --recode --vcf Filtered_Hspopulation.recode.vcf

## OR13
vcftools --chr Scaffold_2038 --from-bp 135542 --to-bp 140490 --out Pi-OR13-Hs --site-pi --vcf Filtered_Hspopulation.recode.vcf
vcftools --chr Scaffold_2038 --from-bp 135542 --to-bp 140490 --out OR13-Hs --recode --vcf Filtered_Hspopulation.recode.vcf

## OR14
vcftools --chr Scaffold_656_new --from-bp 173309 --to-bp 182119 --out Pi-OR14-Hs --site-pi --vcf Filtered_Hspopulation.recode.vcf
vcftools --chr Scaffold_656_new --from-bp 173309 --to-bp 182119 --out OR14-Hs --recode --vcf Filtered_Hspopulation.recode.vcf

## OR16
vcftools --chr Scaffold_656_new --from-bp 199741 --to-bp 206825 --out Pi-OR16-Hs --site-pi --vcf Filtered_Hspopulation.recode.vcf
vcftools --chr Scaffold_656_new --from-bp 199741 --to-bp 206825 --out OR16-Hs --recode --vcf Filtered_Hspopulation.recode.vcf

## Orco
vcftools --chr  Scaffold_8 --from-bp 38994 --to-bp 54701 --out Pi-Orco-Hs --site-pi --vcf Filtered_Hspopulation.recode.vcf
vcftools --chr  Scaffold_8 --from-bp 38994 --to-bp 54701 --out Orco-Hs --recode --vcf Filtered_Hspopulation.recode.vcf


