## the script was used to calculate pi values for OR genes in Hv population

## calculate pi values for OR genes (OR6, OR14-16, Orco)
vcftools --chr NWSH01000007.1  --from-bp 26260 --to-bp 30528 --out Pi-OR6 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
vcftools --chr NWSH01000007.1  --from-bp 78506 --to-bp 83663 --out Pi-OR14 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
vcftools --chr NWSH01000007.1  --from-bp 68291 --to-bp 77211 --out Pi-OR15 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
vcftools --chr NWSH01000007.1  --from-bp 58307 --to-bp 64441 --out Pi-OR16 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
vcftools --chr NWSH01000897.1  --from-bp 13533 --to-bp 34538 --out Pi-Orco --site-pi --vcf Filtered_Hvpopulation.recode.vcf
## OR11 and OR13
vcftools --chr NWSH01001866.1  --from-bp 45142 --to-bp 49299 --out Pi-OR11 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
vcftools --chr NWSH01002174.1  --from-bp 97 --to-bp 4097 --out Pi-OR13 --site-pi --vcf Filtered_Hvpopulation.recode.vcf
