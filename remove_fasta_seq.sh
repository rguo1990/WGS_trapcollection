## the script was used to remove multi fasta sequences from the original genome assembly
## headers.txt contains ">" and the headers

awk '(NR==FNR) { toRemove[$1]; next } /^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 } p' headers.txt file.fasta
