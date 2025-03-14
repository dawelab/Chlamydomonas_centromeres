
## Calculate ZeppL proportion in identified centromere regions ##
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
mkdir ZeppL_proportion_cc1690 && cd ZeppL_proportion_cc1690
grep '>' $DATADIR/cc1690_HiFi_chr.fa | sed 's/>//g' > chr_name.txt
# File: cc1690_centromere_upper.bed  cc1690_centromere_lower.bed  cc1690_ZeppL_v2.bed  ZeppL_Crei.fa
ml BEDTools/2.30.0-GCC-12.2.0
bedtools getfasta -fi $DATADIR/cc1690_HiFi_chr.fa -bed $DATADIR/cc1690_ZeppL_v2.bed -fo cc1690_ZeppL_v2.fa
bedtools getfasta -fi $DATADIR/cc1690_HiFi_chr.fa -bed $DATADIR/cc1690_centromere_upper.bed -fo cc1690_centromere_upper.fa
bedtools getfasta -fi $DATADIR/cc1690_HiFi_chr.fa -bed $DATADIR/cc1690_centromere_lower.bed -fo cc1690_centromere_lower.fa

ml BLAST+/2.13.0-gompi-2022a
makeblastdb -in cc1690_ZeppL_v2.fa -input_type fasta -dbtype nucl -out ZeppL_blastdb
makeblastdb -in cc1690_centromere_upper.fa -input_type fasta -dbtype nucl -out upper_blastdb
makeblastdb -in cc1690_centromere_lower.fa -input_type fasta -dbtype nucl -out lower_blastdb

blastn -task blastn -evalue 1e-5 -query ZeppL_Crei.fa -db ZeppL_blastdb -out ZeppL_enriched.txt -outfmt 6
blastn -task blastn -evalue 1e-5 -query ZeppL_Crei.fa -db lower_blastdb -out cen_lower.txt -outfmt 6
blastn -task blastn -evalue 1e-5 -query ZeppL_Crei.fa -db upper_blastdb -out cen_upper.txt -outfmt 6

awk '$4>300' ZeppL_enriched.txt|awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' |sort -k 1,1 -k 2,2n> ZeppL_enriched.bed && bedtools merge -i ZeppL_enriched.bed > ZeppL_enriched_merged.bed
awk '$4>300' cen_lower.txt|awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' |sort -k 1,1 -k 2,2n > cen_lower.bed && bedtools merge -i cen_lower.bed > cen_lower_merged.bed
awk '$4>300' cen_upper.txt |awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' |sort -k 1,1 -k 2,2n > cen_upper.bed && bedtools merge -i cen_upper.bed > cen_upper_merged.bed

# ZeppL proportion: sum / region_size
cat chr_name.txt|while read id; do grep "$id" ZeppL_enriched_merged.bed| awk '{sum+=$3-$2} END {print sum}';done  
cat chr_name.txt|while read id; do grep "$id" cen_lower_merged.bed| awk '{sum+=$3-$2} END {print sum}';done
cat chr_name.txt|while read id; do grep "$id" cen_upper_merged.bed| awk '{sum+=$3-$2} END {print sum}';done