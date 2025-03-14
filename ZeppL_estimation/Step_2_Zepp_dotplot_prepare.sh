#!/bin/bash
DATADIR="/scratch/mw36149/check_Chlamy"
cd $DATADIR && mkdir ZeppL_compare && cd ZeppL_compare 
# Prepare the sequences files
cat chr_name.txt | while read id; do mkdir ZeppL_"${id%%*:}"; done
cat chr_name.txt | while read id; do grep -A1 $id $DATADIR/cc1690_HiFi_chr.fa > ZeppL_"${id%%*:}"/cc1690_"$id".fa; done
cat chr_name.txt | while read id; do grep -A1 $id $DATADIR/cw15_HiFi_chr.fa > ZeppL_"${id%%*:}"/cw15_"$id".fa; done
cat chr_name.txt | while read id; do grep -A1 $id $DATADIR/cc1690_Nanopore.fa > ZeppL_"${id%%*:}"/cc1690_nanopore_"$id".fa; done
cat chr_name.txt | while read id; do grep -A1 $id $DATADIR/cc5816_Dutcher.fa > ZeppL_"${id%%*:}"/cc5816_"$id".fa; done
cat chr_name.txt | while read id; do grep -A1 $id $DATADIR/cc4532_v6.1_chr.fa > ZeppL_"${id%%*:}"/cc4532_"$id".fa; done
# *_ZeppL_v2.bed is generated in ZeppL identification step.
ml BEDTools/2.30.0-GCC-12.2.0
bedtools makewindows -b $DATADIR/cc1690_ZeppL_v2.bed -w 500 > cc1690.500b.bed
bedtools makewindows -b $DATADIR/cc5816_ZeppL_v2.bed -w 500 > cc5816.500b.bed
bedtools makewindows -b $DATADIR/cc1690_nanopore_ZeppL_v2.bed -w 500 > cc1690_nanopore.500b.bed
bedtools makewindows -b $DATADIR/cw15_ZeppL_v2.bed -w 500 > cw15.500b.bed
bedtools makewindows -b $DATADIR/cc4532_ZeppL_v2.bed -w 500 > cc4532.500b.bed
# 500bp split Zepp region
cat chr_name.txt | while read id; do grep $id cc1690.500b.bed > ZeppL_"${id%%*:}"/cc1690_ZeppL_"$id".bed; done
cat chr_name.txt | while read id; do grep $id cc5816.500b.bed > ZeppL_"${id%%*:}"/cc5816_ZeppL_"$id".bed; done
cat chr_name.txt | while read id; do grep $id cc1690_nanopore.500b.bed > ZeppL_"${id%%*:}"/cc1690_nanopore_ZeppL_"$id".bed; done
cat chr_name.txt | while read id; do grep $id cw15.500b.bed > ZeppL_"${id%%*:}"/cw15_ZeppL_"$id".bed; done
cat chr_name.txt | while read id; do grep $id cc4532.500b.bed > ZeppL_"${id%%*:}"/cc4532_ZeppL_"$id".bed; done

cat chr_name.txt | while read id; do bedtools getfasta -fi ZeppL_"${id%%*:}"/cc1690_"$id".fa -bed ZeppL_"${id%%*:}"/cc1690_ZeppL_"$id".bed -fo ZeppL_"${id%%*:}"/cc1690_ZeppL_"$id".500bp.fa; done # Seperate the secondary Zepp region as cc1690_ZeppL_chr_02T.500bp.fa
cat chr_name.txt | while read id; do bedtools getfasta -fi ZeppL_"${id%%*:}"/cc1690_nanopore_"$id".fa -bed ZeppL_"${id%%*:}"/cc1690_nanopore_ZeppL_"$id".bed -fo ZeppL_"${id%%*:}"/cc1690_nanopore_ZeppL_"$id".500bp.fa; done
cat chr_name.txt | while read id; do bedtools getfasta -fi ZeppL_"${id%%*:}"/cc5816_"$id".fa -bed ZeppL_"${id%%*:}"/cc5816_ZeppL_"$id".bed -fo ZeppL_"${id%%*:}"/cc5816_ZeppL_"$id".500bp.fa; done
cat chr_name.txt | while read id; do bedtools getfasta -fi ZeppL_"${id%%*:}"/cw15_"$id".fa -bed ZeppL_"${id%%*:}"/cw15_ZeppL_"$id".bed -fo ZeppL_"${id%%*:}"/cw15_ZeppL_"$id".500bp.fa; done
cat chr_name.txt | while read id; do bedtools getfasta -fi ZeppL_"${id%%*:}"/cc4532_"$id".fa -bed ZeppL_"${id%%*:}"/cc4532_ZeppL_"$id".bed -fo ZeppL_"${id%%*:}"/cc4532_ZeppL_"$id".500bp.fa; done

#  Whole ZeppL regions
bedtools getfasta -fi $DATADIR/cc5816_Dutcher.fa -bed $DATADIR/cc5816_ZeppL_v2.bed -fo $DATADIR/cc5816_ZeppL_total.fa
cat chr_name.txt | while read id; do grep -A1 "$id" $DATADIR/cc5816_ZeppL_total.fa > ZeppL_"$id"/cc5816_ZeppL_"$id".fa ; done
bedtools getfasta -fi $DATADIR/cc4532_v6.1_chr.fa -bed $DATADIR/cc4532_ZeppL_v2.bed -fo $DATADIR/cc4532_ZeppL_total.fa
cat chr_name.txt | while read id; do grep -A1 "$id" $DATADIR/cc4532_ZeppL_total.fa > ZeppL_"$id"/cc4532_ZeppL_"$id".fa ; done
bedtools getfasta -fi $DATADIR/cc1690_Nanopore.fa -bed $DATADIR/cc1690_nanopore_ZeppL_v2.bed -fo $DATADIR/cc1690_nanopore_ZeppL_total.fa
cat chr_name.txt | while read id; do grep -A1 "$id" $DATADIR/cc1690_nanopore_ZeppL_total.fa > ZeppL_"$id"/cc1690_nanopore_ZeppL_"$id".fa ; done
bedtools getfasta -fi $DATADIR/cc1690_HiFi_chr.fa -bed $DATADIR/cc1690_ZeppL_v2.bed -fo $DATADIR/cc1690_ZeppL_total.fa
cat chr_name.txt | while read id; do grep -A1 "$id" $DATADIR/cc1690_ZeppL_total.fa > ZeppL_"$id"/cc1690_ZeppL_"$id".fa ; done     
#chr02 has two Zepp regions. Separate the second one as cc1690_ZeppL_chr02T.fa

cd $DATADIR/ZeppL_compare
cp ../chr_name.txt ./
ml BLAST+/2.14.1-gompi-2023a
mkdir blast_result_cc1690_cc5816
mkdir blast_result_cc1690_nanopore_cc5816
mkdir blast_result_cw15_cc5816
mkdir blast_result_cc4532_cc5816

cat chr_name.txt | while read id; do makeblastdb -in ../ZeppL_"$id"/cc5816_ZeppL_"$id".fa -input_type fasta -dbtype nucl -out ../ZeppL_"$id"/cc5816_ZeppL_"$id"_blastdb; done
cat chr_name.txt | while read id; do blastn -task blastn -evalue 1e-5 -query ../ZeppL_"$id"/cc1690_ZeppL_"$id".500bp.fa -db ../ZeppL_"$id"/cc5816_ZeppL_"$id"_blastdb -out blast_result_cc1690_cc5816/cc1690_cc5816_"$id".blast.txt -outfmt 6; done
cat chr_name.txt | while read id; do blastn -task blastn -evalue 1e-5 -query ../ZeppL_"$id"/cc1690_nanopore_ZeppL_"$id".500bp.fa -db ../ZeppL_"$id"/cc5816_ZeppL_"$id"_blastdb -out blast_result_cc1690_nanopore_cc5816/cc1690_nanopore_cc5816_"$id".blast.txt -outfmt 6; done
cat chr_name.txt | while read id; do blastn -task blastn -evalue 1e-5 -query ../ZeppL_"$id"/cw15_ZeppL_"$id".500bp.fa -db ../ZeppL_"$id"/cc5816_ZeppL_"$id"_blastdb -out blast_result_cw15_cc5816/cw15_cc5816_"$id".blast.txt -outfmt 6; done
cat chr_name.txt | while read id; do blastn -task blastn -evalue 1e-5 -query ../ZeppL_"$id"/cc4532_ZeppL_"$id".500bp.fa -db ../ZeppL_"$id"/cc5816_ZeppL_"$id"_blastdb -out blast_result_cc4532_cc5816/cc4532_cc5816_"$id".blast.txt -outfmt 6; done
#cc1690 chr02 distal Zepp
cd $DATADIR/ZeppL_compare/ZeppL_chr_02
nano cc1690_chr_02_8405143_8477643.bed #cc1690: chr_02 8405143 8477643
nano cc5816_chr_02_8424698_8497198.bed #cc5816: chr_02 8424698 8497198

bedtools makewindows -b cc1690_chr_02_8405143_8477643.bed -w 500 > cc1690_chr_02_8405143_8477643.500bp.bed
bedtools getfasta -fi $DATADIR/cc1690_HiFi_chr.fa -bed cc1690_chr_02_8405143_8477643.500bp.bed -fo cc1690_chr_02_8405143_8477643.500bp.fa
bedtools getfasta -fi $DATADIR/cc5816_Dutcher.fa -bed cc5816_chr_02_8424698_8497198.bed -fo cc5816_chr_02_8424698_8497198.fa
makeblastdb -in cc5816_chr_02_8424698_8497198.fa -input_type fasta -dbtype nucl -out cc5816_chr_02_8424698_8497198_blastdb
blastn -task blastn -evalue 1e-5 -query cc1690_chr_02_8405143_8477643.500bp.fa -db cc5816_chr_02_8424698_8497198_blastdb -out cc1690_cc5816_chr_02T.blast.txt -outfmt 6
awk -F "\t" -v OFS="\t" '{ gsub(/-/, "\t", $1); gsub(/-/,"\t",$2);print }' cc1690_cc5816_chr_02T.blast.txt > cc1690_cc5816_chr_02T.blast.txt.sub
sed -i 's/:/\t/g' cc1690_cc5816_chr_02T.blast.txt.sub

cat chr_name.txt | while read id; do awk -F "\t" -v OFS="\t" '{ gsub(/-/, "\t", $1); gsub(/-/,"\t",$2);print }' blast_result_cc1690_cc5816/cc1690_cc5816_"$id".blast.txt > blast_result_cc1690_cc5816/cc1690_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do sed -i 's/:/\t/g' blast_result_cc1690_cc5816/cc1690_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do awk -F "\t" -v OFS="\t" '{ gsub(/-/, "\t", $1); gsub(/-/,"\t",$2);print }' blast_result_cc1690_nanopore_cc5816/cc1690_nanopore_cc5816_"$id".blast.txt > blast_result_cc1690_nanopore_cc5816/cc1690_nanopore_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do sed -i 's/:/\t/g' blast_result_cc1690_nanopore_cc5816/cc1690_nanopore_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do awk -F "\t" -v OFS="\t" '{ gsub(/-/, "\t", $1); gsub(/-/,"\t",$2);print }' blast_result_cw15_cc5816/cw15_cc5816_"$id".blast.txt > blast_result_cw15_cc5816/cw15_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do sed -i 's/:/\t/g' blast_result_cw15_cc5816/cw15_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do awk -F "\t" -v OFS="\t" '{ gsub(/-/, "\t", $1); gsub(/-/,"\t",$2);print }' blast_result_cc4532_cc5816/cc4532_cc5816_"$id".blast.txt > blast_result_cc4532_cc5816/cc4532_cc5816_"$id".blast.txt.sub; done
cat chr_name.txt | while read id; do sed -i 's/:/\t/g' blast_result_cc4532_cc5816/cc4532_cc5816_"$id".blast.txt.sub; done

#Generate ZeppL position files for the R plot. 
#cc1690_ZeppL_region_bar.bed
#cc4532_ZeppL_region_bar.bed
#cw15_ZeppL_region_bar.bed
#cc1690_nanopore_ZeppL_region_bar.bed

bedtools subtract -b $DATADIR/cc1690_ZeppL/ZeppL_pos_cc1690_filtered_overlap.bed -a $DATADIR/cc1690_ZeppL/cc1690_ZeppL_v2.bed > $DATADIR/cc1690_ZeppL/ZeppL_pos_cc1690_comp.bed
bedtools subtract -b $DATADIR/cc5816_ZeppL/ZeppL_pos_cc5816_filtered_overlap.bed -a $DATADIR/cc5816_ZeppL/cc5816_ZeppL_v2.bed > $DATADIR/cc5816_ZeppL/ZeppL_pos_cc5816_comp.bed
bedtools subtract -b $DATADIR/cc4532_ZeppL/ZeppL_pos_cc4532_filtered_overlap.bed -a $DATADIR/cc4532_ZeppL/cc4532_ZeppL_v2.bed > $DATADIR/cc4532_ZeppL/ZeppL_pos_cc4532_comp.bed
bedtools subtract -b $DATADIR/cw15_ZeppL/ZeppL_pos_cw15_filtered_overlap.bed -a $DATADIR/cw15_ZeppL/cw15_ZeppL_v2.bed > $DATADIR/cw15_ZeppL/ZeppL_pos_cw15_comp.bed
bedtools subtract -b $DATADIR/cc1690_nanopore_ZeppL/ZeppL_pos_cc1690_nanopore_filtered_overlap.bed -a $DATADIR/cc1690_nanopore_ZeppL/cc1690_nanopore_ZeppL_v2.bed > $DATADIR/cc1690_nanopore_ZeppL/ZeppL_pos_cc1690_nanopore_comp.bed

awk '{print$1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/cc1690_ZeppL/ZeppL_pos_cc1690_filtered_overlap.bed >> cc1690_ZeppL_region_bar.bed && awk '{print$1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/cc1690_ZeppL/ZeppL_pos_cc1690_comp.bed >> cc1690_ZeppL_region_bar.bed
awk '{print$1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/cc5816_ZeppL/ZeppL_pos_cc5816_filtered_overlap.bed >> cc5816_ZeppL_region_bar.bed && awk '{print$1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/cc5816_ZeppL/ZeppL_pos_cc5816_comp.bed >> cc5816_ZeppL_region_bar.bed
awk '{print$1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/cc4532_ZeppL/ZeppL_pos_cc4532_filtered_overlap.bed >> cc4532_ZeppL_region_bar.bed && awk '{print$1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/cc4532_ZeppL/ZeppL_pos_cc4532_comp.bed >> cc4532_ZeppL_region_bar.bed
awk '{print$1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/cw15_ZeppL/ZeppL_pos_cw15_filtered_overlap.bed >> cw15_ZeppL_region_bar.bed && awk '{print$1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/cw15_ZeppL/ZeppL_pos_cw15_comp.bed >> cw15_ZeppL_region_bar.bed
awk '{print$1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/cc1690_nanopore_ZeppL/ZeppL_pos_cc1690_nanopore_filtered_overlap.bed >> cc1690_nanopore_ZeppL_region_bar.bed && awk '{print$1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/cc1690_nanopore_ZeppL/ZeppL_pos_cc1690_nanopore_comp.bed >> cc1690_nanopore_ZeppL_region_bar.bed

awk '{print$1"\t"$3"\t"$2+310000"\t""Empty"}' $DATADIR/cc1690_ZeppL/cc1690_ZeppL_v2.bed | grep -v '8449477' >> cc1690_ZeppL_region_bar.bed
awk '{print$1"\t"$3"\t"$2+310000"\t""Empty"}' $DATADIR/cc5816_ZeppL/cc5816_ZeppL_v2.bed >> cc5816_ZeppL_region_bar.bed
awk '{print$1"\t"$3"\t"$2+310000"\t""Empty"}' $DATADIR/cc1690_nanopore_ZeppL/cc1690_nanopore_ZeppL_v2.bed >> cc1690_nanopore_ZeppL_region_bar.bed
awk '{print$1"\t"$3"\t"$2+310000"\t""Empty"}' $DATADIR/cc4532_ZeppL/cc4532_ZeppL_v2.bed >> cc4532_ZeppL_region_bar.bed
awk '{print$1"\t"$3"\t"$2+310000"\t""Empty"}' $DATADIR/cw15_ZeppL/cw15_ZeppL_v2.bed >> cw15_ZeppL_region_bar.bed