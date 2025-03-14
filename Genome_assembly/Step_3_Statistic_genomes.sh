#!/bin/bash
####== Information for chromosome length and gaps ==####
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
cd $DATADIR
grep -A1 'chr' cc1690_HiFi.fa > cc1690_HiFi_chr.fa
grep -A1 'chr' cw15_HiFi.fa > cw15_HiFi_chr.fa
#cc5816_Dutcher.fa is gap-free genome
#cc1690_Nanopore.fa is downloaded from NCBI JABWPN000000000
# cc4532_v6.1 is downloded from Phytozome
grep -A1 'chr' cc4532_v6.1.fa > cc4532_v6.1_chr.fa
# Chromosome length 
module load bioawk/1.0-GCC-11.2.0
bioawk -c fastx '{print $name,length($seq)}' cc1690_HiFi_chr.fa > cc1690_HiFi_chr_length.txt
bioawk -c fastx '{print $name,length($seq)}' cw15_HiFi_chr.fa > cw15_HiFi_chr_length.txt
bioawk -c fastx '{print $name,length($seq)}' cc5816_Dutcher.fa > cc5816_Dutcher_length.txt
bioawk -c fastx '{print $name,length($seq)}' cc4532_v6.1_chr.fa > cc4532_v6.1_chr_length.txt
bioawk -c fastx '{print $name,length($seq)}' cc1690_Nanopore.fa > cc1690_nanopore_length.txt

# Gap locations
./getsgap.py cc1690_HiFi_chr.fa > cc1690_HiFi_chr_gap.txt
awk '{print$1"\t"$4"\t"$5}' cc1690_HiFi_chr_gap.txt > cc1690_HiFi_gap.bed
./getsgap.py cw15_HiFi_chr.fa > cw15_HiFi_chr_gap.txt
awk '{print$1"\t"$4"\t"$5}' cw15_HiFi_chr_gap.txt > cw15_HiFi_gap.bed
./getsgap.py cc4532_v6.1_chr.fa > cc4532_v6.1_chr_gap.txt
awk '{print$1"\t"$4"\t"$5}' cc4532_v6.1_chr_gap.txt > cc4532_gap.bed
./getsgap.py cc1690_Nanopore.fa > cc1690_nanopore_gap.txt
awk '{print$1"\t"$4"\t"$5}' cc1690_nanopore_gap.txt > cc1690_nanopore_gap.bed