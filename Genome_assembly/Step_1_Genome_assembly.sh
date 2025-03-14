#!/bin/bash
####== Genome assembly ==####
#Generate HiFi contigs for cc1690 (UL-1690) and cw15(CC-400)
#Data: HiFireads_cc1690_21gr.fastq.gz HiFireads_cw15.fastq.gz
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
module load Hifiasm/0.19.4-r575-foss-2019b

mkdir $DATADIR/cc1690_HiFi && cd $DATADIR/cc1690_HiFi
hifiasm -t 10 -l 0 --write-ec $DATADIR/HiFireads_cc1690_21gr.fastq.gz
awk '/^S/{print ">"$2;print $3}' hifiasm.asm.bp.p_ctg.gfa > $DATADIR/cc1690_HiFi_contig.fa

mkdir $DATADIR/cw15_HiFi && cd $DATADIR/cw15_HiFi
hifiasm -t 10 -l 0 --write-ec $DATADIR/HiFireads_cw15.fastq.gz
awk '/^S/{print ">"$2;print $3}' hifiasm.asm.bp.p_ctg.gfa > $DATADIR/cw15_HiFi_contig.fa

#Ragtag correction and scaffolding for HiFi genome
#Data: cc5816_Dutcher.fa as reference genome; cc1690_HiFi_contig.fa and cw15_HiFi_contig.fa as query contig file
module load RagTag/2.1.0
cd $DATADIR/RagTag_correct
ragtag.py correct -u -o cc1690_ragtag_correct $DATADIR/cc5816_Dutcher.fa $DATADIR/cc1690_HiFi_contig.fa
cp cc1690_ragtag_correct/ragtag.correct.fasta $DATADIR/cc1690_HiFi_contig_corrected.fa
ragtag.py correct -u -o cw15_ragtag_correct $DATADIR/cc5816_Dutcher.fa $DATADIR/cw15_HiFi_contig.fa
cp cw15_ragtag_correct/ragtag.correct.fasta $DATADIR/cw15_HiFi_contig_corrected.fa

cd $DATADIR/RagTag_scaffold
ragtag.py scaffold -r -t 30 -o cc1690_ragtag_scaffold $DATADIR/cc5816_Dutcher.fa $DATADIR/cc1690_HiFi_contig_corrected.fa
cp cc1690_ragtag_scaffold/ragtag.scaffold.fasta $DATADIR/cc1690_HiFi_scaffold_corrected.fa
ragtag.py scaffold -r -t 30 -o cw15_ragtag_scaffold $DATADIR/cc5816_Dutcher.fa $DATADIR/cw15_HiFi_contig_corrected.fa
cp cw15_ragtag_scaffold/ragtag.scaffold.fasta $DATADIR/cw15_HiFi_scaffold_corrected.fa
sed 's/_RagTag//g' $DATADIR/cc1690_HiFi_scaffold_corrected.fa > $DATADIR/cc1690_HiFi.fa
sed 's/_RagTag//g' $DATADIR/cw15_HiFi_scaffold_corrected.fa > $DATADIR/cw15_HiFi.fa

#Generate sequence files only containing chromosome sequences, cc1690_HiFi_chr.fa and cw15_HiFi_chr.fa
cd $DATADIR
grep -A1 'chr' cc1690_HiFi.fa > cc1690_HiFi_chr.fa
grep -A1 'chr' cw15_HiFi.fa > cw15_HiFi_chr.fa
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