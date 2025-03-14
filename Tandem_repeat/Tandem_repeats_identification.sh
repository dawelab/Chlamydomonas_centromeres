#!/bin/bash
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
mkdir $DATADIR/Tandem_repeats && cd $DATADIR/Tandem_repeats
ml TRF/4.09.1-GCCcore-11.3.0
trf $DATADIR/cc1690_HiFi_chr.fa 2 7 7 80 10 50 2000 -d -h
awk '/^Sequence:/ {chrom = substr($2, 1); next}  /^[0-9]/ {print chrom, $0}' OFS="\t" cc1690_HiFi_chr.fa.2.7.7.80.10.50.2000.dat > cc1690_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed
sed -i 's/ /\t/g' cc1690_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed

trf $DATADIR/cw15_HiFi_chr.fa 2 7 7 80 10 50 2000 -d -h
awk '/^Sequence:/ {chrom = substr($2, 1); next}  /^[0-9]/ {print chrom, $0}' OFS="\t" cw15_HiFi_chr.fa.2.7.7.80.10.50.2000.dat > cw15_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed
sed -i 's/ /\t/g' cw15_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed

trf $DATADIR/cc5816_Dutcher.fa 2 7 7 80 10 50 2000 -d -h
awk '/^Sequence:/ {chrom = substr($2, 1); next}  /^[0-9]/ {print chrom, $0}' OFS="\t" cc5816_Dutcher.fa.2.7.7.80.10.50.2000.dat > cc5816_Dutcher.fa.2.7.7.80.10.50.2000.dat.bed
sed -i 's/ /\t/g' cc5816_Dutcher.fa.2.7.7.80.10.50.2000.dat.bed

trf $DATADIR/cc1690_nanopore.fa 2 7 7 80 10 50 2000 -d -h
awk '/^Sequence:/ {chrom = substr($2, 1); next}  /^[0-9]/ {print chrom, $0}' OFS="\t" cc1690_nanopore.fa.2.7.7.80.10.50.2000.dat > cc1690_nanopore.fa.2.7.7.80.10.50.2000.dat.bed
sed -i 's/ /\t/g' cc1690_nanopore.fa.2.7.7.80.10.50.2000.dat.bed

trf $DATADIR/cc4532_v6.1.fa 2 7 7 80 10 50 2000 -d -h
awk '/^Sequence:/ {chrom = substr($2, 1); next}  /^[0-9]/ {print chrom, $0}' OFS="\t" cc4532_v6.1.fa.2.7.7.80.10.50.2000.dat > cc4532_v6.1.fa.2.7.7.80.10.50.2000.dat.bed
sed -i 's/ /\t/g' cc4532_v6.1.fa.2.7.7.80.10.50.2000.dat.bed