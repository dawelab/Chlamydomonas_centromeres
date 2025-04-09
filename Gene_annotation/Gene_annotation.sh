#!/bin/bash
#Gene annotation
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
mkdir $DATADIR/Annotation && cd $DATADIR/Annotation
# CreinhardtiiCC_4532_707_v6.1.gene.gff3 downloaded from Phytozome
# cc4532_v6.1_chr.fa only contains the chromosomes sequences
# chr_pairs.txt is required for Liftoff. It contains the chromosome pairs in the format of "chr01,chr01"
module load Liftoff/1.6.3
liftoff -exclude_partial -polish -g CreinhardtiiCC_4532_707_v6.1.gene.gff3 -o cc1690.gene.gff3 -p 28 -chroms chr_pairs.txt cc1690_HiFi_chr.fa cc4532_v6.1_chr.fa

# Only use the longest isoform to Figure3
module load AGAT/1.1.0
agat_sp_keep_longest_isoform.pl -gff cc1690.gene.gff3_polished -o cc1690.gene.gff3_polished_longest
grep -w 'gene' cc1690.gene.gff3_polished_longest| grep 'chr_' |awk '{print $1"\t"$4"\t"$5"\t"$5-$4"\t"$7}' > cc1690_gene_figure.bed
cp cc1690.gene.gff3_polished_longest $DATADIR
cp cc1690_gene_figure.bed $DATADIR
