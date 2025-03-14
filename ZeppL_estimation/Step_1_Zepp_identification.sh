#!/bin/bash
####== ZeppL identification ==####
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
mkdir -p $DATADIR/ZeppL_identify && cd $DATADIR/ZeppL_identify

# Load required modules
ml BLAST+/2.13.0-gompi-2022a 
ml BEDTools/2.30.0-GCC-12.2.0

# Define genomes
GENOMES=("cc1690_HiFi_chr.fa" "cc1690_Nanopore.fa" "cw15_HiFi_chr.fa" "cc5816_Dutcher.fa" "cc4532_v6.1_chr.fa")
PREFIXES=("cc1690" "cc1690_nanopore" "cw15" "cc5816" "cc4532")

# Function for BLAST search
run_blast() {
    local genome=$1
    local prefix=$2
    mkdir -p $DATADIR/ZeppL_identify/${prefix}_ZeppL
    cd $DATADIR/ZeppL_identify/${prefix}_ZeppL

    makeblastdb -in $DATADIR/$genome -input_type fasta -dbtype nucl -out ${prefix}_blastdb
    blastn -task blastn -evalue 1e-5 -query $DATADIR/ZeppL_Crei.fa -db ${prefix}_blastdb -out ZeppL_pos_${prefix}.txt -outfmt 6
    awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' ZeppL_pos_${prefix}.txt > ZeppL_pos_${prefix}.bed
    cp ZeppL_pos_${prefix}.bed $DATADIR/ZeppL_pos_${prefix}.bed
}

# Run BLAST for all genomes
for i in ${!GENOMES[@]}; do
    run_blast "${GENOMES[$i]}" "${PREFIXES[$i]}"
done

# Function to filter BLAST results
filter_blast_hits() {
    local prefix=$1
    awk '$3>90 && $4>=300' $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}.txt |
        awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' | sort -k1,2 > $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_filtered.bed
}

# Filter BLAST hits for all genomes
for prefix in "${PREFIXES[@]}"; do
    filter_blast_hits "$prefix"
done

# Function to merge and extract ZeppL-enriched regions
merge_zeppL_regions() {
    local prefix=$1
    bedtools merge -d 40000 -i $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_filtered.bed |
        awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' > $DATADIR/ZeppL_identify/${prefix}_ZeppL/merge_ZeppL_${prefix}_40kb.txt

    awk '$4>8000' $DATADIR/ZeppL_identify/${prefix}_ZeppL/merge_ZeppL_${prefix}_40kb.txt > $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_enriched_${prefix}_40kb.bed
    awk '{print$1"\t"$2"\t"$3}' $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_enriched_${prefix}_40kb.bed > $DATADIR/ZeppL_identify/${prefix}_ZeppL/${prefix}_ZeppL_v1.bed
}

# Process ZeppL-enriched regions for all genomes
for prefix in "${PREFIXES[@]}"; do
    merge_zeppL_regions "$prefix"
done

# Manually pick the largest region as ZeppL-enriched region
# Generate final files: *_ZeppL_v2.bed (manual step needed)
# copy *_ZeppL_v2.bed to $DATADIR

# Get ZeppL coordinates overlapping enriched regions
extract_zeppL_8kb() {
    local prefix=$1
    bedtools intersect -wa -a $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_filtered.bed \
        -b $DATADIR/ZeppL_identify/${prefix}_ZeppL/${prefix}_ZeppL_v2.bed > $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_filtered_overlap.bed

    awk '$3-$2>=8000' $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_filtered_overlap.bed > $DATADIR/ZeppL_identify/${prefix}_ZeppL/ZeppL_pos_${prefix}_8kb.bed
}

# Extract ZeppL > 8kb for all genomes
for prefix in "${PREFIXES[@]}"; do
    extract_zeppL_8kb "$prefix"
done

# Extract ZeppL > 8kb sequences
mkdir -p $DATADIR/ZeppL_8kb
for i in ${!GENOMES[@]}; do
    bedtools getfasta -fi $DATADIR/${GENOMES[$i]} -bed $DATADIR/ZeppL_identify/${PREFIXES[$i]}_ZeppL/ZeppL_pos_${PREFIXES[$i]}_8kb.bed \
        -fo $DATADIR/ZeppL_8kb/${PREFIXES[$i]}_8kb_ZeppL.fa
done

# Compute ZeppL density for Figure 3
mkdir -p $DATADIR/ZeppL_density && cd $DATADIR/ZeppL_density
ml bioawk/1.0-GCC-11.2.0

# Generate genome size file
bioawk -c fastx '{print $name,length($seq)}' $DATADIR/cc1690_HiFi_chr.fa > cc1690_HiFi_chr_size.bed

# Compute ZeppL coverage in 5kb windows
ml BEDTools/2.30.0-GCC-12.2.0
bedtools makewindows -g cc1690_HiFi_chr_size.bed -w 5000 > cc1690.fasta.5kb.bed
bedtools coverage -a cc1690.fasta.5kb.bed -b $DATADIR/ZeppL_pos_cc1690.bed > ZeppL.coverage.5kb.bed