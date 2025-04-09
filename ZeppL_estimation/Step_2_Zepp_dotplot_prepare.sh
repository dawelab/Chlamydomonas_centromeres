#!/bin/bash
DATADIR="/scratch/mw36149/Chlamydomonas/pre_dotplot"
cd "$DATADIR" && mkdir -p ZeppL_compare && cd ZeppL_compare
# Define sequence sources and filename prefixes
declare -A SOURCES=(
  ["cc1690_HiFi_chr.fa"]="cc1690"
  ["cw15_HiFi_chr.fa"]="cw15"
  ["cc1690_Nanopore.fa"]="cc1690_nanopore"
  ["cc5816_Dutcher.fa"]="cc5816"
  ["cc4532_v6.1_chr.fa"]="cc4532"
)
samples=("cc1690" "cc5816" "cc1690_nanopore" "cw15" "cc4532")

# Read chromosome names and create folders
while read -r id; do
  chr="${id%%:*}"
  mkdir -p "ZeppL_${chr}"

  for file in "${!SOURCES[@]}"; do
    prefix=${SOURCES[$file]}
    grep -A1 "$id" "$DATADIR/$file" > "ZeppL_${chr}/${prefix}_${id}.fa"
  done
done < $DATADIR/chr_name.txt

# Make 500bp windows for each genome
cd $DATADIR
ml BEDTools/2.30.0-GCC-12.2.0
for sample in "${samples[@]}"; do
  bedtools makewindows -b "$DATADIR/${sample}_ZeppL_v2.bed" -w 500 > "${sample}.500b.bed"
done

# Split Zepp regions (500bp) per chromosome
while read -r id; do
  chr="${id%%:*}"
  for sample in "${samples[@]}"; do
    grep "$id" "${sample}.500b.bed" > $DATADIR/ZeppL_compare/"ZeppL_${chr}/${sample}_ZeppL_${id}.bed"
  done
done < $DATADIR/chr_name.txt

# Extract 500bp sequences using bedtools getfasta
while read -r id; do
  chr="${id%%:*}"
  for sample in "${samples[@]}"; do
    fasta_in="ZeppL_${chr}/${sample}_${id}.fa"
    bed_in="ZeppL_${chr}/${sample}_ZeppL_${id}.bed"
    fasta_out="ZeppL_${chr}/${sample}_ZeppL_${id}.500bp.fa"
    bedtools getfasta -fi $DATADIR/ZeppL_compare/"$fasta_in" -bed $DATADIR/ZeppL_compare/"$bed_in" -fo $DATADIR/ZeppL_compare/"$fasta_out"
  done
done < $DATADIR/chr_name.txt

# Extract Whole ZeppL Regions for each genome
cd $DATADIR
bedtools getfasta -fi cc5816_Dutcher.fa -bed cc5816_ZeppL_v2.bed -fo cc5816_ZeppL_total.fa
for sample in "${samples[@]}"; do
  [[ "$sample" == "cc4532" ]] && fasta="$DATADIR/cc4532_v6.1_chr.fa"
  [[ "$sample" == "cc1690" ]] && fasta="$DATADIR/cc1690_HiFi_chr.fa"
  [[ "$sample" == "cw15" ]] && fasta="$DATADIR/cw15_HiFi_chr.fa"
  [[ "$sample" == "cc1690_nanopore" ]] && fasta="$DATADIR/cc1690_Nanopore.fa"
  bed="$DATADIR/${sample}_ZeppL_v2.bed"
  out="$DATADIR/${sample}_ZeppL_total.fa"
  bedtools getfasta -fi "$fasta" -bed "$bed" -fo "$out"
done

while read -r id; do
  chr="${id%%:*}"
  for sample in "${samples[@]}"; do
    grep -A1 "$id" "$DATADIR/${sample}_ZeppL_total.fa" > $DATADIR/ZeppL_compare/"ZeppL_${id}/${sample}_ZeppL_${id}.fa"
  done
done < $DATADIR/chr_name.txt


# Create output folders for BLAST
ml BLAST+/2.14.1-gompi-2023a
pairs=(
  "cc1690:cc5816"
  "cc1690_nanopore:cc5816"
  "cw15:cc5816"
  "cc4532:cc5816"
)

for pair in "${pairs[@]}"; do
  IFS=":" read -r q t <<< "$pair"
  mkdir -p "$DATADIR/ZeppL_compare/blast_result_${q}_${t}"
done

while read -r id; do
  chr_dir="$DATADIR/ZeppL_compare/ZeppL_${id}"
  cd "$chr_dir" || exit 1

  # Make BLAST database for cc5816 reference
  makeblastdb -in cc5816_ZeppL_"$id".fa -input_type fasta -dbtype nucl -out cc5816_ZeppL_"$id"_blastdb
  # Loop through query:target pairs
  for pair in "${pairs[@]}"; do
    IFS=":" read -r q t <<< "$pair"
    
    query_file="${q}_ZeppL_${id}.500bp.fa"
    db_prefix="cc5816_ZeppL_${id}_blastdb"
    out_file="$DATADIR/ZeppL_compare/blast_result_${q}_${t}/${q}_${t}_${id}.blast.txt"
    
    blastn -task blastn -evalue 1e-5 \
      -query "$query_file" \
      -db "$db_prefix" \
      -out "$out_file" \
      -outfmt 6
  done
done < "$DATADIR/chr_name.txt"

#### Special Analysis for distal region of chr02
cd "$DATADIR/ZeppL_compare/ZeppL_chr_02"
# Manual BED files (cc1690_chr_02_8405143_8477643.bed  cc5816_chr_02_8424698_8497198.bed) are assumed to be prepared by nano command
bedtools makewindows -b cc1690_chr_02_8405143_8477643.bed -w 500 > cc1690_chr_02T.500bp.bed
bedtools getfasta -fi "$DATADIR/cc1690_HiFi_chr.fa" -bed cc1690_chr_02T.500bp.bed -fo cc1690_chr_02T.500bp.fa
bedtools getfasta -fi "$DATADIR/cc5816_Dutcher.fa" -bed cc5816_chr_02T_8424698_8497198.bed -fo cc5816_chr_02T.fa

makeblastdb -in cc5816_chr_02T.fa -dbtype nucl -out cc5816_chr_02T_blastdb
blastn -task blastn -evalue 1e-5 -query cc1690_chr_02T.500bp.fa -db cc5816_chr_02T_blastdb -out cc1690_cc5816_chr_02T.blast.txt -outfmt 6
awk -F"\t" -v OFS="\t" '{gsub(/-/, "\t", $1); gsub(/-/, "\t", $2); print}' cc1690_cc5816_chr_02T.blast.txt | sed 's/:/\t/g' > cc1690_cc5816_chr_02T.blast.txt.sub
####

# Process the BLAST ouputs
cd $DATADIR/ZeppL_compare
output_dirs=(
  "blast_result_cc1690_cc5816"
  "blast_result_cc1690_nanopore_cc5816"
  "blast_result_cw15_cc5816"
  "blast_result_cc4532_cc5816"
)
# Clean and reformat
for dir in "${output_dirs[@]}"; do
  while read -r id; do
    file="${dir}/${dir#blast_result_}_${id}.blast.txt"
    awk -F"\t" -v OFS="\t" '{gsub(/-/, "\t", $1); gsub(/-/, "\t", $2); print}' $DATADIR/ZeppL_compare/"$file" | sed 's/:/\t/g' > $DATADIR/ZeppL_compare/"${file}.sub"
  done < $DATADIR/chr_name.txt
done


# Prepare bar plot representing Zepp elements for dotplot
cd "$DATADIR" && mkdir -p ZeppL_bar
cp $DATADIR/ZeppL_identify/*_ZeppL/ZeppL_pos_*_filtered_overlap.bed ZeppL_bar
cp $DATADIR/*ZeppL_v2.bed ZeppL_bar

# Subtract overlap regions to get complement regions
for sample in "${samples[@]}"; do
  bedtools subtract \
    -b "$DATADIR/ZeppL_bar/ZeppL_pos_${sample}_filtered_overlap.bed" \
    -a "$DATADIR/ZeppL_bar/${sample}_ZeppL_v2.bed" \
    > "$DATADIR/ZeppL_bar/ZeppL_pos_${sample}_comp.bed"
done

# Create region bar BED files: ZeppL + Non-ZeppL
for sample in "${samples[@]}"; do
  bar_file="${sample}_ZeppL_region_bar.bed"
  awk '{print $1"\t"$2"\t"$3"\t""ZeppL"}' "$DATADIR/ZeppL_bar/ZeppL_pos_${sample}_filtered_overlap.bed" > $DATADIR/ZeppL_bar/"$bar_file"
  awk '{print $1"\t"$2"\t"$3"\t""Non ZeppL"}' "$DATADIR/ZeppL_bar/ZeppL_pos_${sample}_comp.bed" >> $DATADIR/ZeppL_bar/"$bar_file"
done

# Append “Empty” regions (shift start + 310,000)
for sample in "${samples[@]}"; do
  bed="$DATADIR/ZeppL_bar/${sample}_ZeppL_v2.bed"
  bar_file="${sample}_ZeppL_region_bar.bed"

  if [[ "$sample" == "cc1690" ]]; then
    awk '{print $1"\t"$3"\t"$2+310000"\t""Empty"}' "$bed" | grep -v '8449477' >> $DATADIR/ZeppL_bar/"$bar_file"
  else
    awk '{print $1"\t"$3"\t"$2+310000"\t""Empty"}' "$bed" >> $DATADIR/ZeppL_bar/"$bar_file"
  fi
done

# cw15 chr15 ZeppL-enriched region is different with others. Some ZeppL fragments are outside of the region in cw15_ZeppL_v2.bed. Here, we generated a cw15_ZeppL_v3.bed Manually (Change: chr_15	3148368	3501191).
cp $DATADIR/ZeppL_identify/cw15_ZeppL/ZeppL_pos_cw15_filtered.bed ./
bedtools intersect -wa -a ZeppL_pos_cw15_filtered.bed \
        -b cw15_ZeppL_v3.bed > ZeppL_pos_cw15_filtered_overlap_v3.bed

bedtools subtract \
    -b $DATADIR/ZeppL_bar/ZeppL_pos_cw15_filtered_overlap_v3.bed \
    -a $DATADIR/ZeppL_bar/cw15_ZeppL_v3.bed \
    > $DATADIR/ZeppL_bar/ZeppL_pos_cw15_comp_v3.bed

awk '{print $1"\t"$2"\t"$3"\t""ZeppL"}' $DATADIR/ZeppL_bar/ZeppL_pos_cw15_filtered_overlap_v3.bed > $DATADIR/ZeppL_bar/cw15_ZeppL_region_bar_v3.bed
awk '{print $1"\t"$2"\t"$3"\t""Non ZeppL"}' $DATADIR/ZeppL_bar/ZeppL_pos_cw15_comp_v3.bed >> $DATADIR/ZeppL_bar/cw15_ZeppL_region_bar_v3.bed
awk '{print $1"\t"$3"\t"$2+310000"\t""Empty"}' cw15_ZeppL_v3.bed >> cw15_ZeppL_region_bar_v3.bed  # delete line : chr_15  3501191 3458368 Empty

