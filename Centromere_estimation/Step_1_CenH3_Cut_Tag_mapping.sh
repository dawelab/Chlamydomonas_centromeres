#!/bin/bash
# Define variables
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
SAMPLES=("Chlamy-A-CenH3-1" "Chlamy-A-CenH3-2" "Chlamy-A-CenH3-c" "Chlamy-A-IgG")

# Load required modules
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4 
module load Bowtie2/2.4.5-GCC-10.2.0 
module load SAMtools/1.16.1-GCC-11.3.0 
module load BEDTools/2.30.0-GCC-10.2.0 
module load bioawk/1.0-foss-2019b

# Step 1: Trim reads
mkdir -p "$DATADIR/Cut_Tag" && cd "$DATADIR/Cut_Tag"
for SAMPLE in "${SAMPLES[@]}"; do
    trim_galore --fastqc --gzip --paired "$DATADIR/${SAMPLE}_R1_001.fastq.gz" "$DATADIR/${SAMPLE}_R2_001.fastq.gz" -o "remove_adapter_${SAMPLE}"
done

# Step 2: Gather trimmed reads
mkdir -p total_trim_reads
find . -name '*val*.fq.gz' -exec cp {} total_trim_reads/ \;

# Step 3: Build Bowtie2 index
bowtie2-build "$DATADIR/cc1690_HiFi_chr.fa" cc1690_HiFi_chr

# Step 4: Align reads and process SAM files
for SAMPLE in "${SAMPLES[@]}"; do
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x cc1690_HiFi_chr \
        -1 total_trim_reads/"${SAMPLE}"_R1_001_val_1.fq.gz \
        -2 total_trim_reads/"${SAMPLE}"_R2_001_val_2.fq.gz \
        -S "${SAMPLE}_bowtie2.sam" &> "${SAMPLE}_bowtie2.txt"

    # Convert SAM to BAM, sort, and filter
    samtools view -bS -F 0x04 "${SAMPLE}_bowtie2.sam" | samtools sort -o "${SAMPLE}_bowtie2.sorted.bam"
    samtools view -q 20 -b "${SAMPLE}_bowtie2.sorted.bam" | samtools sort -o "${SAMPLE}_bowtie2.q20.sorted.bam"

    # Extract fragment lengths
    samtools view -F 0x04 "${SAMPLE}_bowtie2.sam" | awk '{print ($9 < 0 ? -$9 : $9)}' | sort | uniq -c | awk '{print $2, $1/2}' > "${SAMPLE}_fragmentLen.txt"

    # Convert BAM to BED
    for QUAL in "" "q20"; do
        bedtools bamtobed -bedpe -i "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.sorted.bam" > "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.bed"
        awk '$1==$4 && $6-$2 < 1000 {print $0}' "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.bed" > "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.clean.bed"
        cut -f 1,2,6 "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.clean.bed" | sort -k1,1 -k2,2n -k3,3n > "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.fragments.bed"
    done
done

# Step 5: Compute genome coverage
bioawk -c fastx '{ print $name,length($seq) }' cc1690_HiFi_chr.fa > cc1690_HiFi_chr_size.bed
for SAMPLE in "${SAMPLES[@]}"; do
    for QUAL in "" "q20"; do
        bedtools genomecov -bg -i "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.fragments.bed" -g cc1690_HiFi_chr_size.bed > "${SAMPLE}${QUAL:+_q20}.fragments.bedgraph"
    done
done

# Step 6: Compute read depth in 5kb bins
bedtools makewindows -g cc1690_HiFi_chr_size.bed -w 5000 > cc1690.fasta.5kb.bed
for SAMPLE in "${SAMPLES[@]}"; do
    for QUAL in "" "q20"; do
        bedtools coverage -a cc1690.fasta.5kb.bed -b "${SAMPLE}_bowtie2${QUAL:+.$QUAL}.fragments.bed" > "${SAMPLE}.coverage.5kb${QUAL:+_q20}.bed"
    done
done