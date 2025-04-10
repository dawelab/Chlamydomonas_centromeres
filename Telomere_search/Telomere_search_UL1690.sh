#!/bin/bash
#SBATCH --job-name=bam2fasta		                        # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=20		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=100gb			                                # Total memory for job
#SBATCH --time=12:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/mw36149/telo_chlamy/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=mw36149@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


### Part I: Identity reads with telomere sequences ###

#set input and output directory variables
DATADIR="/scratch/mw36149/telo_chlamy"
cd $DATADIR
#Preparation HiFi reads in fasta 
ml BamTools/2.5.2-GCC-12.2.0
bamtools convert -format fastq -in *.bam -out cc1690_HiFi.fastq
module load bioawk/1.0-GCC-11.2.0
bioawk -c fastx '{print ">"$name; print $seq}' cc1690_HiFi.fastq > cc1690_HiFi.fa

#Blast search telomere in reads
ml BLAST+/2.13.0-gompi-2022a
makeblastdb -in cc1690_HiFi.fa -input_type fasta -dbtype nucl -out cc1690_HiFi_reads_blastdb
blastn -num_threads 10 -task blastn-short -word_size 8 -evalue 1e-5 -query telomere_blast.fa -db cc1690_HiFi_reads_blastdb -out cc1690_reads_telo_blast.txt -outfmt 6
awk '$3>=96 && $4>=27' cc1690_reads_telo_blast.txt | awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' | sort -k1,1 -k2,2n > cc1690_telo_blast.bed

ml BEDTools/2.31.0-GCC-12.3.0
bedtools merge -d 1000 -i cc1690_telo_blast.bed > cc1690_telo_1000bp.bed
awk '{print$1"\t"$2"\t"$3"\t"$3-$2}' cc1690_telo_1000bp.bed > cc1690_telo_1000bp.bed.sub

#Extract reads with telomere
awk '{print$1}' cc1690_telo_1000bp.bed | sort | uniq > cc1690_telo_reads_name.txt
bioawk -cfastx 'BEGIN{while((getline k <"cc1690_telo_reads_name.txt")>0)i[k]=1}{if(i[$name])print ">"$name"\n"$seq}' cc1690_HiFi.fa > cc1690_telo_reads.fa
bioawk -c fastx '{print $name,length($seq)}' cc1690_telo_reads.fa > cc1690_telo_reads_length.txt


### Part II: Process telo-reads and remap the trimmed telo-reads to the genome ###

# telo_reads_final_HiFi.txt was generated using code in telo_cc1690.R
# Trim reads and then blast
bedtools getfasta -fi cc1690_telo_reads.fa -bed telo_reads_final_HiFi.txt -fo telo_reads_filtered_trimmed.fa
ml BLAST+/2.13.0-gompi-2022a
#makeblastdb -in cc1690_HiFi_chr.fa -input_type fasta -dbtype nucl -out cc1690_HiFi_chr_blastdb
blastn -num_threads 20 -task blastn -evalue 1e-5 -query telo_reads_filtered_trimmed.fa -db cc1690_HiFi_chr_blastdb -out cc1690_telo_reads_filtered_trimmed.txt -outfmt 6
awk '$3>99 && $4>8000' cc1690_telo_reads_filtered_trimmed.txt > cc1690_telo_reads_filtered_trimmed.txt.sub
awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' cc1690_telo_reads_filtered_trimmed.txt.sub > cc1690_telo_reads_filtered_trimmed.txt.sub.bed
awk '{if($9>$10) {print $2"\t"$10"\t"$9"\t"$1} else {print $2"\t"$9"\t"$10"\t"$1}}' cc1690_telo_reads_filtered_trimmed.txt.sub > cc1690_telo_reads_filtered_trimmed.txt.sub.bed2
sed -i 's/:/\t/g' cc1690_telo_reads_filtered_trimmed.txt.sub.bed2
sed -i 's/-/\t/g' cc1690_telo_reads_filtered_trimmed.txt.sub.bed2

awk '{print$4}' check_reads.txt > check_reads_name.txt 
bioawk -cfastx 'BEGIN{while((getline k <"check_reads_name.txt")>0)i[k]=1}{if(i[$name])print ">"$name"\n"$seq}' cc1690_telo_reads2.fa > check_telo_reads.fa
awk '{print$4"\t"$7"\t"$8}' check_reads.txt > check_reads_telo.bed 
bedtools getfasta -fi check_telo_reads.fa -bed check_reads_telo.bed -fo check_reads_telo.fa

bedtools getfasta -fi cc1690_telo_reads.fa -bed telo_reads_final_HiFi200.txt -fo telo_reads_filtered_trimmed_200.fa
blastn -num_threads 10 -task blastn -evalue 1e-5 -query telo_reads_filtered_trimmed_200.fa -db cc1690_HiFi_chr_blastdb -out cc1690_telo_reads_filtered_trimmed_200.txt -outfmt 6
awk '$3>99 && $4>7000' cc1690_telo_reads_filtered_trimmed_200.txt > cc1690_telo_reads_filtered_trimmed_200.txt.sub
awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' cc1690_telo_reads_filtered_trimmed_200.txt.sub > cc1690_telo_reads_filtered_trimmed_200.txt.sub.bed

