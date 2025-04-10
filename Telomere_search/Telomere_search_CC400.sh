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

#set input and output directory variables
DATADIR="/scratch/mw36149/telo_chlamy"
cd $DATADIR
#Preparation HiFi reads in fasta 
ml BamTools/2.5.2-GCC-12.2.0
bamtools convert -format fastq -in *.bam -out cw15_HiFi.fastq
module load bioawk/1.0-GCC-11.2.0
bioawk -c fastx '{print ">"$name; print $seq}' cw15_HiFi.fastq > cw15_HiFi.fa

#Blast search telomere in reads
ml BLAST+/2.13.0-gompi-2022a
makeblastdb -in cw15_HiFi.fa -input_type fasta -dbtype nucl -out cw15_HiFi_reads_blastdb
blastn -task blastn-short -word_size 8 -evalue 1e-5 -query telomere_blast.fa -db cw15_HiFi_reads_blastdb -out cw15_reads_telo_blast.txt -outfmt 6
awk '$3>=96 && $4>=27' cw15_reads_telo_blast.txt | awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' | sort -k1,1 -k2,2n > cw15_telo_blast.bed

ml BEDTools/2.31.0-GCC-12.3.0
bedtools merge -d 1000 -i cw15_telo_blast.bed > cw15_telo_1000bp.bed
awk '{print$1"\t"$2"\t"$3"\t"$3-$2}' cw15_telo_1000bp.bed > cw15_telo_1000bp.bed.sub
#bioawk -c fastx '{print $name,length($seq)}' cw15_HiFi.fa > cw15_HiFi_reads.length.txt

#Extract reads with telomere
awk '{print$1}' cw15_telo_1000bp.bed | sort | uniq > cw15_telo_reads_name.txt
bioawk -cfastx 'BEGIN{while((getline k <"cw15_telo_reads_name.txt")>0)i[k]=1}{if(i[$name])print ">"$name"\n"$seq}' cw15_HiFi.fa > cw15_telo_reads.fa
bioawk -c fastx '{print $name,length($seq)}' cw15_telo_reads.fa > cw15_telo_reads_length.txt

# Trim reads and then blast
bedtools getfasta -fi cw15_telo_reads2.fa -bed telo_reads_final_HiFi.txt -fo telo_reads_filtered_trimmed.fa
ml BLAST+/2.13.0-gompi-2022a
#makeblastdb -in cw15_HiFi_chr.fa -input_type fasta -dbtype nucl -out cw15_HiFi_chr_blastdb
blastn -task blastn -evalue 1e-5 -query telo_reads_filtered_trimmed.fa -db cw15_HiFi_chr_blastdb -out cw15_telo_reads_filtered_trimmed.txt -outfmt 6
awk '$3>99 && $4>8000' cw15_telo_reads_filtered_trimmed.txt > cw15_telo_reads_filtered_trimmed.txt.sub
awk '{if($9>$10) {print $2"\t"$10"\t"$9} else {print $2"\t"$9"\t"$10}}' cw15_telo_reads_filtered_trimmed.txt.sub > cw15_telo_reads_filtered_trimmed.txt.sub.bed

