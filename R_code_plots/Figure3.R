setwd("/Users/mingyuwang/Desktop/Chlamy_2024_10")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)
library(viridis)
library(ggsci)
library(scales)
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
options(scipen = 999)


# Process the TRF output file
name_list <- c("Chr","Start", "End", "Period", "Copies", "Consensus_Match", "Percent_Matches", "Percent_Indels", "Score", "A","C","G","T","Entropy","Repeats","Align_sequences")
cc1690_TRF <- setNames(read.table('cc1690_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
cc1690_TRF_filtered <- cc1690_TRF %>% filter(Copies>400 & Period>125)

#Process the gene bed file
bed_file <- setNames(read.table('cc1690_gene_figure.bed', sep="\t", fill=TRUE),c("Chr","Start","End","Length","Strain"))
bed_filtered <- bed_file %>% select("Chr","Start","End")
#genes <- import(bed_filtered, format = "BED")
genes <- GRanges(
  seqnames = bed_filtered$Chr,  # Chromosome column in data frame
  ranges = IRanges(start = bed_filtered$Start, end = bed_filtered$End),
  strand = "*"
)
mcols(genes) <- bed_filtered[, !(names(bed_filtered) %in% c("Chr", "Start", "End"))]
seqlengths(genes) <- c("chr_01"=8172048,"chr_02"=8574586,"chr_03"=9186790,"chr_04"=4123640,"chr_05"=3655838,"chr_06"=8952341,"chr_07"=6481895,"chr_08"=4602485,"chr_09"=6734950,"chr_10"=6719112,"chr_11"=5287831,"chr_12"=9908619,"chr_13"=5384543,"chr_14"=4156167,"chr_15"=6618926,"chr_16"=7879789,"chr_17"=6764197)
window_size <- 5000  # 5 kb
genome_windows <- tileGenome(seqlengths(genes), tilewidth = window_size, cut.last.tile.in.chrom = TRUE)
gene_density <- countOverlaps(genome_windows, genes)
genome_windows$gene_density <- gene_density
gene_density_df <- as.data.frame(genome_windows)

# process CUT&Tag bed file for three CENH3 antibodies plus IgG antibody
cenh3_1_5k_ori <- read.table('A_CenH3_1.coverage.5kb_noq20.bed', sep = "\t", fill=TRUE)
cenh3_2_5k_ori <- read.table('A_CenH3_2.coverage.5kb_noq20.bed', sep = "\t", fill=TRUE)
cenh3_c_5k_ori <- read.table('A_CenH3_c.coverage.5kb_noq20.bed', sep = "\t", fill=TRUE)
IgG_5k_ori <- read.table('A_IgG.coverage.5kb_noq20.bed', sep = "\t", fill=TRUE)

df_cenh3_1_5k <- cenh3_1_5k_ori %>% select(V1,V2,V3,V4) %>% mutate(V1,V2,V3,V4,V4/sum(V4))
df_cenh3_2_5k <- cenh3_2_5k_ori %>% select(V1,V2,V3,V4)%>% mutate(V1,V2,V3,V4,V4/sum(V4))
df_cenh3_c_5k <- cenh3_c_5k_ori %>% select(V1,V2,V3,V4)%>% mutate(V1,V2,V3,V4,V4/sum(V4))
df_IgG_5k <- IgG_5k_ori %>% select(V1,V2,V3,V4)%>% mutate(V1,V2,V3,V4,V4/sum(V4))

df_peak_1 <- data.frame(cbind(df_cenh3_1_5k$V1,df_cenh3_1_5k$V2,df_cenh3_1_5k$V3,df_cenh3_1_5k$V4,df_cenh3_1_5k$`V4/sum(V4)`,df_IgG_5k$V4,df_IgG_5k$`V4/sum(V4)`)) %>% mutate(across(2:7, as.numeric)) %>% mutate(X1,X2,X3,X4,X5,X6,X7,X5/X7)
colnames(df_peak_1)=c('chr','start_pos','end_pos','CenH3_1','CenH3_1_norm','IgG','IgG_norm','ratio')

df_peak_2 <- data.frame(cbind(df_cenh3_2_5k$V1,df_cenh3_2_5k$V2,df_cenh3_2_5k$V3,df_cenh3_2_5k$V4,df_cenh3_2_5k$`V4/sum(V4)`,df_IgG_5k$V4,df_IgG_5k$`V4/sum(V4)`)) %>% mutate(across(2:7, as.numeric)) %>% mutate(X1,X2,X3,X4,X5,X6,X7,X5/X7)
colnames(df_peak_2)=c('chr','start_pos','end_pos','CenH3_2','CenH3_2_norm','IgG','IgG_norm','ratio')

df_peak_c <- data.frame(cbind(df_cenh3_c_5k$V1,df_cenh3_c_5k$V2,df_cenh3_c_5k$V3,df_cenh3_c_5k$V4,df_cenh3_c_5k$`V4/sum(V4)`,df_IgG_5k$V4,df_IgG_5k$`V4/sum(V4)`)) %>% mutate(across(2:7, as.numeric)) %>% mutate(X1,X2,X3,X4,X5,X6,X7,X5/X7)
colnames(df_peak_c)=c('chr','start_pos','end_pos','CenH3_c','CenH3_c_norm','IgG','IgG_norm','ratio')

fi_peak_1 <- df_peak_1 %>% filter(IgG >10) %>% select('chr','start_pos','end_pos','ratio')
fi_peak_2 <- df_peak_2 %>% filter(IgG >10) %>% select('chr','start_pos','end_pos','ratio')
fi_peak_c <- df_peak_c %>% filter(IgG >10) %>% select('chr','start_pos','end_pos','ratio')


ZeppL_Crei <- read.table('ZeppL.coverage.5kb.bed', sep = "\t", fill=TRUE)
chr_length <- read.table('cc1690_co5816_chr_length.txt', sep = "\t", fill=TRUE)
cc1690gene <- read.table('cc1690_gene_figure.bed', sep = "\t", fill=TRUE)
colnames(cc1690gene)=c('seqnames','start','end','width','strand')
genes <- toGRanges(cc1690gene)

Chlamy <- toGRanges(data.frame(chr=c("chr_01", "chr_02","chr_03", "chr_04","chr_05", "chr_06","chr_07", "chr_08","chr_09", "chr_10","chr_11", "chr_12","chr_13", "chr_14","chr_15", "chr_16","chr_17"), start=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), end=c(8172048,8574586,9186790,4123640,3655838,8952341,6481895,4602485,6734950,6719112,5287831,9908619,5384543,4156167,6618926,7879789,6764197)))

#all chromosome
pp <- getDefaultPlotParams(plot.type=1) 
pp$topmargin <- 1000
pp$data1height <- 1100
pp$data1outmargin <- 4600
pp$data1inmargin <- 500
pp$ideogramheight <- 1000
png(file="Figure3_Chlamy_genome.png", width = 8500, height = 15000,pointsize = 30, bg = "white")
plot_Chlamy <- plotKaryotype(genome = Chlamy,plot.type=1,cex=8,plot.params = pp)
#kpAddBaseNumbers(plot_Chlamy,tick.dist=1000000,units = 'Mb',cex=4,tick.len = 2)
kpDataBackground(plot_Chlamy, col="white",r0=0,r1=6)
#kpPlotDensity(plot_Chlamy,data=genes,window.size = 0.005e6, col="#9A32CD",r0=0,r1=0.95,)
for(chromosome in seqlevels(plot_Chlamy$genome)) {
  peak1_df <- fi_peak_1[which(fi_peak_1$chr==chromosome),]
  peak1_x<-peak1_df$start_pos
  peak1_y<-peak1_df$ratio
  peak2_df <- fi_peak_2[which(fi_peak_2$chr==chromosome),]
  peak2_x<-peak2_df$start_pos
  peak2_y<-peak2_df$ratio
  peakc_df <- fi_peak_c[which(fi_peak_c$chr==chromosome),]
  peakc_x<-peakc_df$start_pos
  peakc_y<-peakc_df$ratio
  ZeppL_df <- ZeppL_Crei[which(ZeppL_Crei$V1==chromosome),]
  ZeppL_x<-ZeppL_df$V2
  ZeppL_y<-ZeppL_df$V4
  gen_df <- gene_density_df[which(gene_density_df$seqnames==chromosome),]
  gen_x<-gen_df$start
  gen_y<-gen_df$gene_density
  kpLines(plot_Chlamy, chr=chromosome, x=ZeppL_x, y=ZeppL_y,ymin=0,ymax=15,r0=4,r1=4.95,lwd=12,col="#534B7A")
  kpLines(plot_Chlamy, chr=chromosome, x=peakc_x, y=peakc_y,ymin=0,ymax=75,r0=3,r1=3.95,lwd=12,col="#F91151")
  kpLines(plot_Chlamy, chr=chromosome, x=peak2_x, y=peak2_y,ymin=0,ymax=75,r0=2,r1=2.95,lwd=12,col="#698B69")
  kpLines(plot_Chlamy, chr=chromosome, x=peak1_x, y=peak1_y,ymin=0,ymax=75,r0=1,r1=1.95,lwd=12,col="#B55151")
  kpArea(plot_Chlamy, chr=chromosome, x=gen_x, y=gen_y,ymin=0,ymax=6,r0=0,r1=0.95,lwd=12,col="#9A32CD")
  kpAxis(plot_Chlamy,ymin = 0, ymax=15,r0=4,r1=4.95,cex=4,side=1,numticks = 2,labels = c("",""),col="black")
  kpAxis(plot_Chlamy,ymin = 0, ymax=75,r0=3,r1=3.95,cex=4,side=1,numticks = 2,labels = c("",""),col="black")
  kpAxis(plot_Chlamy,ymin = 0, ymax=75,r0=2,r1=2.95,cex=4,side=1,numticks = 2,labels = c("",""),col="black")
  kpAxis(plot_Chlamy,ymin = 0, ymax=75,r0=1,r1=1.95,cex=4,side=1,numticks = 2,labels = c("",""),col="black")
  kpAxis(plot_Chlamy,ymin = 0, ymax=6,r0=0,r1=0.95,cex=4,side=1,numticks = 2,labels = c("",""),col="black")
  length_df <- chr_length[which(chr_length$V1==chromosome),]
  length_y <- length_df$V2
  repeats_df <- cc1690_TRF_filtered[which(cc1690_TRF_filtered$Chr==chromosome),]
  repeats_start <- repeats_df$Start
  repeats_end <- repeats_df$End
  kpRect(plot_Chlamy, chr=chromosome, x0=1, x1=length_y, y0=0, y1=1, col="#424242", data.panel="ideogram",angle = 45,border=TRUE)
  kpRect(plot_Chlamy, chr=chromosome, x0=repeats_start, x1=repeats_end, y0=0, y1=1, col="grey", data.panel="ideogram",angle = 45,border=TRUE)
}
dev.off()