setwd("/Users/mingyuwang/Desktop/2024_07Chlamy/blast_result")
library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

ZeppL_1690 <- read.table('cc1690_nanopore_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_1690) <- c("Chr","Start","End","Type")
ZeppL_5816 <- read.table('cc5816_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_5816) <- c("Chr","Start","End","Type")

chr_01 <- read.table('cc1690_nanopore_cc5816_chr_01.blast.txt.sub', sep="\t", fill=TRUE)
chr_02 <- read.table('cc1690_nanopore_cc5816_chr_02.blast.txt.sub', sep="\t", fill=TRUE)
chr_03 <- read.table('cc1690_nanopore_cc5816_chr_03.blast.txt.sub', sep="\t", fill=TRUE)
chr_04 <- read.table('cc1690_nanopore_cc5816_chr_04.blast.txt.sub', sep="\t", fill=TRUE)
chr_05 <- read.table('cc1690_nanopore_cc5816_chr_05.blast.txt.sub', sep="\t", fill=TRUE)
chr_06 <- read.table('cc1690_nanopore_cc5816_chr_06.blast.txt.sub', sep="\t", fill=TRUE)
chr_07 <- read.table('cc1690_nanopore_cc5816_chr_07.blast.txt.sub', sep="\t", fill=TRUE)
chr_08 <- read.table('cc1690_nanopore_cc5816_chr_08.blast.txt.sub', sep="\t", fill=TRUE)
chr_09 <- read.table('cc1690_nanopore_cc5816_chr_09.blast.txt.sub', sep="\t", fill=TRUE)
chr_10 <- read.table('cc1690_nanopore_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_11 <- read.table('cc1690_nanopore_cc5816_chr_11.blast.txt.sub', sep="\t", fill=TRUE)
chr_12 <- read.table('cc1690_nanopore_cc5816_chr_12.blast.txt.sub', sep="\t", fill=TRUE)
chr_13 <- read.table('cc1690_nanopore_cc5816_chr_13.blast.txt.sub', sep="\t", fill=TRUE)
chr_14 <- read.table('cc1690_nanopore_cc5816_chr_14.blast.txt.sub', sep="\t", fill=TRUE)
chr_15 <- read.table('cc1690_nanopore_cc5816_chr_15.blast.txt.sub', sep="\t", fill=TRUE)
chr_16 <- read.table('cc1690_nanopore_cc5816_chr_16.blast.txt.sub', sep="\t", fill=TRUE)
chr_17 <- read.table('cc1690_nanopore_cc5816_chr_17.blast.txt.sub', sep="\t", fill=TRUE)

data_list <- c("chr_01","chr_02","chr_03","chr_04","chr_05","chr_06","chr_07","chr_08","chr_09","chr_10","chr_11","chr_12","chr_13","chr_14","chr_15","chr_16","chr_17")
column_names <- c("qchr","qbin_start","qbin_end","schr","sbin_start","sbin_end","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
for (name in data_list) {
  df <- get(name)
  colnames(df) <- column_names
  df <- filter(df,df$bitscore>800)
  assign(name, df)
}

#chr_01
chr_01_dotplot <- ggplot() +
  geom_point(data=chr_01, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_01"),aes(xmin = min(chr_01$sbin_start)-4000*2, xmax = min(chr_01$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01$qbin_start)-4000*2, ymax =min(chr_01$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_01  ') + ylab(' UL1690.1 Chr_01  ') +
  scale_x_continuous(breaks =seq(min(chr_01$sbin_start),max(chr_01$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_01$qbin_start),max(chr_01$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_01_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_01.png", plot=chr_01_dotplot, width=10, height=8) 

#chr_02
chr_02_dotplot <- ggplot() +
  geom_point(data=filter(chr_02,chr_02$bitscore>800), aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_02"),aes(xmin = min(chr_02$sbin_start)-4000*2, xmax = min(chr_02$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_02",ZeppL_5816$Start<8400000),aes(xmin = Start, xmax = End, ymin=min(chr_02$qbin_start)-4000*2, ymax =min(chr_02$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_02  ') + ylab(' UL1690.1 Chr_02  ') +
  scale_x_continuous(breaks =seq(min(chr_02$sbin_start),max(chr_02$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_02$qbin_start),max(chr_02$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_02_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_02.png", plot=chr_02_dotplot, width=10, height=8) 

#chr_03
chr_03_dotplot <- ggplot() +
  geom_point(data=chr_03, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_03"),aes(xmin = min(chr_03$sbin_start)-4000*2, xmax = min(chr_03$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03$qbin_start)-4000*2, ymax =min(chr_03$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_03  ') + ylab(' UL1690.1 Chr_03  ') +
  scale_x_continuous(breaks =seq(min(chr_03$sbin_start),max(chr_03$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_03$qbin_start),max(chr_03$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_03_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_03.png", plot=chr_03_dotplot, width=10, height=8) 

#chr_04
chr_04_dotplot <- ggplot() +
  geom_point(data=chr_04, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_04"),aes(xmin = min(chr_04$sbin_start)-4000*2, xmax = min(chr_04$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_04"),aes(xmin = Start, xmax = End, ymin=min(chr_04$qbin_start)-4000*2, ymax =min(chr_04$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_04  ') + ylab(' UL1690.1 Chr_04  ') +
  scale_x_continuous(breaks =seq(min(chr_04$sbin_start),max(chr_04$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_04$qbin_start),max(chr_04$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_04_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_04.png", plot=chr_04_dotplot, width=10, height=8) 

#chr_05
chr_05_dotplot <- ggplot() +
  geom_point(data=chr_05, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_05"),aes(xmin = min(chr_05$sbin_start)-4000*2, xmax = min(chr_05$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_05"),aes(xmin = Start, xmax = End, ymin=min(chr_05$qbin_start)-4000*2, ymax =min(chr_05$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_05  ') + ylab(' UL1690.1 Chr_05  ') +
  scale_x_continuous(breaks =seq(min(chr_05$sbin_start),max(chr_05$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_05$qbin_start),max(chr_05$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_05_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_05.png", plot=chr_05_dotplot, width=10, height=8) 

#chr_06
chr_06_dotplot <- ggplot() +
  geom_point(data=chr_06, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_06"),aes(xmin = min(chr_06$sbin_start)-4000*2, xmax = min(chr_06$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06$qbin_start)-4000*2, ymax =min(chr_06$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_06  ') + ylab(' UL1690.1 Chr_06  ') +
  scale_x_continuous(breaks =seq(min(chr_06$sbin_start),max(chr_06$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_06$qbin_start),max(chr_06$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_06_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_06.png", plot=chr_06_dotplot, width=10, height=8) 

#chr_07
chr_07_dotplot <- ggplot() +
  geom_point(data=chr_07, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_07"),aes(xmin = min(chr_07$sbin_start)-4000*2, xmax = min(chr_07$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_07"),aes(xmin = Start, xmax = End, ymin=min(chr_07$qbin_start)-4000*2, ymax =min(chr_07$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_07  ') + ylab(' UL1690.1 Chr_07  ') +
  scale_x_continuous(breaks =seq(min(chr_07$sbin_start),max(chr_07$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_07$qbin_start),max(chr_07$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_07_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_07.png", plot=chr_07_dotplot, width=10, height=8)  

#chr_08
chr_08_dotplot <- ggplot() +
  geom_point(data=chr_08, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_08"),aes(xmin = min(chr_08$sbin_start)-4000*2, xmax = min(chr_08$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_08"),aes(xmin = Start, xmax = End, ymin=min(chr_08$qbin_start)-4000*2, ymax =min(chr_08$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_08  ') + ylab(' UL1690.1 Chr_08  ') +
  scale_x_continuous(breaks =seq(min(chr_08$sbin_start),max(chr_08$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_08$qbin_start),max(chr_08$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_08_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_08.png", plot=chr_08_dotplot, width=10, height=8) 

#chr_09
chr_09_dotplot <- ggplot() +
  geom_point(data=chr_09, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_09"),aes(xmin = min(chr_09$sbin_start)-4000*2, xmax = min(chr_09$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09$qbin_start)-4000*2, ymax =min(chr_09$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_09  ') + ylab(' UL1690.1 Chr_09  ') +
  scale_x_continuous(breaks =seq(min(chr_09$sbin_start),max(chr_09$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_09$qbin_start),max(chr_09$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_09_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_09.png", plot=chr_09_dotplot, width=10, height=8) 

#chr_10
chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_10"),aes(xmin = min(chr_10$sbin_start)-4000*2, xmax = min(chr_10$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10$qbin_start)-4000*2, ymax =min(chr_10$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_10  ') + ylab(' UL1690.1 Chr_10  ') +
  scale_x_continuous(breaks =seq(min(chr_10$sbin_start),max(chr_10$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10$qbin_start),max(chr_10$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_10_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_10.png", plot=chr_10_dotplot, width=10, height=8) 

#chr_11
chr_11_dotplot <- ggplot() +
  geom_point(data=chr_11, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_11"),aes(xmin = min(chr_11$sbin_start)-4000*2, xmax = min(chr_11$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_11"),aes(xmin = Start, xmax = End, ymin=min(chr_11$qbin_start)-4000*2, ymax =min(chr_11$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_11  ') + ylab(' UL1690.1 Chr_11  ') +
  scale_x_continuous(breaks =seq(min(chr_11$sbin_start),max(chr_11$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_11$qbin_start),max(chr_11$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_11_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_11.png", plot=chr_11_dotplot, width=10, height=8) 

#chr_12
chr_12_dotplot <- ggplot() +
  geom_point(data=chr_12, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_12"),aes(xmin = min(chr_12$sbin_start)-4000*2, xmax = min(chr_12$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_12"),aes(xmin = Start, xmax = End, ymin=min(chr_12$qbin_start)-4000*2, ymax =min(chr_12$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_12  ') + ylab(' UL1690.1 Chr_12  ') +
  scale_x_continuous(breaks =seq(min(chr_12$sbin_start),max(chr_12$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_12$qbin_start),max(chr_12$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_12_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_12.png", plot=chr_12_dotplot, width=10, height=8) 

#chr_13
chr_13_dotplot <- ggplot() +
  geom_point(data=chr_13, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_13"),aes(xmin = min(chr_13$sbin_start)-4000*2, xmax = min(chr_13$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13$qbin_start)-4000*2, ymax =min(chr_13$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_13  ') + ylab(' UL1690.1 Chr_13  ') +
  scale_x_continuous(breaks =seq(min(chr_13$sbin_start),max(chr_13$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_13$qbin_start),max(chr_13$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_13_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_13.png", plot=chr_13_dotplot, width=10, height=8) 

#chr_14
chr_14_dotplot <- ggplot() +
  geom_point(data=chr_14, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_14"),aes(xmin = min(chr_14$sbin_start)-4000*2, xmax = min(chr_14$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14$qbin_start)-4000*2, ymax =min(chr_14$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_14  ') + ylab(' UL1690.1 Chr_14  ') +
  scale_x_continuous(breaks =seq(min(chr_14$sbin_start),max(chr_14$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_14$qbin_start),max(chr_14$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_14_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_14.png", plot=chr_14_dotplot, width=10, height=8) 

#chr_15
chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_15"),aes(xmin = min(chr_15$sbin_start)-4000*2, xmax = min(chr_15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15$qbin_start)-4000*2, ymax =min(chr_15$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_15  ') + ylab(' UL1690.1 Chr_15  ') +
  scale_x_continuous(breaks =seq(min(chr_15$sbin_start),max(chr_15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_15$qbin_start),max(chr_15$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_15_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_15.png", plot=chr_15_dotplot, width=10, height=8) 

#chr_16
chr_16_dotplot <- ggplot() +
  geom_point(data=chr_16, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_16"),aes(xmin = min(chr_16$sbin_start)-4000*2, xmax = min(chr_16$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16$qbin_start)-4000*2, ymax =min(chr_16$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_16  ') + ylab(' UL1690.1 Chr_16  ') +
  scale_x_continuous(breaks =seq(min(chr_16$sbin_start),max(chr_16$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_16$qbin_start),max(chr_16$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_16_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_16.png", plot=chr_16_dotplot, width=10, height=8) 

#chr_17
chr_17_dotplot <- ggplot() +
  geom_point(data=chr_17, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_17"),aes(xmin = min(chr_17$sbin_start)-4000*2, xmax = min(chr_17$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17$qbin_start)-4000*2, ymax =min(chr_17$qbin_start)-4000,fill=Type)) +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  xlab('CC5816 Chr_17  ') + ylab(' UL1690.1 Chr_17  ') +
  scale_x_continuous(breaks =seq(min(chr_17$sbin_start),max(chr_17$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_17$qbin_start),max(chr_17$qbin_end),50000))+
  coord_fixed() +
  theme()
chr_17_dotplot
ggsave(file="Figure_cc1690_nanopore_cc5816/chr_17.png", plot=chr_17_dotplot, width=10, height=8) 

scale <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
p_test <- plot_grid(chr_01_dotplot,chr_02_dotplot,chr_02T_dotplot,chr_03_dotplot,chr_04_dotplot,chr_05_dotplot,chr_06_dotplot,chr_07_dotplot,chr_08_dotplot,chr_09_dotplot,chr_10_dotplot,chr_11_dotplot,chr_12_dotplot,chr_13_dotplot,chr_14_dotplot,chr_15_dotplot,chr_16_dotplot,chr_17_dotplot,align="hv",ncol=3,nrow=6,scale = scale)
ggsave(file="Figure_cc1690_nanopore_cc5816/cc1690_cc5816_dotplot.png", plot=p_test, width=50, height=100,limitsize = FALSE) 
ggsave(file="Figure_cc1690_nanopore_cc5816/cc1690_cc5816_dotplot.pdf", plot=p_test, width=50, height=100,limitsize = FALSE)