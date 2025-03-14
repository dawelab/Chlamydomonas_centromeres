setwd("/Users/mingyuwang/Desktop/2024_08Chlamy/ZeppL_8kb")
library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

ZeppL_1690 <- read.table('cc1690_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_1690) <- c("Chr","Start","End","Type")
ZeppL_1690_nanopore <- read.table('cc1690_nanopore_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_1690_nanopore) <- c("Chr","Start","End","Type")
ZeppL_cw15 <- read.table('cw15_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_cw15) <- c("Chr","Start","End","Type")
ZeppL_4532 <- read.table('cc4532_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_4532) <- c("Chr","Start","End","Type")
ZeppL_5816 <- read.table('cc5816_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_5816) <- c("Chr","Start","End","Type")

chr_01_cc1690 <- read.table('cc1690_cc5816_chr_01.blast.txt.sub', sep="\t", fill=TRUE)
chr_01_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_01.blast.txt.sub', sep="\t", fill=TRUE)
chr_01_cc4532 <- read.table('cc4532_cc5816_chr_01.blast.txt.sub', sep="\t", fill=TRUE)
chr_01_cw15 <- read.table('cw15_cc5816_chr_01.blast.txt.sub', sep="\t", fill=TRUE)

chr_03_cc1690 <- read.table('cc1690_cc5816_chr_03.blast.txt.sub', sep="\t", fill=TRUE)
chr_03_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_03.blast.txt.sub', sep="\t", fill=TRUE)
chr_03_cc4532 <- read.table('cc4532_cc5816_chr_03.blast.txt.sub', sep="\t", fill=TRUE)
chr_03_cw15 <- read.table('cw15_cc5816_chr_03.blast.txt.sub', sep="\t", fill=TRUE)

chr_06_cc1690 <- read.table('cc1690_cc5816_chr_06.blast.txt.sub', sep="\t", fill=TRUE)
chr_06_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_06.blast.txt.sub', sep="\t", fill=TRUE)
chr_06_cc4532 <- read.table('cc4532_cc5816_chr_06.blast.txt.sub', sep="\t", fill=TRUE)
chr_06_cw15 <- read.table('cw15_cc5816_chr_06.blast.txt.sub', sep="\t", fill=TRUE)

chr_09_cc1690 <- read.table('cc1690_cc5816_chr_09.blast.txt.sub', sep="\t", fill=TRUE)
chr_09_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_09.blast.txt.sub', sep="\t", fill=TRUE)
chr_09_cc4532 <- read.table('cc4532_cc5816_chr_09.blast.txt.sub', sep="\t", fill=TRUE)
chr_09_cw15 <- read.table('cw15_cc5816_chr_09.blast.txt.sub', sep="\t", fill=TRUE)

chr_10_cc1690 <- read.table('cc1690_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_10_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_10_cc4532 <- read.table('cc4532_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_10_cw15 <- read.table('cw15_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)

chr_13_cc1690 <- read.table('cc1690_cc5816_chr_13.blast.txt.sub', sep="\t", fill=TRUE)
chr_13_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_13.blast.txt.sub', sep="\t", fill=TRUE)
chr_13_cc4532 <- read.table('cc4532_cc5816_chr_13.blast.txt.sub', sep="\t", fill=TRUE)
chr_13_cw15 <- read.table('cw15_cc5816_chr_13.blast.txt.sub', sep="\t", fill=TRUE)

chr_14_cc1690 <- read.table('cc1690_cc5816_chr_14.blast.txt.sub', sep="\t", fill=TRUE)
chr_14_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_14.blast.txt.sub', sep="\t", fill=TRUE)
chr_14_cc4532 <- read.table('cc4532_cc5816_chr_14.blast.txt.sub', sep="\t", fill=TRUE)
chr_14_cw15 <- read.table('cw15_cc5816_chr_14.blast.txt.sub', sep="\t", fill=TRUE)

chr_15_cc1690 <- read.table('cc1690_cc5816_chr_15.blast.txt.sub', sep="\t", fill=TRUE)
chr_15_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_15.blast.txt.sub', sep="\t", fill=TRUE)
chr_15_cc4532 <- read.table('cc4532_cc5816_chr_15.blast.txt.sub', sep="\t", fill=TRUE)
chr_15_cw15 <- read.table('cw15_cc5816_chr_15.blast.txt.sub', sep="\t", fill=TRUE)

chr_16_cc1690 <- read.table('cc1690_cc5816_chr_16.blast.txt.sub', sep="\t", fill=TRUE)
chr_16_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_16.blast.txt.sub', sep="\t", fill=TRUE)
chr_16_cc4532 <- read.table('cc4532_cc5816_chr_16.blast.txt.sub', sep="\t", fill=TRUE)
chr_16_cw15 <- read.table('cw15_cc5816_chr_16.blast.txt.sub', sep="\t", fill=TRUE)

chr_17_cc1690 <- read.table('cc1690_cc5816_chr_17.blast.txt.sub', sep="\t", fill=TRUE)
chr_17_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_17.blast.txt.sub', sep="\t", fill=TRUE)
chr_17_cc4532 <- read.table('cc4532_cc5816_chr_17.blast.txt.sub', sep="\t", fill=TRUE)
chr_17_cw15 <- read.table('cw15_cc5816_chr_17.blast.txt.sub', sep="\t", fill=TRUE)


data_list <- c("chr_01_cc1690","chr_01_cc1690N","chr_01_cc4532","chr_01_cw15",
               "chr_03_cc1690","chr_03_cc1690N","chr_03_cc4532","chr_03_cw15",
               "chr_06_cc1690","chr_06_cc1690N","chr_06_cc4532","chr_06_cw15",
               "chr_09_cc1690","chr_09_cc1690N","chr_09_cc4532","chr_09_cw15",
               "chr_10_cc1690","chr_10_cc1690N","chr_10_cc4532","chr_10_cw15",
               "chr_13_cc1690","chr_13_cc1690N","chr_13_cc4532","chr_13_cw15",
               "chr_14_cc1690","chr_14_cc1690N","chr_14_cc4532","chr_14_cw15",
               "chr_15_cc1690","chr_15_cc1690N","chr_15_cc4532","chr_15_cw15",
               "chr_16_cc1690","chr_16_cc1690N","chr_16_cc4532","chr_16_cw15",
               "chr_17_cc1690","chr_17_cc1690N","chr_17_cc4532","chr_17_cw15")
column_names <- c("qchr","qbin_start","qbin_end","schr","sbin_start","sbin_end","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
for (name in data_list) {
  df <- get(name)
  colnames(df) <- column_names
  df <- filter(df,df$bitscore>800)
  assign(name, df)
}
#label
label <- read.table('ZeppL_coord_8kb.bed', sep="\t", fill=TRUE)
colnames(label) <- c("Strain","Chr","Start","End")
#chr_01
cc1690_cc5816_chr_01_dotplot <- ggplot() +
  geom_point(data=chr_01_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_01"),aes(xmin = min(chr_01_cc1690$sbin_start)-4000*2, xmax = min(chr_01_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc1690$qbin_start)-4000*2, ymax =min(chr_01_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_01"),aes(xmin = min(chr_01_cc1690$sbin_start)-7000*2, xmax = min(chr_01_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc1690$qbin_start)-7000*2, ymax =min(chr_01_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_01  ') + ylab(' UL1690.1 Chr_01  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_01_cc1690$sbin_start),max(chr_01_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_01_cc1690$qbin_start),max(chr_01_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_01_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_01.png", plot=cc1690_cc5816_chr_01_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_01_dotplot <- ggplot() +
  geom_point(data=chr_01_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_01"),aes(xmin = min(chr_01_cc1690N$sbin_start)-4000*2, xmax = min(chr_01_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc1690N$qbin_start)-4000*2, ymax =min(chr_01_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_01"),aes(xmin = min(chr_01_cc1690N$sbin_start)-7000*2, xmax = min(chr_01_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc1690N$qbin_start)-7000*2, ymax =min(chr_01_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_01  ') + ylab(' cc1690 Nanopore Chr_01  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_01_cc1690N$sbin_start),max(chr_01_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_01_cc1690N$qbin_start),max(chr_01_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_01_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_01.png", plot=Ncc1690_cc5816_chr_01_dotplot, width=10, height=10) 

cw15_cc5816_chr_01_dotplot <- ggplot() +
  geom_point(data=chr_01_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_01"),aes(xmin = min(chr_01_cw15$sbin_start)-4000*2, xmax = min(chr_01_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cw15$qbin_start)-4000*2, ymax =min(chr_01_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_01"),aes(xmin = min(chr_01_cw15$sbin_start)-7000*2, xmax = min(chr_01_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cw15$qbin_start)-7000*2, ymax =min(chr_01_cw15$qbin_start)-9000),fill="red") +
  geom_rect(aes(xmin = min(chr_01_cw15$sbin_start)-7000*2, xmax = min(chr_01_cw15$sbin_start)-9000, ymin=1190152, ymax = 1196948),fill="blue") + 	
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_01  ') + ylab(' CC400 Chr_01  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_01_cw15$sbin_start),max(chr_01_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_01_cw15$qbin_start),max(chr_01_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_01_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_01.png", plot=cw15_cc5816_chr_01_dotplot, width=10, height=10) 

cc4532_cc5816_chr_01_dotplot <- ggplot() +
  geom_point(data=chr_01_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_01"),aes(xmin = min(chr_01_cc4532$sbin_start)-4000*2, xmax = min(chr_01_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc4532$qbin_start)-4000*2, ymax =min(chr_01_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_01"),aes(xmin = min(chr_01_cc4532$sbin_start)-7000*2, xmax = min(chr_01_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_01"),aes(xmin = Start, xmax = End, ymin=min(chr_01_cc4532$qbin_start)-7000*2, ymax =min(chr_01_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_01  ') + ylab(' CC4532 chr_01  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_01_cc4532$sbin_start),max(chr_01_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_01_cc4532$qbin_start),max(chr_01_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_01_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_01.png", plot=cc4532_cc5816_chr_01_dotplot, width=10, height=10) 

#chr_03
cc1690_cc5816_chr_03_dotplot <- ggplot() +
  geom_point(data=chr_03_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_03"),aes(xmin = min(chr_03_cc1690$sbin_start)-4000*2, xmax = min(chr_03_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc1690$qbin_start)-4000*2, ymax =min(chr_03_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_03"),aes(xmin = min(chr_03_cc1690$sbin_start)-7000*2, xmax = min(chr_03_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc1690$qbin_start)-7000*2, ymax =min(chr_03_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_03  ') + ylab(' UL1690.1 Chr_03  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_03_cc1690$sbin_start),max(chr_03_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_03_cc1690$qbin_start),max(chr_03_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_03_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_03.png", plot=cc1690_cc5816_chr_03_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_03_dotplot <- ggplot() +
  geom_point(data=chr_03_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_03"),aes(xmin = min(chr_03_cc1690N$sbin_start)-4000*2, xmax = min(chr_03_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc1690N$qbin_start)-4000*2, ymax =min(chr_03_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_03"),aes(xmin = min(chr_03_cc1690N$sbin_start)-7000*2, xmax = min(chr_03_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc1690N$qbin_start)-7000*2, ymax =min(chr_03_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_03  ') + ylab(' cc1690 Nanopore Chr_03  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_03_cc1690N$sbin_start),max(chr_03_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_03_cc1690N$qbin_start),max(chr_03_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_03_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_03.png", plot=Ncc1690_cc5816_chr_03_dotplot, width=10, height=10) 

cw15_cc5816_chr_03_dotplot <- ggplot() +
  geom_point(data=chr_03_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_03"),aes(xmin = min(chr_03_cw15$sbin_start)-4000*2, xmax = min(chr_03_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cw15$qbin_start)-4000*2, ymax =min(chr_03_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_03"),aes(xmin = min(chr_03_cw15$sbin_start)-7000*2, xmax = min(chr_03_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cw15$qbin_start)-7000*2, ymax =min(chr_03_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_03  ') + ylab(' CC400 Chr_03  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_03_cw15$sbin_start),max(chr_03_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_03_cw15$qbin_start),max(chr_03_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_03_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_03.png", plot=cw15_cc5816_chr_03_dotplot, width=10, height=10) 

cc4532_cc5816_chr_03_dotplot <- ggplot() +
  geom_point(data=chr_03_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_03"),aes(xmin = min(chr_03_cc4532$sbin_start)-4000*2, xmax = min(chr_03_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc4532$qbin_start)-4000*2, ymax =min(chr_03_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_03"),aes(xmin = min(chr_03_cc4532$sbin_start)-7000*2, xmax = min(chr_03_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_03"),aes(xmin = Start, xmax = End, ymin=min(chr_03_cc4532$qbin_start)-7000*2, ymax =min(chr_03_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_03  ') + ylab(' CC4532 Chr_03  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_03_cc4532$sbin_start),max(chr_03_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_03_cc4532$qbin_start),max(chr_03_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_03_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_03.png", plot=cc4532_cc5816_chr_03_dotplot, width=10, height=10) 

#chr_06
cc1690_cc5816_chr_06_dotplot <- ggplot() +
  geom_point(data=chr_06_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_06"),aes(xmin = min(chr_06_cc1690$sbin_start)-4000*2, xmax = min(chr_06_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc1690$qbin_start)-4000*2, ymax =min(chr_06_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_06"),aes(xmin = min(chr_06_cc1690$sbin_start)-7000*2, xmax = min(chr_06_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc1690$qbin_start)-7000*2, ymax =min(chr_06_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_06  ') + ylab(' UL1690.1 Chr_06  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_06_cc1690$sbin_start),max(chr_06_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_06_cc1690$qbin_start),max(chr_06_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_06_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_06.png", plot=cc1690_cc5816_chr_06_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_06_dotplot <- ggplot() +
  geom_point(data=chr_06_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_06"),aes(xmin = min(chr_06_cc1690N$sbin_start)-4000*2, xmax = min(chr_06_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc1690N$qbin_start)-4000*2, ymax =min(chr_06_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_06"),aes(xmin = min(chr_06_cc1690N$sbin_start)-7000*2, xmax = min(chr_06_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc1690N$qbin_start)-7000*2, ymax =min(chr_06_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_06  ') + ylab(' cc1690 Nanopore Chr_06  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_06_cc1690N$sbin_start),max(chr_06_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_06_cc1690N$qbin_start),max(chr_06_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_06_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_06.png", plot=Ncc1690_cc5816_chr_06_dotplot, width=10, height=10) 

cw15_cc5816_chr_06_dotplot <- ggplot() +
  geom_point(data=chr_06_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_06"),aes(xmin = min(chr_06_cw15$sbin_start)-4000*2, xmax = min(chr_06_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cw15$qbin_start)-4000*2, ymax =min(chr_06_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_06"),aes(xmin = min(chr_06_cw15$sbin_start)-7000*2, xmax = min(chr_06_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cw15$qbin_start)-7000*2, ymax =min(chr_06_cw15$qbin_start)-9000),fill="red") +
  geom_rect(aes(xmin = min(chr_06_cw15$sbin_start)-7000*2, xmax = min(chr_06_cw15$sbin_start)-9000, ymin=4565244, ymax = 4572973),fill="blue") + 
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_06  ') + ylab(' CC400 Chr_06  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_06_cw15$sbin_start),max(chr_06_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_06_cw15$qbin_start),max(chr_06_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_06_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_06.png", plot=cw15_cc5816_chr_06_dotplot, width=10, height=10) 

cc4532_cc5816_chr_06_dotplot <- ggplot() +
  geom_point(data=chr_06_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_06"),aes(xmin = min(chr_06_cc4532$sbin_start)-4000*2, xmax = min(chr_06_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc4532$qbin_start)-4000*2, ymax =min(chr_06_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_06"),aes(xmin = min(chr_06_cc4532$sbin_start)-7000*2, xmax = min(chr_06_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_06"),aes(xmin = Start, xmax = End, ymin=min(chr_06_cc4532$qbin_start)-7000*2, ymax =min(chr_06_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_06  ') + ylab(' CC4532 Chr_06  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_06_cc4532$sbin_start),max(chr_06_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_06_cc4532$qbin_start),max(chr_06_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_06_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_06.png", plot=cc4532_cc5816_chr_06_dotplot, width=10, height=10) 

#chr_09
cc1690_cc5816_chr_09_dotplot <- ggplot() +
  geom_point(data=chr_09_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_09"),aes(xmin = min(chr_09_cc1690$sbin_start)-4000*2, xmax = min(chr_09_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc1690$qbin_start)-4000*2, ymax =min(chr_09_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_09"),aes(xmin = min(chr_09_cc1690$sbin_start)-7000*2, xmax = min(chr_09_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc1690$qbin_start)-7000*2, ymax =min(chr_09_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_rect(aes(xmin = min(chr_09_cc1690$sbin_start)-7000*2, xmax = min(chr_09_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_09  ') + ylab(' UL1690.1 Chr_09  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_09_cc1690$sbin_start),max(chr_09_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_09_cc1690$qbin_start),max(chr_09_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_09_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_09.png", plot=cc1690_cc5816_chr_09_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_09_dotplot <- ggplot() +
  geom_point(data=chr_09_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_09"),aes(xmin = min(chr_09_cc1690N$sbin_start)-4000*2, xmax = min(chr_09_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc1690N$qbin_start)-4000*2, ymax =min(chr_09_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_09"),aes(xmin = min(chr_09_cc1690N$sbin_start)-7000*2, xmax = min(chr_09_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc1690N$qbin_start)-7000*2, ymax =min(chr_09_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_09  ') + ylab(' cc1690 Nanopore Chr_09  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_09_cc1690N$sbin_start),max(chr_09_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_09_cc1690N$qbin_start),max(chr_09_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_09_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_09.png", plot=Ncc1690_cc5816_chr_09_dotplot, width=10, height=10) 

cw15_cc5816_chr_09_dotplot <- ggplot() +
  geom_point(data=chr_09_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_09"),aes(xmin = min(chr_09_cw15$sbin_start)-4000*2, xmax = min(chr_09_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cw15$qbin_start)-4000*2, ymax =min(chr_09_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_09"),aes(xmin = min(chr_09_cw15$sbin_start)-7000*2, xmax = min(chr_09_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cw15$qbin_start)-7000*2, ymax =min(chr_09_cw15$qbin_start)-9000),fill="red") +
  geom_rect(aes(xmin = min(chr_09_cw15$sbin_start)-7000*2, xmax = min(chr_09_cw15$sbin_start)-9000, ymin=4394981, ymax = 4389091),fill="blue") + 	
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_09  ') + ylab(' CC400 Chr_09  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_09_cw15$sbin_start),max(chr_09_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_09_cw15$qbin_start),max(chr_09_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_09_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_09.png", plot=cw15_cc5816_chr_09_dotplot, width=10, height=10) 

cc4532_cc5816_chr_09_dotplot <- ggplot() +
  geom_point(data=chr_09_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_09"),aes(xmin = min(chr_09_cc4532$sbin_start)-4000*2, xmax = min(chr_09_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc4532$qbin_start)-4000*2, ymax =min(chr_09_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_09"),aes(xmin = min(chr_09_cc4532$sbin_start)-7000*2, xmax = min(chr_09_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_09"),aes(xmin = Start, xmax = End, ymin=min(chr_09_cc4532$qbin_start)-7000*2, ymax =min(chr_09_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_09  ') + ylab(' CC4532 Chr_09  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_09_cc4532$sbin_start),max(chr_09_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_09_cc4532$qbin_start),max(chr_09_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_09_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_09.png", plot=cc4532_cc5816_chr_09_dotplot, width=10, height=10) 

#chr_10
cc1690_cc5816_chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690$sbin_start)-4000*2, xmax = min(chr_10_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690$qbin_start)-4000*2, ymax =min(chr_10_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690$sbin_start)-7000*2, xmax = min(chr_10_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690$qbin_start)-7000*2, ymax =min(chr_10_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  #geom_rect(aes(xmin = min(chr_10_cc1690$sbin_start)-7000*2, xmax = min(chr_10_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_10  ') + ylab(' UL1690.1 Chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cc1690$sbin_start),max(chr_10_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cc1690$qbin_start),max(chr_10_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_10_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_10.png", plot=cc1690_cc5816_chr_10_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690N$sbin_start)-4000*2, xmax = min(chr_10_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690N$qbin_start)-4000*2, ymax =min(chr_10_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690N$sbin_start)-7000*2, xmax = min(chr_10_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690N$qbin_start)-7000*2, ymax =min(chr_10_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_10  ') + ylab(' cc1690 Nanopore Chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cc1690N$sbin_start),max(chr_10_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cc1690N$qbin_start),max(chr_10_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_10_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_10.png", plot=Ncc1690_cc5816_chr_10_dotplot, width=10, height=10) 

cw15_cc5816_chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_10"),aes(xmin = min(chr_10_cw15$sbin_start)-4000*2, xmax = min(chr_10_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cw15$qbin_start)-4000*2, ymax =min(chr_10_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_10"),aes(xmin = min(chr_10_cw15$sbin_start)-7000*2, xmax = min(chr_10_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cw15$qbin_start)-7000*2, ymax =min(chr_10_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_10  ') + ylab(' CC400 Chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cw15$sbin_start),max(chr_10_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cw15$qbin_start),max(chr_10_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_10_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_10.png", plot=cw15_cc5816_chr_10_dotplot, width=10, height=10) 

cc4532_cc5816_chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_10"),aes(xmin = min(chr_10_cc4532$sbin_start)-4000*2, xmax = min(chr_10_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc4532$qbin_start)-4000*2, ymax =min(chr_10_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_10"),aes(xmin = min(chr_10_cc4532$sbin_start)-7000*2, xmax = min(chr_10_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc4532$qbin_start)-7000*2, ymax =min(chr_10_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_10  ') + ylab(' CC4532 Chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cc4532$sbin_start),max(chr_10_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cc4532$qbin_start),max(chr_10_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_10_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_10.png", plot=cc4532_cc5816_chr_10_dotplot, width=10, height=10) 

#chr_13
cc1690_cc5816_chr_13_dotplot <- ggplot() +
  geom_point(data=chr_13_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_13"),aes(xmin = min(chr_13_cc1690$sbin_start)-4000*2, xmax = min(chr_13_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc1690$qbin_start)-4000*2, ymax =min(chr_13_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_13"),aes(xmin = min(chr_13_cc1690$sbin_start)-7000*2, xmax = min(chr_13_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc1690$qbin_start)-7000*2, ymax =min(chr_13_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  #geom_rect(aes(xmin = min(chr_13_cc1690$sbin_start)-7000*2, xmax = min(chr_13_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_13  ') + ylab(' UL1690.1 Chr_13  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_13_cc1690$sbin_start),max(chr_13_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_13_cc1690$qbin_start),max(chr_13_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_13_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_13.png", plot=cc1690_cc5816_chr_13_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_13_dotplot <- ggplot() +
  geom_point(data=chr_13_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_13"),aes(xmin = min(chr_13_cc1690N$sbin_start)-4000*2, xmax = min(chr_13_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc1690N$qbin_start)-4000*2, ymax =min(chr_13_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_13"),aes(xmin = min(chr_13_cc1690N$sbin_start)-7000*2, xmax = min(chr_13_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc1690N$qbin_start)-7000*2, ymax =min(chr_13_cc1690N$qbin_start)-9000),fill="red") +
  geom_rect(aes(xmin = min(chr_13_cc1690N$sbin_start)-7000*2, xmax = min(chr_13_cc1690N$sbin_start)-9000, ymin=4384339, ymax = 4384388),fill="black") + 	
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_13  ') + ylab(' cc1690 Nanopore Chr_13  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_13_cc1690N$sbin_start),max(chr_13_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_13_cc1690N$qbin_start),max(chr_13_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_13_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_13.png", plot=Ncc1690_cc5816_chr_13_dotplot, width=10, height=10) 

cw15_cc5816_chr_13_dotplot <- ggplot() +
  geom_point(data=chr_13_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_13"),aes(xmin = min(chr_13_cw15$sbin_start)-4000*2, xmax = min(chr_13_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cw15$qbin_start)-4000*2, ymax =min(chr_13_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_13"),aes(xmin = min(chr_13_cw15$sbin_start)-7000*2, xmax = min(chr_13_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cw15$qbin_start)-7000*2, ymax =min(chr_13_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_13  ') + ylab(' CC400 Chr_13  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_13_cw15$sbin_start),max(chr_13_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_13_cw15$qbin_start),max(chr_13_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_13_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_13.png", plot=cw15_cc5816_chr_13_dotplot, width=10, height=10) 

cc4532_cc5816_chr_13_dotplot <- ggplot() +
  geom_point(data=chr_13_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_13"),aes(xmin = min(chr_13_cc4532$sbin_start)-4000*2, xmax = min(chr_13_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc4532$qbin_start)-4000*2, ymax =min(chr_13_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_13"),aes(xmin = min(chr_13_cc4532$sbin_start)-7000*2, xmax = min(chr_13_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_13"),aes(xmin = Start, xmax = End, ymin=min(chr_13_cc4532$qbin_start)-7000*2, ymax =min(chr_13_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_13  ') + ylab(' CC4532 Chr_13  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_13_cc4532$sbin_start),max(chr_13_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_13_cc4532$qbin_start),max(chr_13_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_13_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_13.png", plot=cc4532_cc5816_chr_13_dotplot, width=10, height=10) 

#chr_14
cc1690_cc5816_chr_14_dotplot <- ggplot() +
  geom_point(data=chr_14_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_14"),aes(xmin = min(chr_14_cc1690$sbin_start)-4000*2, xmax = min(chr_14_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc1690$qbin_start)-4000*2, ymax =min(chr_14_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_14"),aes(xmin = min(chr_14_cc1690$sbin_start)-7000*2, xmax = min(chr_14_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc1690$qbin_start)-7000*2, ymax =min(chr_14_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  #geom_rect(aes(xmin = min(chr_14_cc1690$sbin_start)-7000*2, xmax = min(chr_14_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_14  ') + ylab(' UL1690.1 Chr_14  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_14_cc1690$sbin_start),max(chr_14_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_14_cc1690$qbin_start),max(chr_14_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_14_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_14.png", plot=cc1690_cc5816_chr_14_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_14_dotplot <- ggplot() +
  geom_point(data=chr_14_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_14"),aes(xmin = min(chr_14_cc1690N$sbin_start)-4000*2, xmax = min(chr_14_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc1690N$qbin_start)-4000*2, ymax =min(chr_14_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_14"),aes(xmin = min(chr_14_cc1690N$sbin_start)-7000*2, xmax = min(chr_14_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc1690N$qbin_start)-7000*2, ymax =min(chr_14_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_14  ') + ylab(' cc1690 Nanopore Chr_14  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_14_cc1690N$sbin_start),max(chr_14_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_14_cc1690N$qbin_start),max(chr_14_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_14_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_14.png", plot=Ncc1690_cc5816_chr_14_dotplot, width=10, height=10) 

cw15_cc5816_chr_14_dotplot <- ggplot() +
  geom_point(data=chr_14_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_14"),aes(xmin = min(chr_14_cw15$sbin_start)-4000*2, xmax = min(chr_14_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cw15$qbin_start)-4000*2, ymax =min(chr_14_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_14"),aes(xmin = min(chr_14_cw15$sbin_start)-7000*2, xmax = min(chr_14_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cw15$qbin_start)-7000*2, ymax =min(chr_14_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_14  ') + ylab(' CC400 Chr_14  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_14_cw15$sbin_start),max(chr_14_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_14_cw15$qbin_start),max(chr_14_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_14_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_14.png", plot=cw15_cc5816_chr_14_dotplot, width=10, height=10) 

cc4532_cc5816_chr_14_dotplot <- ggplot() +
  geom_point(data=chr_14_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_14"),aes(xmin = min(chr_14_cc4532$sbin_start)-4000*2, xmax = min(chr_14_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc4532$qbin_start)-4000*2, ymax =min(chr_14_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_14"),aes(xmin = min(chr_14_cc4532$sbin_start)-7000*2, xmax = min(chr_14_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_14"),aes(xmin = Start, xmax = End, ymin=min(chr_14_cc4532$qbin_start)-7000*2, ymax =min(chr_14_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_14  ') + ylab(' CC4532 Chr_14  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_14_cc4532$sbin_start),max(chr_14_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_14_cc4532$qbin_start),max(chr_14_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_14_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_14.png", plot=cc4532_cc5816_chr_14_dotplot, width=10, height=10) 

#chr_15
cc1690_cc5816_chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690$sbin_start)-4000*2, xmax = min(chr_15_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690$qbin_start)-4000*2, ymax =min(chr_15_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690$sbin_start)-7000*2, xmax = min(chr_15_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690$qbin_start)-7000*2, ymax =min(chr_15_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  #geom_rect(aes(xmin = min(chr_15_cc1690$sbin_start)-7000*2, xmax = min(chr_15_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_15  ') + ylab(' UL1690.1 Chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cc1690$sbin_start),max(chr_15_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_15_cc1690$qbin_start),max(chr_15_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_15_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_15.png", plot=cc1690_cc5816_chr_15_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690N$sbin_start)-4000*2, xmax = min(chr_15_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690N$qbin_start)-4000*2, ymax =min(chr_15_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690N$sbin_start)-7000*2, xmax = min(chr_15_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690N$qbin_start)-7000*2, ymax =min(chr_15_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_15  ') + ylab(' cc1690 Nanopore Chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cc1690N$sbin_start),max(chr_15_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_15_cc1690N$qbin_start),max(chr_15_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_15_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_15.png", plot=Ncc1690_cc5816_chr_15_dotplot, width=10, height=10) 

cw15_cc5816_chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_15"),aes(xmin = min(chr_15_cw15$sbin_start)-4000*2, xmax = min(chr_15_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cw15$qbin_start)-4000*2, ymax =min(chr_15_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_15"),aes(xmin = min(chr_15_cw15$sbin_start)-7000*2, xmax = min(chr_15_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cw15$qbin_start)-7000*2, ymax =min(chr_15_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_15  ') + ylab(' CC400 Chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cw15$sbin_start),max(chr_15_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_15_cw15$qbin_start),max(chr_15_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_15_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_15.png", plot=cw15_cc5816_chr_15_dotplot, width=10, height=10) 

cc4532_cc5816_chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_15"),aes(xmin = min(chr_15_cc4532$sbin_start)-4000*2, xmax = min(chr_15_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc4532$qbin_start)-4000*2, ymax =min(chr_15_cc4532$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_15"),aes(xmin = min(chr_15_cc4532$sbin_start)-7000*2, xmax = min(chr_15_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc4532$qbin_start)-7000*2, ymax =min(chr_15_cc4532$qbin_start)-9000),fill="red") +
  geom_rect(aes(xmin = min(chr_15_cc4532$sbin_start)-7000*2, xmax = min(chr_15_cc4532$sbin_start)-9000, ymin=3527513, ymax = 3521336),fill="blue") + 
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_15  ') + ylab(' CC4532 Chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cc4532$sbin_start),max(chr_15_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_15_cc4532$qbin_start),max(chr_15_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_15_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_15.png", plot=cc4532_cc5816_chr_15_dotplot, width=10, height=10) 

#chr_16
cc1690_cc5816_chr_16_dotplot <- ggplot() +
  geom_point(data=chr_16_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_16"),aes(xmin = min(chr_16_cc1690$sbin_start)-4000*2, xmax = min(chr_16_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc1690$qbin_start)-4000*2, ymax =min(chr_16_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_16"),aes(xmin = min(chr_16_cc1690$sbin_start)-7000*2, xmax = min(chr_16_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc1690$qbin_start)-7000*2, ymax =min(chr_16_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  #geom_rect(aes(xmin = min(chr_16_cc1690$sbin_start)-7000*2, xmax = min(chr_16_cc1690$sbin_start)-9000, ymin=4167491, ymax = 4173496),fill="blue") + 
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_16  ') + ylab(' UL1690.1 Chr_16  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_16_cc1690$sbin_start),max(chr_16_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_16_cc1690$qbin_start),max(chr_16_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_16_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_16.png", plot=cc1690_cc5816_chr_16_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_16_dotplot <- ggplot() +
  geom_point(data=chr_16_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_16"),aes(xmin = min(chr_16_cc1690N$sbin_start)-4000*2, xmax = min(chr_16_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc1690N$qbin_start)-4000*2, ymax =min(chr_16_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_16"),aes(xmin = min(chr_16_cc1690N$sbin_start)-7000*2, xmax = min(chr_16_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc1690N$qbin_start)-7000*2, ymax =min(chr_16_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_16  ') + ylab(' cc1690 Nanopore Chr_16  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_16_cc1690N$sbin_start),max(chr_16_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_16_cc1690N$qbin_start),max(chr_16_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_16_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_16.png", plot=Ncc1690_cc5816_chr_16_dotplot, width=10, height=10) 

cw15_cc5816_chr_16_dotplot <- ggplot() +
  geom_point(data=chr_16_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_16"),aes(xmin = min(chr_16_cw15$sbin_start)-4000*2, xmax = min(chr_16_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cw15$qbin_start)-4000*2, ymax =min(chr_16_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_16"),aes(xmin = min(chr_16_cw15$sbin_start)-7000*2, xmax = min(chr_16_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cw15$qbin_start)-7000*2, ymax =min(chr_16_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_16  ') + ylab(' CC400 Chr_16  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_16_cw15$sbin_start),max(chr_16_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_16_cw15$qbin_start),max(chr_16_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_16_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_16.png", plot=cw15_cc5816_chr_16_dotplot, width=10, height=10) 

cc4532_cc5816_chr_16_dotplot <- ggplot() +
  geom_point(data=chr_16_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_16"),aes(xmin = min(chr_16_cc4532$sbin_start)-4000*2, xmax = min(chr_16_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc4532$qbin_start)-4000*2, ymax =min(chr_16_cc4532$qbin_start)-4000,fill=Type)) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_16"),aes(xmin = min(chr_16_cc4532$sbin_start)-7000*2, xmax = min(chr_16_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_16"),aes(xmin = Start, xmax = End, ymin=min(chr_16_cc4532$qbin_start)-7000*2, ymax =min(chr_16_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 Chr_16  ') + ylab(' CC4532 Chr_16  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_16_cc4532$sbin_start),max(chr_16_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_16_cc4532$qbin_start),max(chr_16_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_16_dotplot
#ggsave(file="png_0808/cc4532_cc5816_chr_16.png", plot=cc4532_cc5816_chr_16_dotplot, width=10, height=10) 

#chr_17
cc1690_cc5816_chr_17_dotplot <- ggplot() +
  geom_point(data=chr_17_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_17"),aes(xmin = min(chr_17_cc1690$sbin_start)-4000*2, xmax = min(chr_17_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc1690$qbin_start)-4000*2, ymax =min(chr_17_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_17"),aes(xmin = min(chr_17_cc1690$sbin_start)-7000*2, xmax = min(chr_16_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc1690$qbin_start)-7000*2, ymax =min(chr_17_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_rect(aes(xmin = min(chr_17_cc1690$sbin_start)-7000*2, xmax = min(chr_17_cc1690$sbin_start)-9000, ymin=5931558, ymax = 5923908),fill="blue") + 	
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_17  ') + ylab(' UL1690.1 chr_17  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_17_cc1690$sbin_start),max(chr_17_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_17_cc1690$qbin_start),max(chr_17_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_17_dotplot
ggsave(file="png_0808/cc1690_cc5816_chr_17.png", plot=cc1690_cc5816_chr_17_dotplot, width=10, height=10) 

Ncc1690_cc5816_chr_17_dotplot <- ggplot() +
  geom_point(data=chr_17_cc1690N, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690_nanopore,ZeppL_1690_nanopore$Chr=="chr_17"),aes(xmin = min(chr_17_cc1690N$sbin_start)-4000*2, xmax = min(chr_17_cc1690N$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc1690N$qbin_start)-4000*2, ymax =min(chr_17_cc1690N$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690_nanopore",label$Chr=="chr_17"),aes(xmin = min(chr_17_cc1690N$sbin_start)-7000*2, xmax = min(chr_17_cc1690N$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc1690N$qbin_start)-7000*2, ymax =min(chr_17_cc1690N$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_17  ') + ylab(' cc1690 Nanopore chr_17  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_17_cc1690N$sbin_start),max(chr_17_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_17_cc1690N$qbin_start),max(chr_17_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_17_dotplot
ggsave(file="png_0808/Nano_1690_cc5816_chr_17.png", plot=Ncc1690_cc5816_chr_17_dotplot, width=10, height=10) 

cw15_cc5816_chr_17_dotplot <- ggplot() +
  geom_point(data=chr_17_cw15, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_cw15,ZeppL_cw15$Chr=="chr_17"),aes(xmin = min(chr_17_cw15$sbin_start)-4000*2, xmax = min(chr_17_cw15$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cw15$qbin_start)-4000*2, ymax =min(chr_17_cw15$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cw15",label$Chr=="chr_17"),aes(xmin = min(chr_17_cw15$sbin_start)-7000*2, xmax = min(chr_17_cw15$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cw15$qbin_start)-7000*2, ymax =min(chr_17_cw15$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_17  ') + ylab(' CC400 chr_17  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_17_cw15$sbin_start),max(chr_17_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_17_cw15$qbin_start),max(chr_17_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_17_dotplot
ggsave(file="png_0808/cc400_cc5816_chr_17.png", plot=cw15_cc5816_chr_17_dotplot, width=10, height=10) 

cc4532_cc5816_chr_17_dotplot <- ggplot() +
  geom_point(data=chr_17_cc4532, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_4532,ZeppL_4532$Chr=="chr_17"),aes(xmin = min(chr_17_cc4532$sbin_start)-4000*2, xmax = min(chr_17_cc4532$sbin_start)-4000, ymin=Start, ymax = End,fill=Type)) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc4532$qbin_start)-4000*2, ymax =min(chr_17_cc4532$qbin_start)-4000,fill=Type)) +
  geom_rect(data=filter(label,label$Strain=="cc4532",label$Chr=="chr_17"),aes(xmin = min(chr_17_cc4532$sbin_start)-7000*2, xmax = min(chr_17_cc4532$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_17"),aes(xmin = Start, xmax = End, ymin=min(chr_17_cc4532$qbin_start)-7000*2, ymax =min(chr_17_cc4532$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_17  ') + ylab(' CC4532 chr_17  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_17_cc4532$sbin_start),max(chr_17_cc4532$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_17_cc4532$qbin_start),max(chr_17_cc4532$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc4532_cc5816_chr_17_dotplot
ggsave(file="png_0808/cc4532_cc5816_chr_17.png", plot=cc4532_cc5816_chr_17_dotplot, width=10, height=10) 
