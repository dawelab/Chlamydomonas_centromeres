setwd("/Users/mingyuwang/Desktop/Chlamydomonas_summary_2024_12_02/2024_08Chlamy/dotplot_update/ZeppL_compare")
library(dplyr)
library(magrittr)
library(knitr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(cowplot)

ZeppL_1690 <- read.table('cc1690_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_1690) <- c("Chr","Start","End","Type")
new_rows <- data.frame(
  Chr = c("chr_10","chr_15"),
  Start = c(3565164,3186546),
  End = c(3575164, 3196546),
  Type = c("Non ZeppL","Non ZeppL")
)
ZeppL_1690 <- rbind(ZeppL_1690, new_rows)
#chr_10	3565164	3575164 Non ZeppL
#chr_15	3186546 3196546 Non ZeppL

ZeppL_1690_nanopore <- read.table('cc1690_nanopore_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_1690_nanopore) <- c("Chr","Start","End","Type")
new_rows <- data.frame(
  Chr = c("chr_10","chr_15"),
  Start = c(3567834,3324962),
  End = c(3577834, 3334962),
  Type = c("Non ZeppL","Non ZeppL")
)
ZeppL_1690_nanopore <- rbind(ZeppL_1690_nanopore, new_rows)
#chr_10	3567834	3577834 Non ZeppL
#chr_15	3324962 3334962 Non ZeppL
ZeppL_cw15 <- read.table('cw15_ZeppL_region_bar_v3.bed', sep="\t", fill=TRUE)
colnames(ZeppL_cw15) <- c("Chr","Start","End","Type")
new_rows <- data.frame(
  Chr = c("chr_10","chr_15"),
  Start = c(3559611,3061428),
  End = c(3569609,3148367),
  Type = c("Non ZeppL","Non ZeppL")
)
ZeppL_cw15 <- rbind(ZeppL_cw15, new_rows)
#chr_10	3559611	3569609 Non ZeppL
#chr_15 3061428 3148367 Non ZeppL
ZeppL_5816 <- read.table('cc5816_ZeppL_region_bar.bed', sep="\t", fill=TRUE)
colnames(ZeppL_5816) <- c("Chr","Start","End","Type")
new_rows <- data.frame(
  Chr = c("chr_10","chr_15"),
  Start = c(3567021,3232113),
  End = c(3577021,3242113),
  Type = c("Non ZeppL","Non ZeppL")
)
ZeppL_5816 <- rbind(ZeppL_5816, new_rows)
#chr_10	3567021	3577021 Non ZeppL
#chr_15	3232113 3242113 Non ZeppL

chr_10_cc1690 <- read.table('cc1690_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_10_cc1690N <- read.table('cc1690_nanopore_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)
chr_10_cw15 <- read.table('cw15_cc5816_chr_10.blast.txt.sub', sep="\t", fill=TRUE)

chr_15_cc1690 <- read.table('cc1690_chr_15_blast.txt.sub', sep="\t", fill=TRUE)
chr_15_cc1690N <- read.table('cc1690_nanopore_chr_15_blast.txt.sub', sep="\t", fill=TRUE)
chr_15_cw15 <- read.table('cw15_chr_15_blast.txt.sub', sep="\t", fill=TRUE)

data_list <- c("chr_10_cc1690","chr_10_cc1690N","chr_10_cw15",
               "chr_15_cc1690","chr_15_cc1690N","chr_15_cw15")
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
#chr_10
cc1690_cc5816_chr_10_dotplot <- ggplot() +
  geom_point(data=chr_10_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690$sbin_start)-4000*2, xmax = min(chr_10_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690$qbin_start)-4000*2, ymax =min(chr_10_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_10"),aes(xmin = min(chr_10_cc1690$sbin_start)-7000*2, xmax = min(chr_10_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_10"),aes(xmin = Start, xmax = End, ymin=min(chr_10_cc1690$qbin_start)-7000*2, ymax =min(chr_10_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_10  ') + ylab(' UL1690.1 chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cc1690$sbin_start),max(chr_10_cc1690$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cc1690$qbin_start),max(chr_10_cc1690$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_10_dotplot
ggsave(file="png_0408/cc1690_cc5816_chr_10.png", plot=cc1690_cc5816_chr_10_dotplot, width=10, height=10) 

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
  #xlab('CC5816 chr_10  ') + ylab(' cc1690 Nanopore chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cc1690N$sbin_start),max(chr_10_cc1690N$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cc1690N$qbin_start),max(chr_10_cc1690N$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_10_dotplot
ggsave(file="png_0408/Nano_1690_cc5816_chr_10.png", plot=Ncc1690_cc5816_chr_10_dotplot, width=10, height=10) 

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
  #xlab('CC5816 chr_10  ') + ylab(' CC400 chr_10  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_10_cw15$sbin_start),max(chr_10_cw15$sbin_end),50000))+
  scale_y_continuous(breaks =seq(min(chr_10_cw15$qbin_start),max(chr_10_cw15$qbin_end),50000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_10_dotplot
ggsave(file="png_0408/cc400_cc5816_chr_10.png", plot=cw15_cc5816_chr_10_dotplot, width=10, height=10) 


#chr_15
cc1690_cc5816_chr_15_dotplot <- ggplot() +
  geom_point(data=chr_15_cc1690, aes(x=sbin_start+sstart, y=qbin_start, colour=bitscore),size=0.5) +
  geom_rect(data=filter(ZeppL_1690,ZeppL_1690$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690$sbin_start)-4000*2, xmax = min(chr_15_cc1690$sbin_start)-4000, ymin=Start, ymax = End,fill=Type) ) +
  geom_rect(data=filter(ZeppL_5816,ZeppL_5816$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690$qbin_start)-4000*2, ymax =min(chr_15_cc1690$qbin_start)-4000,fill=Type) ) +
  geom_rect(data=filter(label,label$Strain=="cc1690",label$Chr=="chr_15"),aes(xmin = min(chr_15_cc1690$sbin_start)-7000*2, xmax = min(chr_15_cc1690$sbin_start)-9000, ymin=Start, ymax = End),fill="red") +
  geom_rect(data=filter(label,label$Strain=="cc5816",label$Chr=="chr_15"),aes(xmin = Start, xmax = End, ymin=min(chr_15_cc1690$qbin_start)-7000*2, ymax =min(chr_15_cc1690$qbin_start)-9000),fill="red") +
  scale_fill_manual(values = c("white","grey","black"),labels=c("","Non ZeppL","ZeppL"))+
  geom_segment() +
  scale_colour_gradient(low="#FDE725FF", high="#440154FF") +
  theme_bw() +
  #xlab('CC5816 chr_15  ') + ylab(' UL1690.1 chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cc1690$sbin_start),max(chr_15_cc1690$sbin_end),100000))+
  scale_y_continuous(breaks =seq(min(chr_15_cc1690$qbin_start),max(chr_15_cc1690$qbin_end),100000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cc1690_cc5816_chr_15_dotplot
ggsave(file="png_0408/cc1690_cc5816_chr_15.png", plot=cc1690_cc5816_chr_15_dotplot, width=10, height=10) 

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
  #xlab('CC5816 chr_15  ') + ylab(' cc1690 Nanopore chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cc1690N$sbin_start),max(chr_15_cc1690N$sbin_end),100000))+
  scale_y_continuous(breaks =seq(min(chr_15_cc1690N$qbin_start),max(chr_15_cc1690N$qbin_end),100000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
Ncc1690_cc5816_chr_15_dotplot
ggsave(file="png_0408/Nano_1690_cc5816_chr_15.png", plot=Ncc1690_cc5816_chr_15_dotplot, width=10, height=10) 

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
  #xlab('CC5816 chr_15  ') + ylab(' CC400 chr_15  ') +
  xlab('') + ylab('') +
  scale_x_continuous(breaks =seq(min(chr_15_cw15$sbin_start),max(chr_15_cw15$sbin_end),100000))+
  scale_y_continuous(breaks =seq(min(chr_15_cw15$qbin_start),max(chr_15_cw15$qbin_end),100000))+
  coord_fixed() +
  theme(legend.position = "none",axis.text = element_text(size = 20))
cw15_cc5816_chr_15_dotplot
ggsave(file="png_0408/cc400_cc5816_chr_15.png", plot=cw15_cc5816_chr_15_dotplot, width=10, height=10) 
