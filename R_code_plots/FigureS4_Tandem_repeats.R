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
UL1690_TRF <- setNames(read.table('cc1690_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
UL1690_TRF_filtered <- UL1690_TRF %>% filter(Copies>400 & Period>100)
UL1690_TRF_filtered <- UL1690_TRF_filtered %>%
  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

cc5816_TRF <- setNames(read.table('cc5816_Dutcher.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
cc5816_TRF_filtered <- cc5816_TRF %>% filter(Copies>400 & Period>125)
cc5816_TRF_filtered <- cc5816_TRF_filtered %>%
  mutate(order = as.numeric(sub("chr_", "", Chr)))  %>% mutate(order=17-order+1)

cc400_TRF <- setNames(read.table('cw15_HiFi_chr.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
cc400_TRF_filtered <- cc400_TRF %>% filter(Copies>400 & Period>125)
cc400_TRF_filtered <- cc400_TRF_filtered %>%
  mutate(order = as.numeric(sub("chr_", "", Chr)))  %>% mutate(order=17-order+1)

cc4532_TRF <- setNames(read.table('cc4532_v6.1_chr.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
cc4532_TRF_filtered <- cc4532_TRF %>% filter(Copies>400 & Period>125)
cc4532_TRF_filtered <- cc4532_TRF_filtered %>%
  mutate(order = as.numeric(sub("chr_", "", Chr)))  %>% mutate(order=17-order+1)

cc1690_TRF <- setNames(read.table('cc1690_nanopore.fa.2.7.7.80.10.50.2000.dat.bed', sep="\t", fill=TRUE),name_list)
cc1690_TRF_filtered <- cc1690_TRF %>% filter(Copies>400 & Period>125)
cc1690_TRF_filtered <- cc1690_TRF_filtered %>%
  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

UL1690_length <- setNames(read.table('UL1690_chr_length.txt', sep="\t", fill=TRUE),c("Chr","Length"))
UL1690_length <- UL1690_length%>% mutate(order = n() - row_number() + 1) 
UL1690_gap <- setNames(read.table('UL1690_gap.bed', sep="\t", fill=TRUE),c("Chr","Start","End"))
UL1690_gap <- UL1690_gap %>%  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

cc5816_length <- setNames(read.table('cc5816_Dutcher_length.txt', sep="\t", fill=TRUE),c("Chr","Length"))
cc5816_length <- cc5816_length%>% mutate(order = n() - row_number() + 1) 

cc4532_length <- setNames(read.table('cc4532_v6.1_chr_length.txt', sep="\t", fill=TRUE),c("Chr","Length"))
cc4532_length <- cc4532_length%>% mutate(order = n() - row_number() + 1) 
cc4532_gap <- setNames(read.table('cc4532_gap.bed', sep="\t", fill=TRUE),c("Chr","Start","End"))
cc4532_gap <- cc4532_gap %>%  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

cc1690_length <- setNames(read.table('cc1690_nanopore_length.txt', sep="\t", fill=TRUE),c("Chr","Length"))
cc1690_length <- cc1690_length%>% mutate(order = n() - row_number() + 1) 
cc1690_gap <- setNames(read.table('cc1690_gap.bed', sep="\t", fill=TRUE),c("Chr","Start","End"))
cc1690_gap <- cc1690_gap %>%  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

cc400_length <- setNames(read.table('cw15_HiFi_chr_length.txt', sep="\t", fill=TRUE),c("Chr","Length"))
cc400_length <- cc400_length%>% mutate(order = n() - row_number() + 1) 
cc400_gap <- setNames(read.table('cw15_gap.bed', sep="\t", fill=TRUE),c("Chr","Start","End"))
cc400_gap <- cc400_gap %>%  mutate(order = as.numeric(sub("chr_", "", Chr))) %>% mutate(order=17-order+1)

UL1690 <- ggplot() +
  geom_rect(data=UL1690_length,aes(xmin = 0, xmax =Length, ymin=order, ymax = order+0.8),fill="#0072B2",color="black") + 
  geom_rect(data=UL1690_TRF_filtered,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="red",color="black") +
  geom_rect(data=UL1690_gap,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="lightblue",color="white") + 
  geom_text(data = UL1690_length, aes(x =-600000, y = order+0.5, label = Chr),color = "black", size = 5) +
  labs(title = "UL1690 Tandem Repeats", x = "", y = "") +
  theme_bw()+
  ylim(0, 18) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40),
        panel.grid = element_blank(),          # Remove grid lines
        axis.text.y = element_blank(),          # Hide y-axis text
        axis.ticks.y = element_blank(),         # Hide y-axis ticks
        axis.title.y = element_blank()          # Hide y-axis title
        )
UL1690
ggsave(file="FigureS4/UL1690_repeats.png", UL1690, width=25, height=8)

cc5816 <- ggplot() +
  geom_rect(data=cc5816_length,aes(xmin = 0, xmax =Length, ymin=order, ymax = order+0.8),fill="#0072B2",color="black") + 
  geom_rect(data=cc5816_TRF_filtered,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="red",color="black") + 
  geom_text(data = cc5816_length, aes(x =-600000, y = order+0.5, label = Chr),color = "black", size = 5) +
  labs(title = "CC5816 Tandem Repeats", x = "", y = "") +
  theme_bw()+
  ylim(0, 18) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40),
        panel.grid = element_blank(),          # Remove grid lines
        axis.text.y = element_blank(),          # Hide y-axis text
        axis.ticks.y = element_blank(),         # Hide y-axis ticks
        axis.title.y = element_blank()          # Hide y-axis title
  )
cc5816
ggsave(file="FigureS4/cc5816_repeats.png", cc5816, width=25, height=8)

cc400 <- ggplot() +
  geom_rect(data=cc400_length,aes(xmin = 0, xmax =Length, ymin=order, ymax = order+0.8),fill="#0072B2",color="black") + 
  geom_rect(data=cc400_TRF_filtered,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="red",color="black") + 
  geom_rect(data=cc400_gap,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="lightblue",color="white") + 
  geom_text(data = cc400_length, aes(x =-600000, y = order+0.5, label = Chr),color = "black", size = 5) +
  labs(title = "CC400 Tandem Repeats", x = "", y = "") +
  theme_bw()+
  ylim(0, 18) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40),
        panel.grid = element_blank(),          # Remove grid lines
        axis.text.y = element_blank(),          # Hide y-axis text
        axis.ticks.y = element_blank(),         # Hide y-axis ticks
        axis.title.y = element_blank()          # Hide y-axis title
  )
cc400
ggsave(file="FigureS4/cc400_repeats.png", cc400, width=25, height=8)

cc4532 <- ggplot() +
  geom_rect(data=cc4532_length,aes(xmin = 0, xmax =Length, ymin=order, ymax = order+0.8),fill="#0072B2",color="black") + 
  geom_rect(data=cc4532_TRF_filtered,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="red",color="black") + 
  geom_rect(data=cc4532_gap,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="lightblue",color="white") + 
  geom_text(data = cc4532_length, aes(x =-600000, y = order+0.5, label = Chr),color = "black", size = 5) +
  labs(title = "CC4532 Tandem Repeats", x = "", y = "") +
  theme_bw()+
  ylim(0, 18) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40),
        panel.grid = element_blank(),          # Remove grid lines
        axis.text.y = element_blank(),          # Hide y-axis text
        axis.ticks.y = element_blank(),         # Hide y-axis ticks
        axis.title.y = element_blank()          # Hide y-axis title
  )
cc4532
ggsave(file="FigureS4/cc4532_repeats.png", cc4532, width=25, height=8)

cc1690 <- ggplot() +
  geom_rect(data=cc1690_length,aes(xmin = 0, xmax =Length, ymin=order, ymax = order+0.8),fill="#0072B2",color="black") + 
  geom_rect(data=cc1690_TRF_filtered,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="red",color="black") + 
  geom_rect(data=cc1690_gap,aes(xmin = Start, xmax =End, ymin=order, ymax = order+0.8),fill="lightblue",color="white") + 
  geom_text(data = cc1690_length, aes(x =-600000, y = order+0.5, label = Chr),color = "black", size = 5) +
  labs(title = "CC1690 Tandem Repeats", x = "", y = "") +
  theme_bw()+
  ylim(0, 18) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40),
        panel.grid = element_blank(),          # Remove grid lines
        axis.text.y = element_blank(),          # Hide y-axis text
        axis.ticks.y = element_blank(),         # Hide y-axis ticks
        axis.title.y = element_blank()          # Hide y-axis title
  )
cc1690
ggsave(file="FigureS4/cc1690_repeats.png", cc1690, width=25, height=8)
# Generate information for TableS4
UL1690_print<- UL1690_TRF_filtered %>% select(Chr, Start,End,Period,Copies,Repeats)
cc1690_print<- cc1690_TRF_filtered %>% select(Chr, Start,End,Period,Copies,Repeats)
cc5816_print<- cc5816_TRF_filtered %>% select(Chr, Start,End,Period,Copies,Repeats)
cc400_print<- cc400_TRF_filtered %>% select(Chr, Start,End,Period,Copies,Repeats)
cc4532_print<- cc4532_TRF_filtered %>% select(Chr, Start,End,Period,Copies,Repeats)

write.table(UL1690_print, "UL1690_TRF_filtered.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cc5816_print, "cc5816_TRF_filtered.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cc400_print, "cc400_TRF_filtered.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cc4532_print, "cc4532_TRF_filtered.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cc1690_print, "cc1690_TRF_filtered.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
