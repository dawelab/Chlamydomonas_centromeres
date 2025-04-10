setwd("/Users/mingyuwang/Desktop/Chlamydomonas_summary_2024_12_02/telo_chlamy_plot")
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
library(ggchicklet)
library(cowplot)
library(patchwork)
library(GenomicRanges)
options(scipen = 999)

read_cov <- setNames(read.table("cc1690_cc5816.bedgraph",sep = "\t",header=FALSE),c("chr","start","end","coverage"))
read_telo <- setNames(read.table("cc1690_telo_reads_filtered_trimmed.txt.sub.bed",sep = "\t",header=FALSE),c("chr","start","end"))
gene_bed <- setNames(read.table("UL1690_gene_allmodel.bed", sep = "\t",header=FALSE) ,c("chr","start","end"))
UL1690_chr_length <- setNames(read.table("cc1690_HiFi_chr_length.txt", sep = "\t",header=FALSE) ,c("chr","length"))


calculate_gene_density <- function(gene_bed, chr_length, tilewidth = 10000) {
  seqlengths_vector <- setNames(chr_length$length, chr_length$chr)
  gene_ranges <- GRanges(seqnames = gene_bed$chr,
                         ranges = IRanges(start = gene_bed$start, end = gene_bed$end))
  sliding_windows <- tileGenome(seqlengths = seqlengths_vector,
                                tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
  gene_density <- countOverlaps(sliding_windows, gene_ranges)
  result <- data.frame(chrom = seqnames(sliding_windows),
                       start = start(sliding_windows),
                       end = end(sliding_windows),
                       gene_density = gene_density)
  return(result)
}
gene_density <- calculate_gene_density(gene_bed,UL1690_chr_length,  tilewidth = 10000)


chr2_df <- read_cov %>% filter(chr=="chr_02")
chr2_telo <- read_telo %>% filter(chr=="chr_02")
chr2_gene <-gene_density %>% filter(chrom=="chr_02")

chr2_plot<- ggplot() +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=0,xmax=8585683,ymin=0,ymax=1,fill="#525252"), r = unit(0.1, "npc")) +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=2555000,xmax=2620000,ymin=0,ymax=1),fill="#66CC66",color="#66CC66", r = unit(0.1, "npc")) +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=8430000,xmax=8460000,ymin=0,ymax=1),fill="#E69F00",color="#E69F00", r = unit(0.1, "npc")) +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=8441143,xmax=8449477,ymin=0,ymax=1),fill="#9966FF", color="#9966FF",r = unit(0.1, "npc")) +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=0,xmax=2000,ymin=0,ymax=1),fill="red", color="red",r = unit(0.1, "npc")) +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=8584683,xmax=8585683,ymin=0,ymax=1),fill="red",color="red", r = unit(0.1, "npc")) +
  scale_fill_manual(
    values = "#525252"
  ) +
  theme_bw() +
  ylim(0,1.3) +
  theme(
    axis.title.x = element_text(size = 14),   # Keep x-axis title
    axis.title.y = element_blank(),           # Remove y-axis title
    axis.text.x = element_text(size = 12),    # Keep x-axis text
    axis.text.y = element_blank(),            # Remove y-axis text
    axis.ticks.x = element_line(),            # Keep x-axis ticks
    axis.ticks.y = element_blank(),           # Remove y-axis ticks
    plot.title = element_blank(),             # Remove plot title
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-10, 10, 10, 0)
  )

chr2_cov_plot<- ggplot() +
  geom_area(data = chr2_df, aes(x = start, y = coverage), alpha = 0.4, fill = "#3399CC", color = "#3399CC") +
  geom_hline(yintercept = mean(chr2_df$coverage), color = "red", linetype = "dashed", size = 0.8) +
  theme_bw() +
  ylab("") +
  scale_y_continuous(
    breaks = c(0,300, 600),
    limits = c(0, 600)) +
  theme(
    axis.title.y  = element_text(size = 10),  
    axis.title.x  = element_blank(), 
    axis.text.x = element_blank(),
    #legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(20, 10, 10, 0),
    axis.line.x = element_line(color = "black"),  # Add bottom border
    axis.line.y = element_line(color = "black")   # Add left border
  )

chr2_gene_plot<- ggplot() +
  #geom_bar(data = chr2_gene, aes(x = start, y = gene_density), alpha = 0.4, fill = "#339966", color = "#339966",stat="identity") +
  geom_area(data = chr2_gene,
    aes(x = start, y = gene_density),
    alpha = 0.4,
    fill = "#339966",
    color = "#339966"
  ) +
  theme_bw() +
  scale_y_continuous(
    breaks = c(0, 10),
    limits = c(0, 10)) +
  ylab("") +
  theme(
    axis.title.y = element_text(size = 10),  
    axis.title.x  = element_blank(), 
    axis.text.x = element_blank(),
    #legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-10, 10, 10, 0),
    axis.line.x = element_line(color = "black"),  # Add bottom border
    axis.line.y = element_line(color = "black")   # Add left border
  )
check_length <- read.table("check_telo_reads_length.txt",sep = "\t",header=FALSE)
UL1690_telo_blast <- read.table("cc1690_telo_reads_filtered_trimmed.txt.sub",sep = "\t",header=FALSE)
chr02_telo_blast <- UL1690_telo_blast %>% filter(V2=="chr_02") %>%
  mutate(
    V1.1 = str_split(V1, "[:\\-]", simplify = TRUE)[, 1],  # First part
    V1.2 = as.numeric(str_split(V1, "[:\\-]", simplify = TRUE)[, 2]),  # Second part
    V1.3 = as.numeric(str_split(V1, "[:\\-]", simplify = TRUE)[, 3])   # Third part
  ) %>%
  filter(V9>=2600000, V10<=3600000) %>% arrange(V1.1) %>% left_join(check_length,by=c("V1.1"="V1"))



chr02_telo_blast_reg1 <- chr02_telo_blast %>% filter(V10<3000000) %>% mutate(order=row_number()) %>% mutate(telo_length=abs(V2.y-V1.3-V1.2))
chr02_telo_blast_reg2 <- chr02_telo_blast %>% filter(V10>3000000) %>% mutate(order=row_number()) %>% mutate(telo_length=abs(V2.y-(V1.3-V1.2)))

chr2_regs <- ggplot()+
  geom_rect(data=chr02_telo_blast_reg1,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  geom_rect(data=chr02_telo_blast_reg2,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  theme_bw() +
  ylim(0,40) +
  xlim(0,8585683) +
  ylab("") +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-90, 10, 10, 10)
  )

zoom_chr2_reg1 <- ggplot()+
  geom_rect(data=chr02_telo_blast_reg1,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  #geom_rect(data=chr02_telo_blast_reg2,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+1),fill="black") +
  geom_rect(data=chr02_telo_blast_reg1,aes(xmin=2782504-telo_length,xmax=2782504,ymin=order,ymax=order+0.9),fill="red",alpha=0.7) +
  #geom_rect(data=chr02_telo_blast_reg2,aes(xmin=3398748,xmax=3398748+V2.y-(V1.3-V1.2),ymin=order,ymax=order+1),fill="red") +
  theme_bw() +
  ylim(0,40) +
  xlim(2775000,2875000) +
  ylab("") +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-10, 10, 10, 0)
  )

zoom_chr2_reg2 <- ggplot()+
  #geom_rect(data=chr02_telo_blast_reg1,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  geom_rect(data=chr02_telo_blast_reg2,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  #geom_rect(data=chr02_telo_blast_reg1,aes(xmin=2782504-telo_length,xmax=2782504,ymin=order,ymax=order+0.9),fill="red",alpha=0.7) +
  geom_rect(data=chr02_telo_blast_reg2,aes(xmin=3398748,xmax=3398748+V2.y-(V1.3-V1.2),ymin=order,ymax=order+0.9),fill="red",alpha=0.7) +
  theme_bw() +
  ylim(0,40) +
  xlim(3370000,3470000) +
  ylab("") +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-10, 10, 10, 0)
  )
zoom_chr2_plot<- ggplot() +
  ggchicklet:::geom_rrect(data=NULL,aes(xmin=2775000,xmax=2810000,ymin=0,ymax=1,fill="grey"), r = unit(0.1, "npc")) +
  scale_fill_manual(
    values = "grey"
  ) +
  theme_bw() +
  ylim(0,1.3) +
  theme(
    axis.title.x = element_text(size = 14),   # Keep x-axis title
    axis.title.y = element_blank(),           # Remove y-axis title
    axis.text.x = element_text(size = 12),    # Keep x-axis text
    axis.text.y = element_blank(),            # Remove y-axis text
    axis.ticks.x = element_blank(),            # Keep x-axis ticks
    axis.ticks.y = element_blank(),           # Remove y-axis ticks
    plot.title = element_blank(),             # Remove plot title
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-20, 10, 10, 10)
  )

zoom_chr2_reg1
zoom_chr2_reg2
zoom_chr2_plot


com_plot <- plot_grid(chr2_cov_plot,chr2_gene_plot,chr2_plot,chr2_regs,nrow = 4,rel_heights = c(1,0.2,0.4,1),align = "v")
zoom_com_plot <- plot_grid(zoom_chr2_reg1,zoom_chr2_reg2,nrow = 2,rel_heights = c(1,1),align = "v")
ggsave(file="chr2.png", com_plot, width=10, height=4, dpi=300) 
ggsave(file="chr2_zoom.png", zoom_com_plot, width=10, height=8, dpi=300) 

chr_02_gene_bed <- gene_bed %>% filter(chr=="chr_02")
chr2_regs <- ggplot()+
  geom_rect(data=chr02_telo_blast_reg1,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  geom_rect(data=chr02_telo_blast_reg2,aes(xmin=V9,xmax=V10,ymin=order,ymax=order+0.9),fill="black") +
  geom_rect(data=chr02_telo_blast_reg1,aes(xmin=2782504-telo_length,xmax=2782504,ymin=order,ymax=order+0.9),alpha=0.7) +
  theme_bw() +
  ylim(0,40) +
  xlim(0,8585683) +
  ylab("") +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    legend.position = "none",                 # Remove legend
    panel.grid = element_blank(),             # Remove grid lines
    panel.background = element_blank(),       # Remove background color
    panel.border = element_blank() ,           # Remove border
    plot.margin = margin(-90, 10, 10, 10)
  )