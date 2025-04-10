setwd("/Users/mingyuwang/Desktop/telo_chlamy")
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
options(scipen = 999)
##### Part I HiFi reads
# Step 1: Data Preparation
telo_hifi_reads <- setNames(read.table('cc1690_telo_1000bp.bed.sub', sep="\t", fill=TRUE), c("name", "start", "end", "size"))
telo_hifi_reads_length <- setNames(read.table('cc1690_telo_reads_length.txt', sep="\t", fill=TRUE), c("name", "length"))
merge_hifi_telo <- telo_hifi_reads %>% 
  inner_join(telo_hifi_reads_length, by = "name") %>%
  filter(size >27) %>%
  group_by(name) %>%
  slice_max(size,n=1) %>%
  ungroup() %>%
  arrange(length) %>%
  mutate( order = row_number(), 
          ratio = round(size / length * 100, 2))

#Plot
telo_plot <- ggplot() +
  geom_rect(data=merge_hifi_telo,aes(xmin = order, xmax =order+1, ymin=0, ymax = length),fill="#0072B2",color="black") + 
  geom_rect(data=merge_hifi_telo,aes(xmin = order, xmax =order+1, ymin=start, ymax = end),fill="red",color="black") + 
  #geom_text(data = merge_hifi_telo, aes(x =order+0.5, y = length+1500, label = ratio),color = "black", size = 4) +
  labs(title = "Telomere size (>27bp) Distribution", x = "read ID", y = "size") +
  theme_bw()+
  ylim(0, 33000) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40))
telo_plot
ggsave(file="png/telo_plot.png", telo_plot, width=25, height=8) 

merge_hifi_telo200 <- merge_hifi_telo %>% filter(size>200)
telo_plot200 <- ggplot() +
  geom_rect(data=merge_hifi_telo200,aes(xmin = order, xmax =order+1, ymin=0, ymax = length),fill="#0072B2",color="black") + 
  geom_rect(data=merge_hifi_telo200,aes(xmin = order, xmax =order+1, ymin=start, ymax = end),fill="red",color="black") + 
  #geom_text(data = merge_hifi_telo, aes(x =order+0.5, y = length+1500, label = ratio),color = "black", size = 4) +
  labs(title = "Telomere size (>200bp) Distribution", x = "read ID", y = "size") +
  theme_bw()+
  ylim(0, 33000) +
  theme(axis.title = element_text(size = 30) ,axis.text = element_text(size = 30),plot.title=element_text(size = 40))
telo_plot200
ggsave(file="png/telo_plot200.png", telo_plot200, width=25, height=8) 

# Prepare Final Data Frame
HiFi_left_telo <- merge_hifi_telo %>% filter(start < 200) %>% select(name, start = end, end = length)
HiFi_right_telo <- merge_hifi_telo %>% filter(start > 200) %>%
  mutate(nstart = 0) %>% 
  select(name, start=nstart, end=start)

HiFi_left_telo200 <- merge_hifi_telo200 %>% filter(start < 200) %>% select(name, start = end, end = length)
HiFi_right_telo200 <- merge_hifi_telo200 %>% filter(start > 200) %>%
  mutate(nstart = 0) %>% 
  select(name, start=nstart, end=start)

# Combine left and right data, and save, using as bed file to trim telomere sequences
HiFi_reads_final <- bind_rows(HiFi_left_telo, HiFi_right_telo )
write.table(HiFi_reads_final, "png/telo_reads_final_HiFi.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

HiFi_reads_final200 <- bind_rows(HiFi_left_telo200, HiFi_right_telo200 )
write.table(HiFi_reads_final200, "png/telo_reads_final_HiFi200.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


telo_reads_blast <- setNames(read.table('cc1690_telo_reads_filtered_trimmed.txt.sub.bed2', sep="\t", fill=TRUE), c("chr", "chr_start", "chr_end", "name","trim_start","trim_end"))
chr2_reads <- telo_reads_blast %>% filter(chr=="chr_02")
chr2_reads_info <- chr2_reads %>%
  inner_join(merge_hifi_telo ,by="name") 
colnames(chr2_reads_info)[c(7, 8)] <- c("telo_start", "telo_end")

chr2_reads_info_filter_region1 <- chr2_reads_info %>% filter(chr_start>2700000, chr_end < 3300000)
chr2_reads_info_filter_region2 <- chr2_reads_info %>% filter(chr_start>3000000, chr_end < 3400000)

check_read <- data.frame(rbind(chr2_reads_info_filter_region1,chr2_reads_info_filter_region2))
write.table(check_read , "png/check_reads.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

