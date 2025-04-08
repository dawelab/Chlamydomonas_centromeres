setwd('/Users/mingyuwang/Desktop/2024_07Chlamy/Figure0422')
library(xlsx)
library(tidyverse)
library(ggplot2)
library(mdthemes)
library(ggrepel)


cc1690_hifi_info <- read.xlsx('cc1690_hifi_info.xlsx', sheetName = "Sheet1")
cc1690_hifi_info <- cc1690_hifi_info[,-c(1:1)]
(ZeppL_length <- cc1690_hifi_info[1,] %>% unlist(use.names = F))
(Max_length <- cc1690_hifi_info[2,] %>% unlist(use.names = F))
(Min_length <- cc1690_hifi_info[3,] %>% unlist(use.names = F))
(chr_name <- c("chr_01", "chr_02","chr_02_2","chr_03", "chr_04","chr_05", "chr_06","chr_07", "chr_08","chr_09", "chr_10","chr_11", "chr_12","chr_13", "chr_14","chr_15", "chr_16","chr_17"))

df<- data.frame(chr = rep(chr_name,times=3),
                length = c(ZeppL_length,
                           Max_length,
                           Min_length),
                group = rep(c("ZeppL-LINE1 region","Upper centromere size estimate","Lower centromere size estimate"),each=18),
                group2=rep(c("ZeppL-LINE1 region","Upper centromere size estimate","Lower centromere size estimate"),each=18))

df[3,4] <- "Neocentromere on Chr_02"
df[21,4] <- "Neocentromere on Chr_02"
df[39,4] <- "Neocentromere on Chr_02"

df$group2 = factor(df$group2
                   ,levels = c("ZeppL-LINE1 region","Upper centromere size estimate","Lower centromere size estimate","Neocentromere on Chr_02"))

png(file="UL1690_size_scatter_plot.png", width = 13000, height = 4000, units = "px", res=700,pointsize = 20,bg = "white")
ggplot(filter(df,chr!="chr_02_2"),aes(x=group,y=length/1000,col = group2)) +
  geom_jitter(size=5) +
  geom_point(data=filter(df,chr=="chr_02_2"),aes(x=group,y=length/1000),size=5)+
  theme_classic()+
  scale_color_brewer(palette = "Set1",name="")+
  theme(text = element_text(size = 20),legend.text = element_text(size = 25))+
  xlab('ZeppL-LINE1 region size') + ylab('Upper centromere size estimate')
dev.off()


df_max <- data.frame(Chromosome = rep(chr_name,times=1),
                     ZeppL_length = ZeppL_length,
                     Cen_max=Max_length)
df_max[3,1] <- "chr_02(Neocentromere)"
df_max$Chromosome = factor(df_max$Chromosome
                           ,levels = c("chr_01", "chr_02","chr_02(Neocentromere)","chr_03", "chr_04","chr_05", "chr_06","chr_07", "chr_08","chr_09", "chr_10","chr_11", "chr_12","chr_13", "chr_14","chr_15", "chr_16","chr_17"))


library(ggpmisc)

my.formula <- y ~ x

png(file="FigureS5a_scatter_plot.png", width =11000, height = 7000, units = "px", res=700,bg = "white")
df_max %>% ggplot(aes(x=ZeppL_length/1000,y=Cen_max/1000)) + 
  geom_jitter(aes(color=Chromosome),size=9) + 
  geom_smooth( method = "lm",se=FALSE,col="black",formula=my.formula) + 
  scale_colour_manual(name="Chromosome",
                      values=c("#c10023","#008e17","#000000","#fb8500","#f60000", "#FE0092", "#bc9000", 
                               "#4ffc00","#00bcac", "#0099cc", "#D35400", "#00eefd", "#cf6bd6", 
                               "#99cc00", "#aa00ff", "#ff00ff", "#00896e", "#f2a287" ))+    
  theme_bw() + 
  theme(text = element_text(size = 20),legend.text = element_text(size = 25))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  xlab('ZeppL-LINE1 region size') + ylab('Upper centromere size estimate')

dev.off()


df_min <- data.frame(Chromosome = rep(chr_name,times=1),
                     ZeppL_length = ZeppL_length,
                     Cen_min=Min_length)
df_min[3,1] <- "chr_02(Neocentromere)"
df_min$Chromosome = factor(df_min$Chromosome
                           ,levels = c("chr_01", "chr_02","chr_02(Neocentromere)","chr_03", "chr_04","chr_05", "chr_06","chr_07", "chr_08","chr_09", "chr_10","chr_11", "chr_12","chr_13", "chr_14","chr_15", "chr_16","chr_17"))

png(file="FigureS5b_scatter_plot.png", width = 11000, height = 7000, units = "px", res=700,bg = "white")
df_min %>% ggplot(aes(x=ZeppL_length/1000,y=Cen_min/1000)) + 
  geom_jitter(aes(color=Chromosome),size=9) + 
  scale_colour_manual(name="Chromosome",
                      values=c("#c10023","#008e17","#000000","#fb8500","#f60000", "#FE0092", "#bc9000", 
                               "#4ffc00","#00bcac", "#0099cc", "#D35400", "#00eefd", "#cf6bd6", 
                               "#99cc00", "#aa00ff", "#ff00ff", "#00896e", "#f2a287" )) +   
  #geom_smooth( method = "lm",se=FALSE,col="white") + 
  theme_bw() + 
  theme(text = element_text(size = 20),legend.text = element_text(size = 25))+
  xlab('ZeppL-LINE1 region size') + ylab('Upper centromere size estimate')
dev.off()

ZeppL_info <- read.xlsx('ZeppL_4strains_plot.xlsx', sheetName = "Sheet1")
ZeppL_info$Strain = factor(ZeppL_info$Strain 
                           ,levels = c("UL1690","CC400","CC5816","CC4532"))
png(file="bar_plot.png", width = 8000, height = 10000, units = "px", res=500,bg = "white")
ZeppL_info %>% ggplot(aes(x=chr,y=ZeppL_length/1000,fill=Strain))+
  geom_bar(stat = "identity",position = position_dodge (width = 0.75),width=0.7)+
  theme_bw() + 
  theme(text = element_text(size = 30),legend.position = "top",legend.text = element_text(size = 25))+
  scale_fill_manual(name="",values = c("red", "steelblue","mediumorchid","limegreen"))+ 
  #scale_fill_brewer(palette="Dark2")+
  xlab("") + ylab('')+
  coord_flip()
dev.off()


