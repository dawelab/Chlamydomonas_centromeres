# Set working directory
setwd('/Users/mingyuwang/Desktop/Chlamydomonas_summary_2024_12_02/Code/Plot_data/FigureS5')

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpmisc)

# Load the data
df_info <- read.table("cc1690_hifi_info.txt", header = FALSE, sep = "\t", fill = TRUE)
colnames(df_info) <- c("chr", "chr_length", "ZeppL_length", "Max_length", "Min_length")

# ---- Figure S5a: ZeppL vs upper centromere ----
df_max <- df_info %>%
  filter(!is.na(Max_length)) %>%
  mutate(Chromosome = ifelse(chr == "chr_02_2", "chr_02(Neocentromere)", chr)) %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome)))

my.formula <- y ~ x

png(file = "FigureS5a_scatter_plot.png", width = 11000, height = 7000, units = "px", res = 700, bg = "white")
ggplot(df_max, aes(x = ZeppL_length / 1000, y = Max_length / 1000)) +
  geom_jitter(aes(color = Chromosome), size = 9) +
  geom_smooth(method = "lm", se = FALSE, col = "black", formula = my.formula) +
  scale_color_manual(
    name = "chr",
    values = c(
      "#c10023", "#008e17", "#000000", "#fb8500", "#f60000", "#FE0092", "#bc9000", 
      "#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#cf6bd6", 
      "#99cc00", "#aa00ff", "#ff00ff", "#00896e", "#f2a287"
    )
  ) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.text = element_text(size = 25)) +
  stat_poly_eq(formula = my.formula, aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE) +
  xlab('ZeppL-LINE1 region size (kb)') + ylab('Upper centromere size estimate (kb)')
dev.off()

# ---- Figure S5b: ZeppL vs lower centromere ----
df_min <- df_info %>%
  filter(!is.na(Min_length)) %>%
  mutate(Chromosome = ifelse(chr == "chr_02_2", "chr_02(Neocentromere)", chr)) %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome)))

png(file = "FigureS5b_scatter_plot.png", width = 11000, height = 7000, units = "px", res = 700, bg = "white")
ggplot(df_min, aes(x = ZeppL_length / 1000, y = Min_length / 1000)) +
  geom_jitter(aes(color = Chromosome), size = 9) +
  scale_color_manual(
    name = "chr",
    values = c(
      "#c10023", "#008e17", "#000000", "#fb8500", "#f60000", "#FE0092", "#bc9000", 
      "#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#cf6bd6", 
      "#99cc00", "#aa00ff", "#ff00ff", "#00896e", "#f2a287"
    )
  ) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.text = element_text(size = 25)) +
  xlab('ZeppL-LINE1 region size (kb)') + ylab('Lower centromere size estimate (kb)')
dev.off()

# ---- Figure S5c: Chromosome size vs ZeppL ----
sorted_levels <- c("chr_01", "chr_02", "chr_03", "chr_04", "chr_05", "chr_06", 
                   "chr_07", "chr_08", "chr_09", "chr_10", "chr_11", "chr_12", "chr_13", 
                   "chr_14", "chr_15", "chr_16", "chr_17")

df_chr <- df_info %>%
  filter(!is.na(chr_length) & !is.na(ZeppL_length)) %>%
  filter(chr != "chr_02_2") %>%
  mutate(
    chr = factor(chr, levels = sorted_levels),
    chr_length_mb = chr_length / 1e6,
    ZeppL_length_kb = ZeppL_length / 1e3
  )

png(file = "FigureS5c_chr_size_vs_ZeppL.png", width = 9000, height = 7000, units = "px", res = 700, bg = "white")
ggplot(df_chr, aes(x = chr_length_mb, y = ZeppL_length_kb, color = chr)) +
  geom_point(size = 9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "black", formula = my.formula) +
  stat_poly_eq(
    formula = my.formula,
    aes(label = paste(..rr.label.., sep = "~~~"), group = 1),
    parse = TRUE,
    size = 6
  ) +
  scale_color_manual(
    name = "chr",
    values = c(
      "#c10023", "#008e17", "#fb8500", "#f60000", "#FE0092", "#bc9000", 
      "#4ffc00", "#00bcac", "#0099cc", "#D35400", "#00eefd", "#cf6bd6", 
      "#99cc00", "#aa00ff", "#ff00ff", "#00896e", "#f2a287"
    )
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = "right",
    legend.text = element_text(size = 20)
  ) +
  xlab("Chromosome size (Mb)") +
  ylab("ZeppL-LINE1 region size (Kb)")
dev.off()


