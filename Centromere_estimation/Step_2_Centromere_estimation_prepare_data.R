# Load required packages
library(tidyverse)

# Set working directory
setwd('/Users/mingyu/Desktop/cc1690_HiFi/centromere_estimation')

# Function to process coverage data
process_coverage <- function(cenh3_file, IgG_file) {
  cenh3 <- read.table(cenh3_file, sep = "\t", fill = TRUE, col.names = c('chr', 'start_pos', 'end_pos', 'CenH3'))
  IgG <- read.table(IgG_file, sep = "\t", fill = TRUE, col.names = c('chr', 'start_pos', 'end_pos', 'IgG'))

  # Adjust background noise
  IgG$IgG <- ifelse(IgG$IgG <= 10, 20, IgG$IgG)

  # Normalize values
  cenh3 <- cenh3 %>% mutate(CenH3_norm = CenH3 / sum(CenH3))
  IgG <- IgG %>% mutate(IgG_norm = IgG / sum(IgG))

  # Merge datasets and compute ratio
  df <- left_join(cenh3, IgG, by = c("chr", "start_pos", "end_pos")) %>%
        mutate(ratio = CenH3_norm / IgG_norm) %>%
        filter(IgG > 0) %>%
        select(chr, start_pos, end_pos, ratio)

  return(df)
}

# Process datasets
fi_peak_1_0 <- process_coverage('A_CenH3_1.coverage.5kb_noq20.bed', 'A_IgG.coverage.5kb_noq20.bed')
fi_peak_1_q20_0 <- process_coverage('A_CenH3_1.coverage.5kb.q20.bed', 'A_IgG.coverage.5kb.q20.bed')

# Save output
write.table(fi_peak_1_0, file = "fi_peak_1_20_0.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(fi_peak_1_q20_0, file = "fi_peak_1_q20_20_0.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Function for plotting
plot_peaks <- function(data, output_file, ylim_max) {
  png(output_file, height = 1000, width = 1000)
  data %>%
    ggplot(aes(start_pos / 1000, ratio)) +
    geom_area(aes(fill = chr)) +
    ylim(0, ylim_max) +
    theme_bw() +
    xlab('Genome position (kb)') + ylab('Ratio') +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 20)) +
    facet_wrap(~chr, ncol = 1, strip.position = "right") +
    theme(strip.text.y = element_text(size = 13, angle = 270))
  dev.off()
}

# Generate plots only for checking
plot_peaks(fi_peak_1_0, 'peaks_5kb_20_0.png', 115)
plot_peaks(fi_peak_1_q20_0, 'q20_peaks_5kb_20_0.png', 235)