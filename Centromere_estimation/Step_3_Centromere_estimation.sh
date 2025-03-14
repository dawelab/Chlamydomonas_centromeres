## Download from R output fi_peak_1_20_0.txt fi_peak_1_q20_20_0.txt
## Calculate the centromere size estimation as code below
DATADIR="scratch/mw36149/Chlamydomonas_centromere"
mkdir $DATADIR/centromere_estimation && cd $DATADIR/centromere_estimation
# File: fi_peak_1_20_0.txt  fi_peak_1_q20_20_0.txt
# Upper centromere estimation
awk '$4>5.5' $DATADIR/fi_peak_1_20_0.txt > enrichment_fold_5.5_noq20.txt
bedtools merge -d 30000 -i enrichment_fold_5.5_noq20.txt | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' > merge_enrichment_5.5_30kb_noq20.txt
awk '$4>10000' merge_enrichment_5.5_30kb_noq20.txt | awk '{print $1"\t"$2"\t"$3}' > centromere_30kb_5.5_noq20.bed
bedtools intersect -wa -a centromere_30kb_5.5_noq20.bed -b cc1690_ZeppL_v2.bed > cc1690_centromere_upper.bed
# lower centromere estimation
awk '$4>10' $DATADIR/fi_peak_1_q20_20_0.txt > enrichment_fold_10.txt
bedtools merge -d 30000 -i enrichment_fold_10.txt | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' > merge_enrichment_10_30kb_q20.txt
awk '$4>10000' merge_enrichment_10_30kb_q20.txt | awk '{print $1"\t"$2"\t"$3}' > centromere_30kb_10_q20.bed
bedtools intersect -wa -a centromere_30kb_10_q20.bed -b cc1690_ZeppL_v2.bed > cc1690_centromere_lower.bed

#Manually add secondary centromere on chr2 tail to files cc1690_centromere_upper.bed and cc1690_centromere_lower.bed
echo 'chr_02\t8430000\t8460000' >> cc1690_centromere_upper.bed 
echo 'chr_02\t8425000\t8460000' >> cc1690_centromere_lower.bed
cp cc1690_centromere_upper.bed  $DATADIR
cp cc1690_centromere_lower.bed  $DATADIR 

