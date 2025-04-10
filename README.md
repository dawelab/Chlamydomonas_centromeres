Chlamydomonas UL1690 Genome Analysis Pipeline

🧬 Genome Assembly

  1) Laboratory strains UL-1690 (a derivative of CC-1690/21gr) and CC-400 (cw15) were selected for PacBio HiFi sequencing.
  2) Genome assemblies were generated using HiFiasm.
  3) RagTag was used to correct assembly errors and generate pseudochromosomes.
  
 
🧬 Gene Annotation

  1) Gene models in UL-1690 were annotated based on the GFF3 file of reference strain CC-4532.
  2) Liftoff v1.6.3 was used to transfer annotations from CC-4532 to the UL-1690 assembly.
  3) The longest isoform per gene was retained using AGAT, supporting visualization in IGV and gene density calculation.
 
🧬 Tandem Repeats Identification

   TRF (Tandem Repeats Finder) was used to annotate tandem repeats across four laboratory strains.
     
 
📌 ZeppL-LINE1 Enrichment Analysis

  1) ZeppL_Crei.fa was used as the query for BLASTN searches to identify ZeppL fragments in all genomes.
  2) Fragments < 300 bp and with < 90% identity were filtered out.
  3) ZeppL-enriched regions along chromosomes were identified using BEDTools, with 40 kb intervals.
 
📌 Centromere Estimation

  1) CENH3 CUT&Tag reads were processed as follows:
     1) Adapter trimming: TrimGalore v0.6.7.
     2) Genome alignment: Bowtie2 v2.4.5.
  2) Mapping strategies:
     1) Multiple + unique mapping: used to estimate upper bounds of CENH3 enrichment.
     2) Unique-only mapping (MAPQ ≥ 20): used to estimate lower bounds of CENH3 enrichment.
  3) Normalization & visualization:
     1) BAMs converted to bedgraph and BED using BEDTools.
     2) Per-bin (5 kb) read depth and fold enrichment computed and visualized with karyoploteR in R.

📌 Telomere Identification

  'CCCTAAAACCCTAAAACCCTAAAACCCTAAAA' was used for BLAST to identify sing reads with telomere sequences. 
        
📊 R Code & Plotting
   1) All R scripts for generating plots are included.
   2) All key input data files are provided.
