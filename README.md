Chlamydomonas UL1690 Genome Analysis Pipeline

🧬 Genome Assembly

	• Laboratory strains UL-1690 (a derivative of CC-1690/21gr) and CC-400 (cw15) were selected for PacBio HiFi sequencing.
	• Genome assemblies were generated using HiFiasm.
	• RagTag was used to correct assembly errors and generate pseudochromosomes.
 
🧬 Gene Annotation

	• Gene models in UL-1690 were annotated based on the GFF3 file of reference strain CC-4532.
	• Liftoff v1.6.3 was used to transfer annotations from CC-4532 to the UL-1690 assembly.
	• The longest isoform per gene was retained using AGAT, supporting visualization in IGV and gene density calculation.
 
🧬 Tandem Repeats Identification

	• TRF (Tandem Repeats Finder) was used to annotate tandem repeats across four laboratory strains.
 
📌 ZeppL-LINE1 Enrichment Analysis

	• ZeppL_Crei.fa was used as the query for BLASTN searches to identify ZeppL fragments in all genomes.
	• Fragments < 300 bp and with < 90% identity were filtered out.
	• ZeppL-enriched regions along chromosomes were identified using BEDTools, with 40 kb intervals.
 
📌 Centromere Estimation

	• CENH3 CUT&Tag reads were processed as follows:
	      Adapter trimming: TrimGalore v0.6.7.
	      Genome alignment: Bowtie2 v2.4.5.
	• Mapping strategies:
	      Multiple + unique mapping: used to estimate upper bounds of CENH3 enrichment.
	      Unique-only mapping (MAPQ ≥ 20): used to estimate lower bounds of CENH3 enrichment.
	• Normalization & visualization:
	      BAMs converted to bedgraph and BED using BEDTools.
	      Per-bin (5 kb) read depth and fold enrichment computed and visualized with karyoploteR in R.
       
📊 R Code & Plotting

	• All R scripts for generating plots are included.
	• All key input data files are provided.
