# Nascent_RNA_manuscript
Manuscript code repository This repository contains the commands and data related to EU_nuclear RNA published manuscript

### Samples list

https://dataview.ncbi.nlm.nih.gov/object/PRJNA978065

Accession  | Sample | Fw | Rv
--- | --- | --- | --- 
SRR24776115 | EU_nascent_BR1 | BR1_EU_nascent_S2_L004_R1_001.fastq.gz | BR1_EU_nascent_S2_L004_R2_001.fastq.gz
SRR24776114 | EU_nascent_BR2 | BR2_EU_nascent_S10_L004_R1_001.fastq.gz | BR2_EU_nascent_S10_L004_R2_001.fastq.gz
SRR24776113 | EU_nascent_BR3 | BR3_EU_nascent_S4_L004_R1_001.fastq.gz | BR3_EU_nascent_S4_L004_R2_001.fastq.gz
SRR24776112 | EU_nascent_BR4 | BR4_EU_nascent_S12_L004_R1_001.fastq.gz | BR4_EU_nascent_S12_L004_R2_001.fastq.gz
SRR24776111 | EU_nascent_BR5 | BR5_EU_nascent_S17_L004_R1_001.fastq.gz | BR5_EU_nascent_S17_L004_R2_001.fastq.gz
SRR24776110 | Nuclear_BR1 | BR1_EU_contol_S1_L004_R1_001.fastq.gz | BR1_EU_contol_S1_L004_R2_001.fastq.gz
SRR24776109 | Nuclear_BR2 | BR2_EU_contol_S9_L004_R1_001.fastq.gz | BR2_EU_contol_S9_L004_R2_001.fastq.gz
SRR24776108 | Nuclear_BR3 | BR3_EU_contol_S3_L004_R1_001.fastq.gz | BR3_EU_contol_S3_L004_R2_001.fastq.gz
SRR24776107 | Nuclear_BR4 | BR4_EU_contol_S11_L004_R1_001.fastq.gz | BR4_EU_contol_S11_L004_R2_001.fastq.gz
SRR24776106 | Nuclear_BR5 | BR5_EU_contol_S13_L004_R1_001.fastq.gz | BR5_EU_contol_S13_L004_R2_001.fastq.gz

https://www.ncbi.nlm.nih.gov/bioproject/PRJNA853797

Accession  | Sample | Fw | Rv
--- | --- | --- | --- 
SRR19889573	| Z.mays_B73_Primary+SeminalRoot_0-1mm_steady-state_BR4 | | 
SRR19889574	| Z.mays_B73_Primary+SeminalRoot_0-1mm_steady-state_BR3 | | 
SRR19889575	| Z.mays_B73_Primary+SeminalRoot_0-1mm_steady-state_BR2 | | 
SRR19889576	| Z.mays_B73_Primary+SeminalRoot_0-1mm_steady-state_BR1 | | 

### Step1 - trim reads
Raw FASTQ files were preprocessed with Trimmomatic v0.39 (Bolger et al., 2014) to remove Illumina sequencing adapters. 5′ and 3′ ends with a quality score below 5 (Phred+33) were trimmed and reads shorter than 20 bp after trimming were discarded. 

### Step2 - sortMeRNA
The reads matching sequences of RNA-Seq contaminants were filtered using sortMeRNA (Kopylova et al., 2012), which compares the reads against a custom database of known contaminants, including ribosomal subunits 28S, 18S, 5.8S and 5S, Internal Transcribed Spacer (ITS), Signal recognition particle (SRP), and chloroplast and mitochondrial genomes. The database was compiled as described on our Github (https://github.com/lconcia/RiboScreen/). 

### Step3 - align reads
The remaining reads were aligned to the Zm-B73-NAM-5.0 reference genome assembly (Hufford et al., 2021) with STAR (Dobin et al., 2013), limiting the output to the uniquely-mapping reads only. The resulting bam files were sorted and indexed, and reads mapped in proper pairs (with SAM flag 0x2) were selected using samtools (Danecek et al., 2021).

### Step4 - deduplication
 Duplicated reads were removed with UMI-Tools v1.1.4 (Smith et al., 2017).
 
### Step5
Read normalization across the genome
Normalized read density over non-overlapping 50-bp bins was calculated using the normalization reads per genomic content (RPGC) in the utility bamCoverage from the bioinformatics suite deeptools2 v3.5.4 (Ramirez et al., 2016).

### Step6
Meta-gene profiles of RNA-Seq signal were generated with deeptools2 computeMatrix and plotted with plotProfile.

### Step7 Read counts in annotations
Raw read counts over exons were obtained using htseq-counts v. 2.0.3 (Putri et al., 2022) and the NCBI RefSeq v5 Zea mays release 103 gene annotation (https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Zea_mays/103/). Reads mapping on multiple genes were assigned fractionally to each gene. To determine the distribution of reads among introns, exons, intergenic space and junctions between these three, we ran htseq-counts using a custom annotation file that included annotations of exons, introns and intergenic spaces. We applied the setting "--nonunique none" to report as "junctions" all the reads overlapping with more than one feature, such as adjacent introns and exons. 

### Step8 - Biological replicate reproducibility and differential expression analysis
Differential expression analysis and principal component analysis were performed with the package DESeq2 v1.40.2 (Love et al., 2014) for the R statistical suite **(RCoreTeam, 2024). All charts were plotted with the R package ggplot2 v3.5.1 (Wickham, 2016). The volcano plots in Figure 5 were drawn with the R package EnhancedVolcano v1.18.0 (Blighe et al., 2024).

### Step 9  - Saturation analysis
To determine the saturating sequencing depth in each RNA type, we took a similar approach as **Tarazona et al, (2011). After pooling all the biological replicates of each RNA type, we randomly sub-sampled the reads at increasing depths between 2 million and 75 million reads, and scored the number of genes with at least 5 read counts in each sub-sample. Using the same data, we established, for each RNA type, the number of newly detected genes at each increment of sequencing depth.

### Step10  - TPM and intron/exons - **work in progress
TPM (transcripts per million) were calculated as in **Wagner et al, 2012.
For each sample, we divided the number of reads aligned on the exons of each gene by the length of the gene, to obtain the "reads per kilobase" (RPK) of each gene. We then divided the sum of all RPKs by 1,000,000, to have a sample-specific "scaling factor" that accounts for the library size. Finally, we divided the RPKs values by the scaling factor to obtain the TPM of each gene.

### Step11
Replication timing analysis
Maize B73 replication timing (Repli-seq) data including G1 and early, mid and late S-phase datasets from sorted nuclei in the mitotic cell cycle of terminal 1-mm root tips (NCBI SRA # PRJNA327875; Wear et al., 2017) were analyzed using the Repliscan application (Zynda et al., 2017; Mickelson-Young et al., 2022)** implemented on the CyVerse Discovery Environment platform (Swetnam et al., 2024). The recommended maize B73 parameters from Mickelson-Young et al. (2024)** were used, except that the full datasets (not downsampled) were mapped to Zm-B73-NAM-5.0 genome assembly (Hufford et al., 2021), and the following minor modifications were made. Read densities were aggregated in 10-kb non-overlapping windows across the genome (--window 10000). Each 10-kb window was assigned a predominant replication time (RT) class using a segmentation threshold of 1.0 (--threshold value, --value 1.0), so that the segmentation algorithm included any window with early, mid or late replication signals that are greater than the G1 reference signal (> ratio of 1.0). 
