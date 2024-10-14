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

### Step1 - Adapter trimming
Raw FASTQ files were preprocessed with Trimmomatic v0.39 (Bolger et al., 2014) to remove Illumina sequencing adapters. 5′ and 3′ ends with a quality score below 5 (Phred+33) were trimmed and reads shorter than 20 bp after trimming were discarded. 

```bash
for x in *_R1_001.fastq.gz
do
trimmomatic PE -threads 16 -phred33 -validatePairs -summary $(basename $x _R1_001.fastq.gz).summary.txt \
/$(basename $x _R1_001.fastq.gz)_R1_001.fastq.gz \
/$(basename $x _R1_001.fastq.gz)_R2_001.fastq.gz \
/trim.$(basename $x _R1_001.fastq.gz)_R1_001.fastq \
/trim.$(basename $x _R1_001.fastq.gz)_R1_001.U.fastq \
/trim.$(basename $x _R1_001.fastq.gz)_R2_001.fastq \
/trim.$(basename $x _R1_001.fastq.gz)_R2_001.U.fastq \
ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:5 TRAILING:5 MINLEN:20  
done
```

### Step2 - contaminant filtering with sortMeRNA
The reads matching sequences of RNA-Seq contaminants were filtered using sortMeRNA (Kopylova et al., 2012), by comparing the reads against a custom database of known contaminants, including ribosomal subunits 28S, 18S, 5.8S and 5S, Internal Transcribed Spacer (ITS), Signal recognition particle (SRP), and chloroplast and mitochondrial genomes. The database was compiled as described on Github (https://github.com/lconcia/RiboScreen/).

```bash
for x in trim*_R1_001.fastq
do
sortmerna --threads 16 --out2 True --zip-out Yes --fastx true --blast '1 cigar qcov qstrand' --paired_in true \
-ref sortMeRNA.B73_AGPv5.db/ChlMt_db/Mt.vs_AGPv5.fasta \
-ref sortMeRNA.B73_AGPv5.db/ChlMt_db/Chl.vs_AGPv5.fasta \
-ref sortMeRNA.B73_AGPv5.db/tRNAs_db/tRNAs.vs_AGPv5.merged.fasta \
-ref sortMeRNA.B73_AGPv5.db/rRNA_db/annotation_vs_AGPv5/5.8S.Zea_Mays_vs_AGPv5.querysize.above_90.merged.fasta \
-ref sortMeRNA.B73_AGPv5.db/rRNA.Zm-B73-REFERENCE-NAM-5.0_genomic.with_scaffolds.ITS.fasta \
-ref sortMeRNA.B73_AGPv5.db/18S.Zea_Mays_vs_AGPv5.querysize.above_90.merged.fasta \
-ref sortMeRNA.B73_AGPv5.db/28S.Zea_Mays_vs_AGPv5.querysize.above_90.merged.fasta \
-ref sortMeRNA.B73_AGPv5.db/5S.Zea_Mays_vs_AGPv5.querysize.above_90.merged.fasta \
-ref sortMeRNA.B73_AGPv5.db/SRP_db/SRP.Zea_Mays_vs_AGPv5.querysize.above_90.merged.fasta \
-reads $(basename $x _R1_001.fastq).contaminant.R1.fastq \
-reads $(basename $x _R1_001.fastq).contaminant.R2.fastq \
--idx-dir sortMeRNA.B73_AGPv5.db/sortmerna.idx-dir \
--workdir output_sortMeRNA.$(basename $x _R1_001.fastq) \
--kvdb    output_sortMeRNA.$(basename $x _R1_001.fastq)/kvdb \
--readb   output_sortMeRNA.$(basename $x _R1_001.fastq)/readb \
--other   output_sortMeRNA.$(basename $x _R1_001.fastq)/clean_$(basename $x _R1_001.fastq).fastq
done
```

### Step3 - reads alignment
The surviving reads were aligned to the Zm-B73-NAM-5.0 reference genome assembly (Hufford et al., 2021) with STAR (Dobin et al., 2013), limiting the output to the uniquely-mapping reads only. The resulting bam files were sorted and indexed, and reads mapped in proper pairs (with SAM flag 0x2) were selected using samtools (Danecek et al., 2021).

### Step4 - de-duplication of aligned reads
Duplicated reads were removed with UMI-Tools v1.1.4 (Smith et al., 2017).

### Step5 - Read counts over annotated genes
Raw read counts over exons were obtained using htseq-counts v. 2.0.3 (Putri et al., 2022) and the NCBI RefSeq v5 Zea mays release 103 gene annotation (https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Zea_mays/103/). Reads mapping on multiple genes were assigned fractionally to each gene. To determine the distribution of reads among introns, exons, intergenic space and junctions between these three, we ran htseq-counts using a custom annotation file modified to annotations of exons, introns and intergenic spaces. We applied the setting "--nonunique none" to report as "junctions" all the reads overlapping with more than one feature, such as adjacent introns and exons. 

#################

