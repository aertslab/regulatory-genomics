# PWM and CRM scanning

In this practical we ask how well a transcription factor motif, written as a position weight matrix (PWM), predicts where that factor actually binds in a genome. We start with the classic example of GAL4 in yeast and then move to CEBPA, a master regulator of the liver, on human chromosome 19.

Scanning a genome with a single motif runs into the **futility theorem** (Wasserman & Sandelin, *Nat Rev Genet* 2004): almost every predicted site is never bound in vivo. We first quantify this with information theory (R<sub>sequence</sub> vs. R<sub>frequency</sub>) and with FIMO, and compare the hits with the CEBPA ChIP-seq peaks of the [ChIP-seq practical](../chip-seq/index.md). We then look at what distinguishes bound from unbound sites (GC content, nearby partner motifs, clustering of sites), which leads to **cis-regulatory module (CRM) scanning**: looking for dense clusters of motif matches instead of single sites.

The motifs come from [JASPAR](https://jaspar.elixir.no/): GAL4 (MA0299.1), CEBPA (MA0102.3) and the liver partners HNF4A (MA0114.5), FOXA2 (MA0047.4), HNF1A (MA0046.3) and ONECUT1/HNF6 (MA0679.3).

The practical consists of two notebooks:

1. [PWM scanning: how many binding sites does a motif predict?](pwm-scanning.ipynb) (bash): the background base composition and information content of a genome, scanning the yeast genome with GAL4 (FIMO) and comparing the hits with the ~10 known GAL4 targets, the expected number of false positives, scanning human chr19 with CEBPA and intersecting the hits with liver CEBPA ChIP-seq peaks (bedtools), and comparing bound and unbound sites (GC content, partner motifs, clustering).
2. [CRM scanning: from single sites to clusters](crm-scanning.ipynb) (bash): how a CRM scanner like Cluster-Buster works, homotypic clusters (CEBPA only) and heterotypic clusters (CEBPA with its liver partners) on chr19, their overlap with the CEBPA peaks compared with randomly placed regions, which motifs make up the clusters, and a control with scrambled motifs.

Run the PWM notebook first: the CRM notebook reuses files it creates (for example `cebpa_peaks.bed`).

## Software

Both notebooks run in a Bash kernel. On the teaching server the tools are loaded as modules at the start of the notebook:

    module load MEME/5.5.5-gompi-2023a
    module load SAMtools/1.18-GCC-12.3.0
    module load BEDTools/2.31.0-GCC-12.3.0

The PWM notebook also calls `Rscript` for a Wilcoxon test. The CRM notebook downloads and compiles [Cluster-Buster](https://github.com/weng-lab/cluster-buster) itself (`git` and `make` are required).

| Notebook | Environment | Kernel |
|---|---|---|
| PWM scanning | modules (MEME, SAMtools, BEDTools, R) | Bash (PWM) |
| CRM scanning | modules (SAMtools, BEDTools) + Cluster-Buster | Bash (PWM) |

## Data

The yeast genome (sacCer3), its gene annotation (`sgdGene`) and the JASPAR matrices are downloaded by the notebooks themselves.

The human part needs two files from earlier steps of the course:

* the human reference genome `../genome/hg38.fa` (only chromosome 19 is extracted from it);
* the filtered CEBPA peaks `../res/CEBPA_peaks.filtered.narrowPeak`, produced in the [ChIP-seq practical](../chip-seq/index.md).
