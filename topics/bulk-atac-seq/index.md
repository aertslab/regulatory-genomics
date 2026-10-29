---
downloads:
  - file: scripts/00_download_refs.sh
    title: 00_download_refs.sh
    filename: 00_download_refs.sh
  - file: scripts/02_download_bulk.sh
    title: 02_download_bulk.sh
    filename: 02_download_bulk.sh
  - file: scripts/03_rna_counts.sh
    title: 03_rna_counts.sh
    filename: 03_rna_counts.sh
  - file: scripts/10_bulk_pipeline.sh
    title: 10_bulk_pipeline.sh
    filename: 10_bulk_pipeline.sh
  - file: scripts/11_bulk_summary.py
    title: 11_bulk_summary.py
    filename: 11_bulk_summary.py
---

# Bulk ATAC-seq

In this practical we analyse the ATAC-seq data of:

> Hansen JL, Loell KJ, Cohen BA. *A test of the pioneer factor hypothesis using ectopic liver gene activation.* eLife (2022) 11:e73358. https://doi.org/10.7554/eLife.73358
> GEO: [GSE182188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182188)

This is the same study as the [bulk RNA-seq practical](../bulk-rna-seq/index.md). The authors expressed the liver transcription factors FOXA1 and HNF4A in K562 cells, a leukemia cell line in which liver genes are silent, by adding doxycycline (dox). We use the ATAC-seq of the line that expresses both factors: 2 uninduced and 2 induced samples (paired-end, 2 × 75 bp, 90-115 million read pairs per sample). FOXA1 is the textbook example of a pioneer factor, a transcription factor that can open closed chromatin; we ask which regions open (and close) when FOXA1 and HNF4A are switched on, which motifs they contain, and whether they lie near the genes that the RNA-seq practical found activated.

**Research questions:**

1. How good are the ATAC-seq libraries?
2. Which regions open (or close) when FOXA1 and HNF4A are expressed?
3. Do the opened regions contain FOXA1 and HNF4A motifs, and do they lie near the genes that are switched on?

The practical consists of three notebooks:

1. [Part 1: from reads to peaks](part_1_reads_to_peaks.ipynb) (bash): read quality (FastQC), trimming (fastp), alignment (bowtie2), duplicate marking (Picard), filtering, ATAC-specific quality control (fragment sizes, TSS enrichment, FRiP), coverage tracks (deepTools) and peak calling (MACS3), on one uninduced and one induced sample.
2. [Part 2: consensus peaks and MACS3 differential peaks](part_2_consensus_peaks_macs_diff.ipynb) (bash): peaks for all four samples, three ways to build a consensus peak set, a count matrix (featureCounts), and MACS3 `bdgdiff`.
3. [Part 3: DiffBind and DESeq2](part_3_diffbind_deseq2.ipynb) (R): differential accessibility with replicates, the effect of normalisation when many regions open, comparison of the methods, annotation of the regions to genes (ChIPseeker), the link with the activated genes of the RNA-seq practical, motif enrichment (HOMER), a comparison with the paper and the conclusion.

To keep the notebooks fast, they use the reads on **chromosome 17** only; the genome-wide results of all four samples were made with the same steps.

## Software

The software for each notebook is installed in a separate mamba environment. The commands to create them and to register them as Jupyter kernels are given at the start of parts 1 and 3.

| Notebook | Environment | Kernel |
|---|---|---|
| Parts 1 and 2 | `atac-cli` | Bash (atac-cli) |
| Part 3 | `atac-r` (includes HOMER) | R (atac-r) |

## Data

You work in the folder `data`, next to the notebooks: each notebook writes everything it makes to `data/part_N`. The input data (the **course data**) are already prepared for you on the teaching server, in the read-only folder `/data/teachers/data/bulk_atac_seq/hansen_foxa1_hnf4a`; the notebooks find them automatically, so you don't need to copy or download anything. You need about 3.5 GB of disk space in `data`.

If you run the notebooks somewhere else, they look for the course data in a folder `course_data` next to the notebooks; see below.

### Preparing the course data yourself

To run the notebooks somewhere else than on the teaching server (for example on your own computer or on the VSC), build the course data yourself, in a folder `course_data` next to the notebooks. The scripts put the downloads and intermediate files (raw reads, full-genome alignments, about 140 GB) in a second folder, `course_data_build`, which the notebooks don't need. Download the scripts below from the download menu at the top of this page, put them in a folder `scripts` next to the notebooks, and run them from the folder of the notebooks with the `atac-cli` environment (each script has a Slurm header, and the usage in its first lines):

The scripts read the sample sheet from `course_data_build/samplesheet.tsv`: run the sample-sheet cell of part 1 (Step 0) and copy `data/part_1/samplesheet.tsv` to `course_data_build/`.

| Script | Does | Time |
|---|---|---|
| `00_download_refs.sh` | chr17 sequence, bowtie2 index (GRCh38 without alt contigs), chrom sizes, ENCODE blacklist, TSS of protein-coding genes from the GENCODE v46 annotation (`course_data/genome/hg38`) | ~1 h |
| `02_download_bulk.sh` | the 4 FASTQ pairs from ENA, md5-checked (`course_data_build/raw`, 47 GB) | ~1 h |
| `03_rna_counts.sh` | RNA-seq gene counts of the same line, from the Salmon quantifications on GEO (`course_data/rna`, used in part 3) | minutes |
| `10_bulk_pipeline.sh` | the steps of part 1 on the full data of each sample (`course_data_build/pipeline/<sample>`), plus the filtered chr17 BAM files (`course_data/data_generated/<sample>`) | ~2.5 h per sample (16 CPU) |
| `11_bulk_summary.py` | genome-wide QC table, chr17 reads of the two part-1 samples (`course_data/fastq_chr17`) | minutes |

The pre-generated results `course_data/nb01_generated` and `course_data/nb02_generated` are only needed to run part 2 or 3 without running the previous part first; they are copies of what parts 1 and 2 write to `data/part_1` (the sorted BAM files) and `data/part_2` (`peaks/`, `consensus_peaks.bed`/`.saf`, `counts_matrix.tsv`, `bdgdiff/`). Rebuilding the course data from the raw reads needs about 150 GB; afterwards `course_data_build` can be deleted.
