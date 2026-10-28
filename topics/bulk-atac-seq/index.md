# Bulk ATAC-seq

In this practical we analyse the liver ATAC-seq data of:

> Jeong BK, Choi WI, Choi W, et al. *A male mouse model for metabolic dysfunction-associated steatotic liver disease and hepatocellular carcinoma.* Nature Communications (2024) 15:6506. https://doi.org/10.1038/s41467-024-50660-y
> GEO: [GSE246213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246213)

The authors built a mouse model of metabolic dysfunction-associated steatotic liver disease (MASLD) that progresses to liver cancer. Male mice received a low dose of streptozotocin (STZ, which destroys the insulin-producing beta cells) and/or a high-fat diet (HFD) instead of standard chow (SCD), in a 2 × 2 design with 3 mice per group. ATAC-seq was done on liver tissue at 20 weeks of age (paired-end, 2 × 151 bp, about 44 million read pairs per sample). We ask which regions change accessibility in the disease model, and how much the answer depends on library quality, normalisation and replicates.

**Research questions:**

1. How good are the ATAC-seq libraries, and how does their quality shape the results?
2. Which regions change accessibility in the livers of STZ + HFD mice, compared with control mice?
3. Which genes lie near these regions?

The practical consists of three notebooks:

1. [Part 1: from reads to peaks](part_1_reads_to_peaks.ipynb) (bash): read quality (FastQC), trimming (fastp), alignment (bowtie2), duplicate marking (Picard), filtering, ATAC-specific quality control (fragment sizes, TSS enrichment, FRiP), coverage tracks (deepTools) and peak calling (MACS3), on one control and one STZ + HFD mouse.
2. [Part 2: consensus peaks and MACS3 differential peaks](part_2_consensus_peaks_macs_diff.ipynb) (bash): peaks for all 12 mice, three ways to build a consensus peak set, a count matrix (featureCounts), and MACS3 `bdgdiff`.
3. [Part 3: DiffBind and DESeq2](part_3_diffbind_deseq2.ipynb) (R): differential accessibility with replicates, the effect of normalisation, a factorial design (STZ × diet), comparison of the methods, annotation of the regions to genes (ChIPseeker), a comparison with the paper and the conclusion.

To keep the notebooks fast, they use the reads on **chromosome 19** only; the genome-wide results of all 12 samples were made with the same steps.

## Software

The software for each notebook is installed in a separate mamba environment, described in the folder `envs`. The commands to create them and to register them as Jupyter kernels are given at the start of parts 1 and 3.

| Notebook | Environment | Kernel |
|---|---|---|
| Parts 1 and 2 | `atac-cli` (`envs/atac-cli.yml`) | Bash (atac-cli) |
| Part 3 | `atac-r` (`envs/atac-r.yml`) | R (atac-r) |

## Data

All data is stored in the folder `data`, next to the notebooks; the notebooks write their results to `data/results/part_N`. The notebooks show the download commands as text (so they are not run by accident); the data are built by the scripts in `scripts/`, run from this folder with the `atac-cli` environment (each script has a Slurm header, and the usage in its first lines):

| Script | Does | Time |
|---|---|---|
| `00_download_refs.sh` | mm10 genome, bowtie2 index, chrom sizes, ENCODE blacklist, GENCODE vM25 annotation, TSS of protein-coding genes (`data/genome/mm10`) | ~1 h |
| `02_download_bulk.sh` | the 12 FASTQ pairs from ENA, md5-checked (`data/raw`, 92 GB) | 1-2 h |
| `10_bulk_pipeline.sh` | the steps of part 1 on the full data of each sample, plus the chr19 subsets (`data/data_generated/<sample>`) | ~1.5 h per sample (16 CPU) |
| `11_bulk_summary.py` | genome-wide QC table, chr19 reads of the two part-1 samples (`data/fastq_chr19`) | minutes |

The notebooks themselves need about 10 GB: the genome and bowtie2 index (7 GB), the chr19 reads and BAM files, and the pre-generated results of parts 1 and 2 (`data/nb01_generated`, `data/nb02_generated`), which let each notebook run on its own. Rebuilding everything from the raw reads needs about 220 GB of disk space.
