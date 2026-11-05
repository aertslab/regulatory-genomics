# Bulk ATAC-seq: preparing the course data

Notes for teachers: how the data of the notebooks are built and installed on the teaching server. Students don't need this; the topic page (`../index.md`) tells them where the data are.

## Layout

Three folders next to the notebooks, all git-ignored:

- `course_data/`: the **course data**, everything the notebooks read and nothing else (about 7.5 GB). On the teaching server the notebooks read it from the shared, read-only folder `/data/teachers/data/bulk_atac_seq/hansen_foxa1_hnf4a`; when that folder doesn't exist, they use `course_data/` next to the notebooks. The setup cell of each notebook chooses between the two (variable `BULK`, or `bulk` in part 3).
- `course_data_build/`: downloads and intermediate files, only needed to (re)build `course_data/` (about 150 GB): the sample sheet, raw reads, the full-genome pipeline outputs, the hg38 FASTA and GENCODE annotation, the RNA-seq quantifications, logs.
- `data/`: the working folder of the notebooks; each notebook writes to `data/part_N` (one per student, writable).

## Installing on the teaching server

Copy `course_data/` as a whole:

```bash
rsync -a course_data/ /data/teachers/data/bulk_atac_seq/hansen_foxa1_hnf4a/
```

What it contains:

| Path | Contents | Size |
|---|---|---|
| `genome/hg38/bowtie2/` | bowtie2 index of GRCh38 without alt contigs | 4.9 GB |
| `genome/hg38/chr17.fa`, `chr17.fa.fai`, `hg38.chrom.sizes`, `hg38-blacklist.v2.bed`, `hg38_tss_protein_coding.bed` | chr17 sequence (all reads and peaks of the notebooks are on chr17), chromosome sizes, ENCODE blacklist, TSS of protein-coding genes | 82 MB |
| `fastq_chr17/` | chr17 reads of the two samples of part 1 | 1 GB |
| `data_generated/qc_summary.tsv`, `data_generated/<sample>/` | genome-wide QC table; per sample the idxstats (chrM fraction) and the filtered chr17 BAM file | 1.1 GB |
| `nb01_generated/`, `nb02_generated/` | pre-generated results of parts 1 and 2, so that each notebook can also run on its own | 0.4 GB |
| `rna/foxa1hnf4a_gene_counts.tsv` | RNA-seq gene counts of the same line (part 3) | 1.4 MB |

## Rebuilding the course data

The steps are on the topic page (`../index.md`, section "Preparing the course data yourself"), so that students can also follow them. The scripts write the files the notebooks use to `course_data/` and everything else to `course_data_build/`. `nb01_generated/` and `nb02_generated/` are not made by the scripts: after running parts 1 and 2, copy the sorted BAM files of `data/part_1` and the `peaks/`, `consensus_peaks.bed`/`.saf`, `counts_matrix.tsv` and `bdgdiff/` of `data/part_2`.

The RNA counts table in `course_data/rna` was made before `03_rna_counts.sh`; the script reproduces it to within ±1 count for a few genes (rounding) and in a different row order.
