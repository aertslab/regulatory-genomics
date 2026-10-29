#!/bin/bash
#SBATCH --job-name=refs_hg38
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# hg38 references for the notebooks: the chr17 sequence (all regions in the notebooks are on chr17), prebuilt bowtie2 index
# (GRCh38 without alt contigs, UCSC names), chrom sizes, blacklist and a TSS bed of protein-coding genes (from the GENCODE v46
# GTF, as in the bulk RNA-seq practical). The files the notebooks use go to course_data/genome/hg38; the downloaded archives
# stay in course_data_build/genome/hg38.
#
# Run from the topic folder (bulk-atac-seq), with the atac-cli environment active:
#     mkdir -p course_data_build/logs && bash scripts/00_download_refs.sh     (or: sbatch scripts/00_download_refs.sh)
set -euo pipefail
G=$PWD/course_data/genome/hg38        # used by the notebooks
X=$PWD/course_data_build/genome/hg38  # downloads only needed to make them
mkdir -p $G $X && cd $X

wget -q -c -P $G https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
[ -s $G/chr17.fa ] || zcat chr17.fa.gz > $G/chr17.fa
samtools faidx $G/chr17.fa

wget -q -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
zcat hg38-blacklist.v2.bed.gz | cut -f1-3 > $G/hg38-blacklist.v2.bed

wget -q -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
# TSS of protein-coding genes (one per gene): chrom, tss, tss+1, gene_name, ., strand
zcat gencode.v46.primary_assembly.annotation.gtf.gz | awk -F'\t' '$3=="gene" && $9 ~ /gene_type "protein_coding"/ {
    match($9, /gene_name "[^"]+"/); g = substr($9, RSTART+11, RLENGTH-12)
    tss = ($7 == "+") ? $4 - 1 : $5 - 1
    print $1"\t"tss"\t"tss+1"\t"g"\t.\t"$7 }' | sort -k1,1 -k2,2n > $G/hg38_tss_protein_coding.bed

# prebuilt bowtie2 index (Langmead lab); the zip unpacks the index files directly into the bowtie2 folder
mkdir -p $G/bowtie2 && cd $G/bowtie2
[ -s GRCh38_noalt_as.1.bt2 ] || wget -q -c https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
if [ ! -s GRCh38_noalt_as.1.bt2 ]; then unzip -o GRCh38_noalt_as.zip && rm GRCh38_noalt_as.zip; fi
ls -la $G $G/bowtie2
echo DONE
