#!/bin/bash
#SBATCH --job-name=refs_mm10
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=data/logs/%x_%j.log
# mm10 references for the notebooks: genome, prebuilt bowtie2 index, chrom sizes, blacklist,
# GENCODE vM25 GTF and a TSS bed of protein-coding genes. Writes data/genome/mm10.
#
# Run from the topic folder (bulk-atac-seq), with the atac-cli environment active:
#     mkdir -p data/logs && bash scripts/00_download_refs.sh     (or: sbatch scripts/00_download_refs.sh)
set -euo pipefail
G=$PWD/data/genome/mm10
mkdir -p $G && cd $G

wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
[ -s mm10.fa ] || gunzip -k mm10.fa.gz
samtools faidx mm10.fa

wget -q -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
zcat mm10-blacklist.v2.bed.gz | cut -f1-3 > mm10-blacklist.v2.bed

wget -q -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
# TSS of protein-coding genes (one per gene): chrom, tss, tss+1, gene_name, ., strand
zcat gencode.vM25.annotation.gtf.gz | awk -F'\t' '$3=="gene" && $9 ~ /gene_type "protein_coding"/ {
    match($9, /gene_name "[^"]+"/); g = substr($9, RSTART+11, RLENGTH-12)
    tss = ($7 == "+") ? $4 - 1 : $5 - 1
    print $1"\t"tss"\t"tss+1"\t"g"\t.\t"$7 }' | sort -k1,1 -k2,2n > mm10_tss_protein_coding.bed

# prebuilt bowtie2 index (UCSC mm10, Langmead lab)
mkdir -p bowtie2 && cd bowtie2
wget -q -c https://genome-idx.s3.amazonaws.com/bt/mm10.zip
[ -s mm10.1.bt2 ] || unzip -o mm10.zip
ls -la
echo DONE
