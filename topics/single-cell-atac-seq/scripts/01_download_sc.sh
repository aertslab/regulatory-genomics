#!/bin/bash
#SBATCH --job-name=dl_sc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Liver multiome (GSE218468): ATAC fragments of Mouse_4 and Mouse_5, the published count matrices,
# and the few mm10 annotation files the notebooks use. Writes course_data/fragments, course_data/genome and the ATAC count
# matrix in course_data/raw (used by part 1); the original GEO fragment files, the RNA count matrix and Mallet go to course_data_build.
#
# Run from the topic folder (single-cell-atac-seq), with an environment that has wget, bgzip and tabix
# (e.g. atac-cli of the bulk ATAC-seq practical):
#     mkdir -p course_data_build/logs && bash scripts/01_download_sc.sh     (or: sbatch scripts/01_download_sc.sh)
set -euo pipefail
D=$PWD/course_data        # used by the notebooks
X=$PWD/course_data_build  # only needed to make them

# mm10 annotation: chrom sizes, blacklist, GENCODE vM25 GTF, TSS of protein-coding genes
mkdir -p $D/genome && cd $D/genome
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
wget -q -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
zcat mm10-blacklist.v2.bed.gz | cut -f1-3 > mm10-blacklist.v2.bed && rm mm10-blacklist.v2.bed.gz
wget -q -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
zcat gencode.vM25.annotation.gtf.gz | awk -F'\t' '$3=="gene" && $9 ~ /gene_type "protein_coding"/ {
    match($9, /gene_name "[^"]+"/); g = substr($9, RSTART+11, RLENGTH-12)
    tss = ($7 == "+") ? $4 - 1 : $5 - 1
    print $1"\t"tss"\t"tss+1"\t"g"\t.\t"$7 }' | sort -k1,1 -k2,2n > mm10_tss_protein_coding.bed

# GEO files
B=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218468/suppl
mkdir -p $D/raw && wget -q -c -P $D/raw $B/GSE218468_multiome_atac_counts.tsv.gz
mkdir -p $X/raw && cd $X/raw
for f in GSE218468_TEW_043783_523d9c_Multiome_Liver_10xprotocol_fragments.tsv.gz \
         GSE218468_TEW_ebb273_b33e6f_Multiome_Liver_CTRL_NSTprotocol_fragments.tsv.gz \
         GSE218468_multiome_rna_counts.tsv.gz; do
    wget -q -c $B/$f
done
ls -la $X/raw $D/raw
# GEO serves plain gzip; tabix needs BGZF. Recompress under short names and index.
F=$D/fragments; mkdir -p $F
zcat GSE218468_TEW_043783_523d9c_Multiome_Liver_10xprotocol_fragments.tsv.gz | bgzip -@ 4 > $F/Mouse_4_10x.fragments.tsv.gz
zcat GSE218468_TEW_ebb273_b33e6f_Multiome_Liver_CTRL_NSTprotocol_fragments.tsv.gz | bgzip -@ 4 > $F/Mouse_5_NST.fragments.tsv.gz
for f in $F/*.fragments.tsv.gz; do tabix -p bed $f; done
ls -la $F

# Mallet, for pycisTopic topic modelling
mkdir -p $X/software && cd $X/software
wget -q -c https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
[ -d Mallet-202108 ] || tar -xzf Mallet-202108-bin.tar.gz
echo DONE
