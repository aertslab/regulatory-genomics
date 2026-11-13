#!/bin/bash
#SBATCH --job-name=refs_mm10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# mm10 references, cell metadata and 10x Multiome barcode lists for the ChromBPNet notebook (Zemke et al. 2023, mouse
# MOp 10x Multiome, GEO GSE229169), in course_data_build: genome, chrom sizes, blacklist and chromap index in
# genome/mm10, the cell metadata, barcode lists and run sheet in meta.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet-prep environment active:
#     mkdir -p course_data_build/logs && sbatch scripts/chrombpnet/00_download_refs.sh
set -euo pipefail
module load SAMtools/1.18-GCC-12.3.0
B=$PWD/course_data_build
G=$B/genome/mm10
mkdir -p $G $B/meta

cd $B/genome/mm10
wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
[ -s $G/mm10.fa ] || zcat mm10.fa.gz > $G/mm10.fa
samtools faidx $G/mm10.fa
wget -q -c -P $G https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
# ChromBPNet (folds, non-peaks, training) only uses the main chromosomes
grep -P '^chr([0-9]+|X|Y)\t' $G/mm10.chrom.sizes | sort -V > $G/mm10.main.chrom.sizes
wget -q -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
zcat mm10-blacklist.v2.bed.gz | cut -f1-3 > $G/mm10-blacklist.v2.bed
[ -s chromap.index ] || chromap -i -t ${SLURM_CPUS_PER_TASK:-8} -r $G/mm10.fa -o chromap.index

cd $B/meta
# cell barcodes (GEX barcodes, as reported by cellranger-arc) with their sample and cell type
wget -q -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229169/suppl/GSE229169_mouse_filtered_metadata.tsv.gz
# 10x Multiome barcode lists: line n of the ATAC list and line n of the GEX list are on the same gel bead
wget -q -c https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/atac_737K-arc-v1.txt.gz
wget -q -c https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/gex_737K-arc-v1.txt.gz

# run sheet: the SRA runs of the 8 mouse ATAC libraries (subregions 2C, 3C, 4B, 5D, 2 replicates each)
GSM="GSM7156248 GSM7156250 GSM7156252 GSM7156254 GSM7156256 GSM7156258 GSM7156260 GSM7156262"
Q=$(for g in $GSM; do printf 'sample_alias="%s" OR ' $g; done | sed 's/ OR $//')
{ printf 'library\tgsm\trun\tspots\n'
  curl -sG https://www.ebi.ac.uk/ena/portal/api/search --data-urlencode result=read_run --data-urlencode "query=$Q" \
      --data-urlencode fields=run_accession,sample_alias,sample_title,read_count --data-urlencode format=tsv \
      | awk -F'\t' 'NR > 1 {sub(/ \(Multiome ATAC\)/, "", $3); print $3"\t"$2"\t"$1"\t"$4}' | sort; } > runs.tsv
echo "$(($(wc -l < runs.tsv) - 1)) runs in runs.tsv"
ls -la $G $B/meta
echo DONE
