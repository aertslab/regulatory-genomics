#!/bin/bash
#SBATCH --job-name=prep_cbp
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# ChromBPNet inputs for the microglia model, in course_data_build/chrombpnet:
#     peaks.narrowPeak      MACS3 summits with q <= 0.01 on the main chromosomes, at least 1,057 bp (half the 2,114 bp
#                           input of the model) away from blacklisted regions
#     fold_0.json           test chromosomes chr9 and chr18, validation chr8 and chr10 (as in the CREsted notebook)
#     nonpeaks_negatives.bed  background regions with the GC content of the peaks (intermediate files in
#                           nonpeaks_auxiliary)
#
# Run from the topic folder (s2f-model-training), with the chrombpnet environment active, after 03_pseudobulk.sh:
#     sbatch scripts/chrombpnet/04_prep_chrombpnet.sh
set -euo pipefail
module load BEDTools/2.31.0-GCC-12.3.0
G=$PWD/course_data_build/genome/mm10
F=$PWD/course_data_build/fragments_mgl
C=$PWD/course_data_build/chrombpnet; mkdir -p $C

bedtools slop -i $G/mm10-blacklist.v2.bed -g $G/mm10.main.chrom.sizes -b 1057 > $C/blacklist_slop1057.bed
awk 'NR == FNR {ok[$1]; next} ($1 in ok) && $9 >= 2' $G/mm10.main.chrom.sizes $F/macs3/microglia_peaks.narrowPeak \
    | bedtools intersect -v -a - -b $C/blacklist_slop1057.bed > $C/peaks.narrowPeak
rm $C/blacklist_slop1057.bed
echo "$(wc -l < $C/peaks.narrowPeak) peaks"

chrombpnet prep splits -c $G/mm10.main.chrom.sizes -tcr chr9 chr18 -vcr chr8 chr10 -op $C/fold_0
chrombpnet prep nonpeaks -g $G/mm10.fa -p $C/peaks.narrowPeak -c $G/mm10.main.chrom.sizes -fl $C/fold_0.json \
    -br $G/mm10-blacklist.v2.bed -o $C/nonpeaks
echo "$(wc -l < $C/nonpeaks_negatives.bed) non-peaks"
echo DONE
