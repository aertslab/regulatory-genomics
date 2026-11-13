#!/bin/bash
#SBATCH --job-name=pseudobulk_mgl
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Microglia pseudobulk: merges the microglia fragments of all runs (02_align_run.sh) into one fragment file, the input of
# ChromBPNet, and calls microglia peaks on it with MACS3. Most libraries were sequenced in several runs, so duplicates
# are removed again (same position and same cell = PCR duplicate). The barcodes get the library as prefix, as in the
# metadata (mop2c1_<barcode>), because the same barcode can occur in different libraries. Only the main chromosomes
# (chr1-19, X, Y) are kept: chrM and the unplaced contigs hold 0.25% of the fragments, and chrombpnet fails on
# chromosomes missing from its chrom sizes file. The per-run fragments are deleted at the end.
# Output in course_data_build/fragments_mgl: microglia.fragments.tsv.gz (+ .tbi) and macs3/microglia_peaks.narrowPeak.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet-prep environment active, after all tasks of
# 02_align_run.sh:
#     sbatch scripts/chrombpnet/03_pseudobulk.sh
set -euo pipefail
G=$PWD/course_data_build/genome/mm10
M=$PWD/course_data_build/meta
F=$PWD/course_data_build/fragments_mgl
T=${SLURM_CPUS_PER_TASK:-8}
TMP=$(mktemp -d -p ${TMPDIR:-/tmp}); trap 'rm -rf $TMP' EXIT

# the runs of the libraries with microglia (02_align_run.sh skips libraries with fewer than 50)
awk 'NR > 1 {print $1, $3}' $M/runs.tsv | while read -r LIB RUN; do
    [ $(wc -l < $M/mgl_barcodes/$LIB.txt) -lt 50 ] && continue
    [ -s $F/runs/$RUN.tsv.gz ] || { echo "missing $F/runs/$RUN.tsv.gz" >&2; exit 1; }
    ID=$(echo $LIB | tr A-Z a-z | sed 's/_rep//; s/_//')
    zcat $F/runs/$RUN.tsv.gz | awk -v id=$ID 'BEGIN {OFS = "\t"} NR == FNR {ok[$1]; next} ($1 in ok) {$4 = id "_" $4; print}' \
        $G/mm10.main.chrom.sizes -
done > $TMP/all.tsv
echo "$(wc -l < $TMP/all.tsv) fragments before removing duplicates between runs"

# one line per fragment and cell; column 5 sums the reads of the duplicates
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 16G --parallel=$T -T $TMP $TMP/all.tsv \
    | awk 'BEGIN {OFS = "\t"} {k = $1 FS $2 FS $3 FS $4}
           k == p {n += $5; next} NR > 1 {print c, s, e, b, n} {p = k; c = $1; s = $2; e = $3; b = $4; n = $5}
           END {print c, s, e, b, n}' \
    | bgzip -@ $T > $F/microglia.fragments.tsv.gz
tabix -f -p bed $F/microglia.fragments.tsv.gz
rm $TMP/all.tsv
zcat $F/microglia.fragments.tsv.gz | awk '{n++; c[$4]++} END {printf "%d fragments, %d cells, %.0f fragments per cell\n", n, length(c), n / length(c)}'

# peaks on the Tn5 insertion sites (both ends of each fragment), as in the ChromBPNet tutorial
zcat $F/microglia.fragments.tsv.gz | awk 'BEGIN {OFS = "\t"} {print $1, $2, $2 + 1; print $1, $3 - 1, $3}' > $TMP/ins.bed
macs3 callpeak -t $TMP/ins.bed -f BED -g mm -n microglia --outdir $F/macs3 \
    -p 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits 2> $F/macs3.log
awk '$9 >= 2' $F/macs3/microglia_peaks.narrowPeak | wc -l | xargs echo "peaks (summits) with q <= 0.01:"

rm $F/runs/*.tsv.gz
echo DONE
