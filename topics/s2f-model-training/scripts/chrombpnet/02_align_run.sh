#!/bin/bash
#SBATCH --job-name=align_atac
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-31%6
#SBATCH --output=course_data_build/logs/%x_%A_%a.log
# One array task per SRA run (row of course_data_build/meta/runs.tsv): download the run from NCBI (the FASTQs of ENA
# lack the barcode read), align it with chromap in chunks of 50M spots, and keep only the fragments of the microglia
# of that library. The download and the FASTQs are deleted as soon as they are used. Libraries with fewer than 50
# microglia (Mop_3C_rep1 and Mop_4B_rep1, 9 and 5) are skipped.
#
# Read layout (see 01_barcodes.sh): reads 1 and 3 are the genomic reads, read 2 has the cell barcode at bases 9-24,
# reverse-complemented relative to the 10x ATAC barcode list. chromap corrects the barcode against the reverse-
# complemented list, removes PCR duplicates per cell, and writes the paired GEX barcode (the one of the metadata).
# Output: course_data_build/fragments_mgl/runs/<run>.tsv.gz (chrom, start, end, GEX barcode, number of reads).
#
# Run from the topic folder (s2f-model-training), with the chrombpnet-prep environment active, after 00 and 01:
#     sbatch scripts/chrombpnet/02_align_run.sh
# A run whose output exists is skipped, so failed tasks can simply be resubmitted. The two runs of Mop_5D_rep1 (tasks 26
# and 27, 1-1.2 billion spots) take longest; start them first, with more CPUs:
#     sbatch --cpus-per-task=16 --mem=96G --array=26,27 scripts/chrombpnet/02_align_run.sh
#     sbatch --array=2-25,28-31%6 scripts/chrombpnet/02_align_run.sh
set -euo pipefail
B=$PWD/course_data_build
G=$B/genome/mm10
M=$B/meta
T=${SLURM_CPUS_PER_TASK:-8}
CHUNK=50000000   # chromap keeps the mappings of a chunk in memory (~64 GB for 89M spots)
read -r LIB RUN SPOTS <<< "$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR == i + 1 {print $1, $3, $4}' $M/runs.tsv)"
O=$B/fragments_mgl/runs; mkdir -p $O
[ -s $O/$RUN.tsv.gz ] && { echo "$RUN: done"; exit 0; }
[ $(wc -l < $M/mgl_barcodes/$LIB.txt) -lt 50 ] && { echo "$RUN: skipped, few microglia in $LIB"; exit 0; }
TMP=$(mktemp -d -p ${TMPDIR:-/tmp}); trap 'rm -rf $TMP' EXIT
echo "== $LIB $RUN ($SPOTS spots) on $(hostname), $(date)"

cd $TMP
prefetch -q --max-size u $RUN
for c in $(seq 0 $(( (SPOTS - 1) / CHUNK ))); do
    from=$(( c * CHUNK + 1 )); to=$(( (c + 1) * CHUNK < SPOTS ? (c + 1) * CHUNK : SPOTS ))
    echo "-- spots $from-$to, $(date)"
    fastq-dump -N $from -X $to --split-files -O fq $RUN > /dev/null
    chromap --preset atac -t $T -x $B/genome/mm10/chromap.index -r $G/mm10.fa \
        -1 fq/${RUN}_1.fastq -2 fq/${RUN}_3.fastq -b fq/${RUN}_2.fastq --read-format r1:0:-1,r2:0:-1,bc:8:23 \
        --barcode-whitelist $M/atac_rc_737K-arc-v1.txt --barcode-translate $M/translate_atac_rc_to_gex.tsv \
        -o frag.tsv 2> $O/$RUN.chunk$c.chromap.log
    awk -F'\t' 'NR == FNR {mgl[$1]; next} ($4 in mgl)' $M/mgl_barcodes/$LIB.txt frag.tsv >> mgl.tsv
    echo "   $(wc -l < frag.tsv) fragments, $(wc -l < mgl.tsv) microglia fragments so far"
    rm -rf fq frag.tsv
done
rm -rf $RUN   # the download
sort -k1,1 -k2,2n -k3,3n -k4,4 mgl.tsv | bgzip > mgl.tsv.gz
mv mgl.tsv.gz $O/$RUN.tsv.gz
echo "DONE $(date)"
