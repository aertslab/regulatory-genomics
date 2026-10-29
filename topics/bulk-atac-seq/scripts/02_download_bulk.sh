#!/bin/bash
#SBATCH --job-name=dl_bulk
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=course_data_build/logs/%x_%A_%a.log
# GSE182188 FASTQs (FOXA1-HNF4A line) from ENA into course_data_build/raw, md5-checked. One array task per sample
# (row of course_data_build/samplesheet.tsv, the sample sheet written in part 1).
#
# Run from the topic folder (bulk-atac-seq):
#     mkdir -p course_data_build/logs && sbatch scripts/02_download_bulk.sh
# or, without Slurm, one sample at a time:
#     for i in $(seq 1 4); do SLURM_ARRAY_TASK_ID=$i bash scripts/02_download_bulk.sh; done
set -euo pipefail
B=$PWD/course_data_build
RUN=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i+1 {print $2}' $B/samplesheet.tsv)
mkdir -p $B/raw && cd $B/raw
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n 1 > $RUN.ena.tsv
IFS=';' read -ra URLS <<< "$(cut -f2 $RUN.ena.tsv)"   # col 1 is run_accession
IFS=';' read -ra MD5S <<< "$(cut -f3 $RUN.ena.tsv)"
for i in "${!URLS[@]}"; do
    f=$(basename ${URLS[$i]})
    ok=0
    for attempt in 1 2 3; do
        if wget -q -c https://${URLS[$i]} -O $f && echo "${MD5S[$i]}  $f" | md5sum -c -; then ok=1; break; fi
        echo "retry $attempt for $f"; rm -f $f
    done
    [ $ok -eq 1 ] || { echo "FAILED $f"; exit 1; }
done
echo DONE $RUN
