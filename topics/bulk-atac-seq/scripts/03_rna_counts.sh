#!/bin/bash
# RNA-seq gene counts of the FOXA1-HNF4A line (GSE182190, used in part 3): the Salmon quantifications of GEO,
# summed per gene symbol (field 6 of the GENCODE transcript names) and rounded. Writes course_data/rna/foxa1hnf4a_gene_counts.tsv;
# the downloads stay in course_data_build/rna.
#
# Run from the topic folder (bulk-atac-seq):
#     bash scripts/03_rna_counts.sh
set -euo pipefail
R=$PWD/course_data_build/rna
OUT=$PWD/course_data/rna
mkdir -p $R/ex $OUT && cd $R
SAMPLES="GSM5520301_RNAseq_FoxA1-Hnf4a_uninduced_rep1:U1 GSM5520302_RNAseq_FoxA1-Hnf4a_uninduced_rep2:U2
GSM5520303_RNAseq_FoxA1-Hnf4a_uninduced_rep3:U3 GSM5520304_RNAseq_FoxA1-Hnf4a_induced_rep1:I1
GSM5520305_RNAseq_FoxA1-Hnf4a_induced_rep2:I2 GSM5520306_RNAseq_FoxA1-Hnf4a_induced_rep3:I3"
for s in $SAMPLES; do
    f=${s%%:*}_quants.tar.gz; gsm=${s:0:10}
    [ -s $f ] || curl -sS -o $f "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/$gsm/suppl/$f"
    mkdir -p ex/${s%%:*} && tar -xzf $f -C ex/${s%%:*}
done
# one column per sample: NumReads summed per gene symbol
for s in $SAMPLES; do
    tail -n +2 ex/${s%%:*}/*/quant.sf | awk -F'\t' '{split($1, n, "|"); c[n[6]] += $5} END {for (g in c) print g "\t" c[g]}' \
        | sort -k1,1 > ex/${s##*:}.tsv
done
( printf "gene"; for s in $SAMPLES; do printf "\tFOXA1HNF4A_%s" ${s##*:}; done; echo
  join -t $'\t' ex/U1.tsv ex/U2.tsv | join -t $'\t' - ex/U3.tsv | join -t $'\t' - ex/I1.tsv | join -t $'\t' - ex/I2.tsv \
      | join -t $'\t' - ex/I3.tsv | awk 'BEGIN {OFS = "\t"} {for (i = 2; i <= NF; i++) $i = sprintf("%.0f", $i); print}'
) > $OUT/foxa1hnf4a_gene_counts.tsv
wc -l $OUT/foxa1hnf4a_gene_counts.tsv
