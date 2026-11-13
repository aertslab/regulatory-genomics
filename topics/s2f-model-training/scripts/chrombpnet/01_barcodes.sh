#!/bin/bash
# Cell barcodes: checks the read layout of the SRA runs on the first 100,000 spots of one run, then writes the barcode
# files that 02_align_run.sh uses, in course_data_build/meta:
#     atac_rc_737K-arc-v1.txt        the 10x ATAC barcode list, reverse-complemented (as the barcode read is)
#     translate_atac_rc_to_gex.tsv   GEX barcode <tab> ATAC barcode (rc): chromap writes the GEX barcode, as in the metadata
#     mgl_barcodes/<library>.txt     GEX barcodes of the microglia (subclass MGC) of each library
#
# Result of the check (SRR24107412): reads 1 and 3 are the genomic reads (101 bp), read 2 is the 24 bp i5 read (8 bp
# spacer + 16 bp cell barcode). 97% of the barcodes match the reverse-complemented ATAC list, none the list as is. All
# 47,404 cell barcodes of the metadata are in the GEX list; 3,129 cells are microglia.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet-prep environment active, after 00_download_refs.sh:
#     bash scripts/chrombpnet/01_barcodes.sh [SRR24107412]
set -euo pipefail
M=$PWD/course_data_build/meta
RUN=${1:-SRR24107412}
D=$PWD/course_data_build/check_barcodes; mkdir -p $D && cd $D
fastq-dump -X 100000 --split-files -O . $RUN > /dev/null   # keeps the technical (barcode) read by default

echo "== read lengths per file"
for f in ${RUN}_*.fastq; do echo "$f $(awk 'NR % 4 == 2 {print length($0)}' $f | sort | uniq -c | sort -k1,1nr | head -1)"; done

echo "== fraction of barcodes (last 16 bp of the 24 bp read) in the ATAC list, as is and reverse-complemented"
BC=$(for f in ${RUN}_*.fastq; do if [ "$(awk 'NR == 2 {print length($0)}' $f)" = 24 ]; then echo $f; fi; done)
awk 'NR % 4 == 2 {print substr($0, 9, 16)}' $BC > bc.txt
zcat $M/atac_737K-arc-v1.txt.gz > wl.txt
rev wl.txt | tr ACGT TGCA > wl_rc.txt
for w in wl.txt wl_rc.txt; do
    awk -v w=$w 'NR == FNR {ok[$1]; next} {n++; if ($1 in ok) m++} END {printf "%s: %.3f\n", w, m / n}' $w bc.txt
done

echo "== GEO cell metadata"
zcat $M/GSE229169_mouse_filtered_metadata.tsv.gz | awk 'NR <= 3'
echo "-- cells per library (orig.ident)"
zcat $M/GSE229169_mouse_filtered_metadata.tsv.gz | awk 'NR > 1 {print $2}' | sort | uniq -c
echo "-- cells per subclass"
zcat $M/GSE229169_mouse_filtered_metadata.tsv.gz | awk 'NR > 1 {print $3}' | sort | uniq -c | sort -k1,1nr
echo "-- fraction of metadata barcodes in the GEX list (should be 1)"
zcat $M/gex_737K-arc-v1.txt.gz > gex.txt
zcat $M/GSE229169_mouse_filtered_metadata.tsv.gz | awk 'NR > 1 {sub(/^.*_/, "", $1); sub(/-1$/, "", $1); print $1}' \
    | awk 'NR == FNR {ok[$1]; next} {n++; if ($1 in ok) m++} END {printf "%.4f of %d\n", m / n, n}' gex.txt -
rm -f ${RUN}_*.fastq bc.txt wl.txt wl_rc.txt gex.txt

echo "== barcode files"
cd $M
# line n of the GEX list and line n of the ATAC list are on the same gel bead
zcat atac_737K-arc-v1.txt.gz | rev | tr ACGT TGCA > atac_rc_737K-arc-v1.txt
paste <(zcat gex_737K-arc-v1.txt.gz) atac_rc_737K-arc-v1.txt > translate_atac_rc_to_gex.tsv
mkdir -p mgl_barcodes
for LIB in $(awk 'NR > 1 {print $1}' runs.tsv | sort -u); do
    ID=$(echo $LIB | tr A-Z a-z | sed 's/_rep//; s/_//')   # Mop_2C_rep1 -> mop2c1, as in orig.ident
    zcat GSE229169_mouse_filtered_metadata.tsv.gz | awk -v id=$ID '$2 == id && $3 == "MGC" {
        sub(/^.*_/, "", $1); sub(/-1$/, "", $1); print $1}' > mgl_barcodes/$LIB.txt
done
wc -l mgl_barcodes/*.txt
