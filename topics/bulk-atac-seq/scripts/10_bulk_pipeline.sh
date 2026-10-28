#!/bin/bash
#SBATCH --job-name=bulk_atac
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-12%6
#SBATCH --output=data/logs/%x_%A_%a.log
# GSE246213 bulk ATAC, one array task per sample: the same steps part 1 teaches, genome-wide.
# Outputs go to data/data_generated/<sample>/; the chr19 subsets feed parts 1 and 2.
#
# Run from the topic folder (bulk-atac-seq), with the atac-cli environment active, after
# 00_download_refs.sh and 02_download_bulk.sh:
#     mkdir -p data/logs && sbatch scripts/10_bulk_pipeline.sh
# Each step is skipped if its output exists, so a failed task can simply be resubmitted.
set -euo pipefail
B=$PWD/data
G=$B/genome/mm10
T=${SLURM_CPUS_PER_TASK:-16}
read -r SAMPLE RUN <<< "$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i+1 {print $1, $2}' $B/samplesheet.tsv)"
O=$B/data_generated/$SAMPLE; mkdir -p $O; cd $O
TMP=$(mktemp -d -p ${TMPDIR:-/tmp}); trap 'rm -rf $TMP' EXIT
echo "== $SAMPLE ($RUN) on $(hostname), $(date)"

# 1. adapter/quality trimming
[ -s $SAMPLE.raw.bam ] || fastp -w $T --detect_adapter_for_pe --length_required 25 \
    -i $B/raw/${RUN}_1.fastq.gz -I $B/raw/${RUN}_2.fastq.gz \
    -o $TMP/R1.fq.gz -O $TMP/R2.fq.gz -j $SAMPLE.fastp.json -h $SAMPLE.fastp.html 2> $SAMPLE.fastp.log

# 2. alignment (fragments up to 2 kb); a read group is required by Picard
if [ ! -s $SAMPLE.raw.bam ]; then
    bowtie2 -p $T --very-sensitive -X 2000 --rg-id $SAMPLE --rg SM:$SAMPLE -x $G/bowtie2/mm10 \
        -1 $TMP/R1.fq.gz -2 $TMP/R2.fq.gz 2> $SAMPLE.bowtie2.log \
        | samtools sort -@ 4 -m 2G -T $TMP/sort -o $SAMPLE.raw.tmp.bam -
    mv $SAMPLE.raw.tmp.bam $SAMPLE.raw.bam
    samtools index $SAMPLE.raw.bam
fi

# 3. mark duplicates
picard -Xmx16g MarkDuplicates -I $SAMPLE.raw.bam -O $SAMPLE.markdup.bam -M $SAMPLE.markdup_metrics.txt \
    --READ_NAME_REGEX null --VALIDATION_STRINGENCY LENIENT --TMP_DIR $TMP --QUIET true
samtools index $SAMPLE.markdup.bam
samtools flagstat -@ $T $SAMPLE.markdup.bam > $SAMPLE.markdup.flagstat
samtools idxstats $SAMPLE.markdup.bam > $SAMPLE.markdup.idxstats

# 4. ATAC QC on the unfiltered (duplicate-marked) BAM
ataqv --threads $T --name $SAMPLE --metrics-file $SAMPLE.ataqv.json.gz \
    --tss-file $G/mm10_tss_protein_coding.bed --excluded-region-file $G/mm10-blacklist.v2.bed \
    --mitochondrial-reference-name chrM mouse $SAMPLE.markdup.bam > $SAMPLE.ataqv.txt

# 5. filter: proper pairs, MAPQ>=30, no unmapped/secondary/QC-fail/duplicate/supplementary, main chromosomes without chrM, no blacklist
CHROMS=$(grep -E '^chr([0-9]+|X|Y)\s' $G/mm10.chrom.sizes | cut -f1 | tr '\n' ' ')
samtools view -@ $T -b -f 2 -F 3852 -q 30 $SAMPLE.markdup.bam $CHROMS \
    | bedtools intersect -v -abam stdin -b $G/mm10-blacklist.v2.bed > $TMP/filt.bam
samtools sort -@ $T -m 2G -T $TMP/sort2 -o $SAMPLE.filtered.bam $TMP/filt.bam
samtools index $SAMPLE.filtered.bam
samtools flagstat -@ $T $SAMPLE.filtered.bam > $SAMPLE.filtered.flagstat
samtools view -b $SAMPLE.filtered.bam chr19 > $SAMPLE.filtered.chr19.bam
samtools index $SAMPLE.filtered.chr19.bam

# 6. peaks (paired-end mode) and coverage
macs3 callpeak -f BAMPE -t $SAMPLE.filtered.bam -n $SAMPLE -g mm --keep-dup all -q 0.01 --outdir macs3 2> $SAMPLE.macs3.log
bamCoverage -p $T -b $SAMPLE.filtered.bam -o $SAMPLE.bw --binSize 10 --normalizeUsing RPGC \
    --effectiveGenomeSize 2652783500 --extendReads 2> /dev/null

# 7. chr19 read pairs as untrimmed FASTQ (input of part 1): pairs with a primary alignment on chr19
samtools view -F 2304 $SAMPLE.raw.bam chr19 | cut -f1 | sort -u > $TMP/chr19.ids
for r in 1 2; do
    seqkit grep -j $T -f $TMP/chr19.ids $B/raw/${RUN}_${r}.fastq.gz -o ${SAMPLE}_chr19_R${r}.fastq.gz
done
echo "chr19 pairs: $(wc -l < $TMP/chr19.ids)"
rm -f $SAMPLE.raw.bam $SAMPLE.raw.bam.bai  # the duplicate-marked BAM contains all reads
echo "== DONE $SAMPLE $(date)"
