#!/bin/bash
#SBATCH --job-name=bundle_cbp
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Packs what the ChromBPNet notebook downloads into course_data_build/chrombpnet_microglia.tar.gz (the notebook gets
# the mm10 sequence of chr9 and chr18 from UCSC itself):
#     peaks.narrowPeak, nonpeaks_negatives.bed, fold_0.json   the ChromBPNet inputs (04_prep_chrombpnet.sh)
#     fragments_chr18.tsv.gz (+ .tbi)  microglia fragments on chr18, to show the format
#     insertions_chr9_chr18.bw         Tn5 insertions per bp on the test chromosomes (made by the chrombpnet pipeline)
#     models/                          bias.h5 (05_bias.sh), and chrombpnet.h5 (with bias), chrombpnet_nobias.h5 and
#                                      bias_model_scaled.h5 (the bias model as used inside chrombpnet.h5) (06_chrombpnet.sh)
#     training/                        loss per epoch and model parameters of both models
#     evaluation/                      predictions of chrombpnet.h5 and bias_model_scaled.h5 on the test peaks, the metrics
#                                      of both pipelines, and the Tn5 motif found in the insertions (bw_shift_qc.png)
#     modisco/                         TF-MoDISco results of the bias model (profile) and of ChromBPNet (profile, counts),
#                                      with their TOMTOM matches (08_modisco.sh)
#     footprints/                      marginal footprints of both ChromBPNet models and the inserted motifs (07_interpret.sh)
# The notebook downloads it from BUNDLE_URL (scripts/build_chrombpnet.py), so upload it there.
#
# The inputs to rerun the training (Steps 4 and 5 of the notebook) go into a second, optional file,
# course_data_build/chrombpnet_microglia_training.tar (TRAINING_URL in scripts/build_chrombpnet.py): the full microglia
# fragment file (+ index) and the sizes of the main chromosomes. The complete mm10 sequence comes from UCSC.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet environment active, after 08_modisco.sh:
#     sbatch scripts/chrombpnet/09_bundle.sh
set -euo pipefail
B=$PWD/course_data_build
C=$B/chrombpnet
CB=$C/bias
CM=$C/microglia
CI=$C/interpret
TMP=$(mktemp -d -p ${TMPDIR:-/tmp}); trap 'rm -rf $TMP' EXIT
S=$TMP/bundle
mkdir -p $S/models $S/training $S/evaluation $S/modisco $S/footprints

cp $C/peaks.narrowPeak $C/nonpeaks_negatives.bed $C/fold_0.json $S/
tabix $B/fragments_mgl/microglia.fragments.tsv.gz chr18 | bgzip > $S/fragments_chr18.tsv.gz
tabix -p bed $S/fragments_chr18.tsv.gz

for chrom in chr9 chr18; do
    bigWigToBedGraph -chrom=$chrom $CM/auxiliary/data_unstranded.bw $TMP/$chrom.bedGraph
done
cat $TMP/chr9.bedGraph $TMP/chr18.bedGraph | LC_COLLATE=C sort -k1,1 -k2,2n > $TMP/test.bedGraph
awk '$1 == "chr9" || $1 == "chr18"' $B/genome/mm10/mm10.main.chrom.sizes > $TMP/test.chrom.sizes
bedGraphToBigWig $TMP/test.bedGraph $TMP/test.chrom.sizes $S/insertions_chr9_chr18.bw

cp $CB/models/bias.h5 $CM/models/chrombpnet.h5 $CM/models/chrombpnet_nobias.h5 $CM/models/bias_model_scaled.h5 $S/models/
cp $CB/logs/bias.log $CB/logs/bias_model_params.tsv $CM/logs/chrombpnet.log $CM/logs/chrombpnet_model_params.tsv $S/training/
cp $CM/evaluation/chrombpnet_predictions.h5 $CM/evaluation/bias_predictions.h5 $CM/evaluation/chrombpnet_metrics.json \
    $CB/evaluation/bias_metrics.json $CB/evaluation/bw_shift_qc.png $S/evaluation/
cp $CB/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 $S/modisco/bias_profile.h5
cp $CM/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 $S/modisco/chrombpnet_profile.h5
cp $CI/modisco_results_counts_scores.h5 $S/modisco/chrombpnet_counts.h5
for name in bias_profile chrombpnet_profile chrombpnet_counts; do cp $CI/tomtom/$name.tsv $S/modisco/${name}_tomtom.tsv; done
cp $CI/footprints/motifs.tsv $CI/footprints/chrombpnet_nobias_footprints.h5 $CI/footprints/chrombpnet_footprints.h5 $S/footprints/
du -sh $S/* | sort -h

tar -czf $B/chrombpnet_microglia.tar.gz -C $S .
# not compressed again: the fragment file is bgzipped already
tar -cf $B/chrombpnet_microglia_training.tar -C $B/fragments_mgl microglia.fragments.tsv.gz microglia.fragments.tsv.gz.tbi \
    -C $B/genome/mm10 mm10.main.chrom.sizes
ls -la $B/chrombpnet_microglia.tar.gz $B/chrombpnet_microglia_training.tar
echo DONE
