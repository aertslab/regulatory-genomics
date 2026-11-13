#!/bin/bash
#SBATCH --job-name=cbp_interpret
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Interpretation of the microglia model, the GPU part (TF-MoDISco on the result runs on the CPU, in 08_modisco.sh). The
# chrombpnet pipeline (06_chrombpnet.sh) already made the profile contribution scores of 30,000 random peaks and their
# TF-MoDISco motifs (evaluation/modisco_profile). This adds, in course_data_build/chrombpnet/interpret:
#     chrombpnet_nobias.counts_scores.h5 (+ .bw)  counts contribution scores of the same 30,000 peaks
#     footprints/                                 marginal footprints: the average predicted profile when a motif is
#                                                 inserted in the middle of the background regions of the test
#                                                 chromosomes, for the bias-free model (chrombpnet_nobias) and for the
#                                                 full model (chrombpnet, Tn5 bias included)
# The inserted motifs are the consensus sequences of the main TF-MoDISco patterns of the profile head.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet environment active, after 06_chrombpnet.sh, on an A100:
#     sbatch -A <account> -p gpu_a100_48C_96T_512GB scripts/chrombpnet/07_interpret.sh
set -euo pipefail
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
python -c "import tensorflow as tf, sys; g = tf.config.list_physical_devices('GPU'); print(g); sys.exit(not g)"
G=$PWD/course_data_build/genome/mm10
C=$PWD/course_data_build/chrombpnet
M=$C/microglia
O=$C/interpret; mkdir -p $O/footprints

chrombpnet contribs_bw -m $M/models/chrombpnet_nobias.h5 -r $M/auxiliary/30K_subsample_peaks.bed -g $G/mm10.fa \
    -c $G/mm10.main.chrom.sizes -op $O/chrombpnet_nobias -pc counts

cat > $O/footprints/motifs.tsv <<'EOF'
PU.1	AAGAGGAAGT
PU.1-IRF	GGAAGTGAAAC
IRF	GAAACTGAAACT
CEBP	GTTGCACAAT
MEF2	CTAAAAATAGC
MAF	GCTGAGTCAGCA
MITF-TFE	GTCACGTG
CTCF	CGCCCCCTGGTGG
SP-KLF	GCCCCGCCC
NFY	CCAATCA
NRF1	GCGCATGCGC
NFI	TGGCAGGCTGCCAG
CRE	TGACGTCA
EOF
for m in chrombpnet_nobias chrombpnet; do
    chrombpnet footprints -m $M/models/$m.h5 -r $C/nonpeaks_negatives.bed -g $G/mm10.fa -fl $C/fold_0.json \
        -op $O/footprints/$m -pwm_f $O/footprints/motifs.tsv
done
echo DONE
