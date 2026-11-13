#!/bin/bash
#SBATCH --job-name=cbp_bias
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Tn5 bias model, trained on the microglia background regions (chrombpnet bias pipeline, -b 0.5 as recommended for
# ATAC). Output in course_data_build/chrombpnet/bias: models/bias.h5 and the QC report evaluation/overall_report.html.
# Check the report before training the ChromBPNet model: the bias model should learn the Tn5 motif and no TF motifs.
#
# ChromBPNet uses TensorFlow 2.8, which only supports GPUs up to the A100. Run from the topic folder (s2f-model-training),
# with the chrombpnet environment active, after 04_prep_chrombpnet.sh, on an A100 (on the VIB cluster):
#     sbatch -A <account> -p gpu_a100_48C_96T_512GB scripts/chrombpnet/05_bias.sh
set -euo pipefail
module load SAMtools/1.18-GCC-12.3.0 BEDTools/2.31.0-GCC-12.3.0
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}   # CUDA 11.2 and cuDNN 8.1 of the environment
export PYTHONUNBUFFERED=1   # progress in the log as it happens
python -c "import tensorflow as tf, sys; g = tf.config.list_physical_devices('GPU'); print(g); sys.exit(not g)"
G=$PWD/course_data_build/genome/mm10
C=$PWD/course_data_build/chrombpnet
O=$C/bias   # chrombpnet makes it, and refuses to run if it exists

chrombpnet bias pipeline -ifrag $PWD/course_data_build/fragments_mgl/microglia.fragments.tsv.gz -d ATAC \
    -g $G/mm10.fa -c $G/mm10.main.chrom.sizes -p $C/peaks.narrowPeak -n $C/nonpeaks_negatives.bed \
    -fl $C/fold_0.json -b 0.5 -o $O
echo DONE
