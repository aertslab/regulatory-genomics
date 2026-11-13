#!/bin/bash
#SBATCH --job-name=cbp_mgl
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# Bias-factorized ChromBPNet model of the microglia (chrombpnet pipeline), with the bias model of 05_bias.sh. Output in
# course_data_build/chrombpnet/microglia: models/chrombpnet.h5 (with bias), models/chrombpnet_nobias.h5 (bias removed,
# the model to interpret), models/bias_model_scaled.h5, and the QC report evaluation/overall_report.html.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet environment active, after 05_bias.sh, on an A100:
#     sbatch -A <account> -p gpu_a100_48C_96T_512GB scripts/chrombpnet/06_chrombpnet.sh
set -euo pipefail
module load SAMtools/1.18-GCC-12.3.0 BEDTools/2.31.0-GCC-12.3.0
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
python -c "import tensorflow as tf, sys; g = tf.config.list_physical_devices('GPU'); print(g); sys.exit(not g)"
G=$PWD/course_data_build/genome/mm10
C=$PWD/course_data_build/chrombpnet
O=$C/microglia   # chrombpnet makes it, and refuses to run if it exists

chrombpnet pipeline -ifrag $PWD/course_data_build/fragments_mgl/microglia.fragments.tsv.gz -d ATAC \
    -g $G/mm10.fa -c $G/mm10.main.chrom.sizes -p $C/peaks.narrowPeak -n $C/nonpeaks_negatives.bed \
    -fl $C/fold_0.json -b $C/bias/models/bias.h5 -o $O
echo DONE
