#!/bin/bash
#SBATCH --job-name=topic_models
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# pycisTopic topic models for part 2 (needs the outputs of part 1 in course_data/data_generated/).
#
# Run from the topic folder (single-cell-atac-seq), with the atac-sc environment active:
#     mkdir -p course_data_build/logs && sbatch scripts/12_topic_models.sh
set -euo pipefail
export JAVA_HOME=${JAVA_HOME:-$CONDA_PREFIX/lib/jvm}
python -I scripts/12_topic_models.py
