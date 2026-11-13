#!/bin/bash
#SBATCH --job-name=modisco_counts
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=course_data_build/logs/%x_%j.log
# TF-MoDISco motifs from the counts contribution scores of 07_interpret.sh, with the settings the chrombpnet pipeline
# uses for the profile scores (at most 50,000 seqlets, the central 500 bp of each peak) and the same motif database for
# the TOMTOM matches. The pipeline only runs TF-MoDISco on the profile head (where in a region the reads fall); the
# counts head predicts how many reads a region gets. Then the TOMTOM matches of the three TF-MoDISco results (bias
# model profile, ChromBPNet profile and counts) as tables, which the notebook loads (run_tomtom takes several minutes).
# Output in course_data_build/chrombpnet/interpret: modisco_results_counts_scores.h5, modisco_counts/motifs.html and
# tomtom/{bias_profile,chrombpnet_profile,chrombpnet_counts}.tsv. TF-MoDISco is skipped if its results exist.
#
# Run from the topic folder (s2f-model-training), with the chrombpnet environment active, after 07_interpret.sh:
#     sbatch scripts/chrombpnet/08_modisco.sh
set -euo pipefail
C=$PWD/course_data_build/chrombpnet
O=$C/interpret
MEME=$(python -c "import chrombpnet, os; print(os.path.dirname(chrombpnet.__file__))" 2> /dev/null)/data/motifs.meme.txt

[ -s $O/modisco_results_counts_scores.h5 ] || \
    modisco motifs -i $O/chrombpnet_nobias.counts_scores.h5 -n 50000 -w 500 -o $O/modisco_results_counts_scores.h5
[ -s $O/modisco_counts/motifs.html ] || \
    modisco report -i $O/modisco_results_counts_scores.h5 -o $O/modisco_counts -m $MEME

mkdir -p $O/tomtom
python - $MEME $O/tomtom $C/bias/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 \
    $C/microglia/auxiliary/interpret_subsample/modisco_results_profile_scores.h5 $O/modisco_results_counts_scores.h5 <<'EOF'
import sys
from modiscolite.report import run_tomtom
meme, out = sys.argv[1:3]
for name, h5 in zip(["bias_profile", "chrombpnet_profile", "chrombpnet_counts"], sys.argv[3:]):
    run_tomtom(h5, None, meme, top_n_matches=3).to_csv(f"{out}/{name}.tsv", sep="\t", index=False)
    print(name, "done")
EOF
echo DONE
