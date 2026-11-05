"""Pre-compute the pycisTopic topic models that part 2 loads.

Uses the cells and consensus peaks written by part 1 (copied to course_data/data_generated/).
Submit with scripts/12_topic_models.sh; writes course_data/data_generated/topic_models.pkl.
"""
import os
import pickle
from pathlib import Path

import pandas as pd
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments, merge
from pycisTopic.lda_models import run_cgs_models_mallet

SC = str(Path(__file__).resolve().parents[1] / "course_data")
GEN = f"{SC}/data_generated"
MALLET = str(Path(__file__).resolve().parents[1] / "course_data_build" / "software" / "Mallet-202108" / "bin" / "mallet")  # from 01_download_sc.sh
N_CPU = int(os.environ.get("SLURM_CPUS_PER_TASK", 16))
TMP = os.path.join(os.environ.get("TMPDIR", "/tmp"), "mallet")
os.makedirs(TMP, exist_ok=True)

cells = pd.read_csv(f"{GEN}/cells.tsv", sep="\t", index_col=0)
objects = []
for sample in ["Mouse_4_10x", "Mouse_5_NST"]:
    objects.append(create_cistopic_object_from_fragments(
        path_to_fragments=f"{SC}/fragments/{sample}.fragments.tsv.gz",
        path_to_regions=f"{GEN}/consensus_peaks.bed",
        path_to_blacklist=f"{SC}/genome/mm10-blacklist.v2.bed",
        valid_bc=cells.index[cells["sample"] == sample].tolist(),
        n_cpu=N_CPU, project=sample))
cistopic_obj = merge(objects)
print(cistopic_obj, flush=True)

os.environ["MALLET_MEMORY"] = "100G"
models = run_cgs_models_mallet(
    cistopic_obj, n_topics=[10, 20, 30, 40, 50], n_cpu=N_CPU, n_iter=500, random_state=555,
    alpha=50, alpha_by_topic=True, eta=0.1, eta_by_topic=False,
    tmp_path=TMP, save_path=None, mallet_path=MALLET)
with open(f"{GEN}/topic_models.pkl", "wb") as f:
    pickle.dump(models, f)
print("saved", [m.n_topic for m in models])
