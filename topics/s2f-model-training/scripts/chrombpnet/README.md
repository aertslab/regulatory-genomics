# s2f model training: making the ChromBPNet data

Notes for teachers. The ChromBPNet notebook downloads chromosomes 9 and 18 of mm10 from UCSC, and one file, `chrombpnet_microglia.tar.gz` (336 MB), from `BUNDLE_URL` in `../build_chrombpnet.py`. The scripts here make that file from the raw data of Zemke et al. (2023): microglia fragments from the SRA reads, a Tn5 bias model, the ChromBPNet model and its interpretation. Students don't need them.

The notebook itself is built with `../build_chrombpnet.py` (see the README at the root of the repository), with the `chrombpnet` kernel.

## Environments

- `chrombpnet-prep` (`../../envs/chrombpnet-prep.yml`): sra-tools, chromap, MACS3, htslib. Used by `00` to `03`.
- `chrombpnet` (`../../envs/chrombpnet.yml`): ChromBPNet 1.0.1 with TensorFlow 2.8. Used by `04` to `09` and by the notebook.

```bash
conda env create --solver=libmamba -f envs/chrombpnet-prep.yml
conda env create --solver=libmamba -f envs/chrombpnet.yml
```

samtools and bedtools come from the cluster modules (`module load SAMtools/1.18-GCC-12.3.0 BEDTools/2.31.0-GCC-12.3.0`; the scripts load them). TensorFlow 2.8 supports GPUs up to the A100; on newer GPUs it runs on the CPU without an error, so the GPU scripts check that a GPU is found and stop otherwise.

## Steps

Run from the topic folder (`s2f-model-training`), with the environment of the step active, e.g. `sbatch scripts/chrombpnet/05_bias.sh`; each script has its exact command at the top. Everything is written to `course_data_build/` (git-ignored); logs go to `course_data_build/logs/`. On the VIB cluster, pass the account to `sbatch` (`-A <account>`). The times below are those of the run that made the current bundle.

| Script | What | Resources | Time |
|---|---|---|---|
| `00_download_refs.sh` | mm10, blacklist, chromap index, cell metadata, 10x barcode lists, the list of SRA runs | 8 CPUs | 5 min |
| `01_barcodes.sh` | checks where the barcode sits in the reads; barcode translation table; microglia barcodes per library | login node | a few min |
| `02_align_run.sh` | per SRA run: download, align with chromap, keep the microglia fragments (array of 31 runs; libraries with fewer than 50 microglia are skipped) | 8 to 16 CPUs, 64 to 96 GB per run | 5 h with 6 runs at a time |
| `03_pseudobulk.sh` | merges the runs, removes duplicates, keeps the main chromosomes, calls peaks with MACS3 | 8 CPUs | 15 min |
| `04_prep_chrombpnet.sh` | peaks, folds, GC-matched background regions | 4 CPUs | 5 min |
| `05_bias.sh` | Tn5 bias model and its report | 1 A100 | 4.5 h, of which 15 min training |
| `06_chrombpnet.sh` | ChromBPNet and its report | 1 A100 | 4.5 h, of which 80 min training |
| `07_interpret.sh` | counts contribution scores, marginal footprints | 1 A100 | 1.2 h |
| `08_modisco.sh` | TF-MoDISco on the counts contribution scores | 8 CPUs | 1.2 h |
| `09_bundle.sh` | packs `course_data_build/chrombpnet_microglia.tar.gz` and `chrombpnet_microglia_training.tar` | 4 CPUs | 1 min |

Worth checking before `06`: the bias report, `course_data_build/chrombpnet/bias/evaluation/overall_report.html`. The TF-MoDISco motifs of the bias model should all be Tn5 motifs. After `06`, the report in `course_data_build/chrombpnet/microglia/evaluation/` shows the test performance and whether the bias correction is complete (the footprints of the Tn5 motifs should stay below 0.003).

## The bundle

Upload `course_data_build/chrombpnet_microglia.tar.gz` to `BUNDLE_URL`. The notebook downloads it into `data/chrombpnet/` next to the notebook and unpacks it there. If the file is already in `data/chrombpnet/`, the notebook uses it, so on a server without internet access it can be copied there instead.

Upload `course_data_build/chrombpnet_microglia_training.tar` (505 MB) to `TRAINING_URL`. It holds the full microglia fragment file (43.4 million fragments, 3,115 nuclei, main chromosomes) and the main chromosome sizes; the notebook only refers to it, in Step 4, for students who want to rerun the training. From this fragment file, the peaks (the MACS3 step of `03`) and everything after them (`04` to `09`, about 10 GPU hours) can be made again without the SRA reads.

After the upload, `course_data_build/` (19 GB, most of it the genome and chromap index) can be deleted. Making it again from the SRA reads takes about 20 hours.
