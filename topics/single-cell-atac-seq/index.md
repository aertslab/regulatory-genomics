# Single-cell ATAC-seq

In this practical we analyse the chromatin accessibility of single nuclei from the mouse liver:

> Bravo González-Blas C, et al. *Single-cell spatial multi-omics and deep learning dissect enhancer-driven gene regulatory networks in liver zonation.* Nature Cell Biology (2024).
> https://doi.org/10.1038/s41556-023-01316-4

We use the ATAC part of two 10x Multiome samples (GEO [GSE218468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218468)): Mouse_4, prepared with the 10x nuclei isolation protocol, and Mouse_5, prepared with the NST protocol, about 11,000 cells in total. Multiome measures ATAC and RNA in the same nucleus, so the cell types and hepatocyte zones annotated from the RNA in the paper give us reliable labels to check the ATAC analysis. The [single-nucleus RNA-seq practical](../single-nucleus-rna-seq/index.md) analyses single-nucleus RNA-seq data of the same study.

**Research questions:**

1. Can the cell types of the liver be recovered from accessibility alone, without RNA?
2. Is hepatocyte zonation (periportal versus pericentral) visible in the chromatin?
3. Which regions differ between pericentral and periportal hepatocytes, and are they near known zonated genes?

The practical consists of four notebooks:

1. [Part 1: quality control and preprocessing](part_1_qc_preprocessing.ipynb): fragment files, per-cell QC (fragments, TSS enrichment, fragment sizes), doublets, cell-type and zone annotations checked with gene activity, pseudobulk peak calling and a consensus peak set (SnapATAC2, MACS3).
2. [Part 2: topic modelling with pycisTopic](part_2_pycistopic.ipynb): topic models, model selection, batch correction (Harmony), clustering and UMAP, cell-type and zonation topics.
3. [Part 3: dimensionality reduction and clustering with SnapATAC2](part_3_snapatac2.ipynb): spectral embedding, batch correction, clustering, and a comparison with pycisTopic.
4. [Part 4: differentially accessible regions](part_4_dars.ipynb): cell-type DARs, and pericentral versus periportal DARs with SnapATAC2 (per mouse) and pycisTopic (imputed accessibility), compared with each other and with known zonated genes, pseudobulk coverage tracks per cell type and zone (pycisTopic) to look at them in IGV; a comparison with the paper and the conclusion.

## Software

All four notebooks use one mamba environment, `atac-sc` (`envs/atac-sc.yml`), with SnapATAC2, pycisTopic, scanpy and MACS3; the kernel is "Python (atac-sc)". The commands to create it and register the kernel are given at the start of part 1. Part 2 also uses [Mallet](https://mimno.github.io/Mallet/) for topic modelling, stored in `data/software`; it needs Java, which is included in the environment.

## Data

All data is stored in the folder `data`, next to the notebooks; the notebooks write their results to `data/results/part_N`. The notebooks show the download commands as text (so they are not run by accident); the data are built by the scripts in `scripts/`, run from this folder (each script has a Slurm header, and the usage in its first lines):

| Script | Does | Time |
|---|---|---|
| `01_download_sc.sh` | fragment files (recompressed with bgzip and indexed), the published count matrices, mm10 annotation files and Mallet | ~30 min |
| `13_annotations.py` | `data/annotations/cell_annotations.tsv` from the public Scope loom of the paper (download instructions in the script) | minutes |
| part 1 | consensus peaks, cell list and peak matrix; copy `liver_peaks.h5ad`, `consensus_peaks.bed` and `cells.tsv` to `data/data_generated` | ~45 min |
| `12_topic_models.sh` | pycisTopic models with 10 to 50 topics (`data/data_generated/topic_models.pkl`) | ~3 h (16 CPU) |

`data/data_generated` also holds the results of each notebook that the next one needs (imported samples, pseudobulk peaks, cisTopic object, SnapATAC2 matrix), so every notebook can also be run on its own. The notebooks need about 12 GB of disk space, plus about 5 GB for the pseudobulk tracks of part 4.
