---
downloads:
  - file: scripts/01_download_sc.sh
    title: 01_download_sc.sh
    filename: 01_download_sc.sh
  - file: scripts/13_annotations.py
    title: 13_annotations.py
    filename: 13_annotations.py
  - file: scripts/12_topic_models.sh
    title: 12_topic_models.sh
    filename: 12_topic_models.sh
  - file: scripts/12_topic_models.py
    title: 12_topic_models.py
    filename: 12_topic_models.py
---

# Single-cell ATAC-seq

In this practical we analyse the chromatin accessibility of single nuclei from the mouse liver:

> Bravo González-Blas C, et al. *Single-cell spatial multi-omics and deep learning dissect enhancer-driven gene regulatory networks in liver zonation.* Nature Cell Biology (2024).
> https://doi.org/10.1038/s41556-023-01316-4

We use the ATAC part of two 10x Multiome samples (GEO [GSE218468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218468)): Mouse_4, prepared with the 10x nuclei isolation protocol, and Mouse_5, prepared with the NST protocol, about 11,000 cells in total. Multiome measures ATAC and RNA in the same nucleus, so the cell types and hepatocyte zones annotated from the RNA in the paper give us reliable labels to check the ATAC analysis. The [single-nucleus RNA-seq practical](../single-nucleus-rna-seq/index.md) analyses single-nucleus RNA-seq data of the same study, and the [MPRA practical](../mpra/index.md) tests which of the accessible regions are active enhancers.

**Research questions:**

1. Can the cell types of the liver be recovered from accessibility alone, without RNA?
2. Is hepatocyte zonation (periportal versus pericentral) visible in the chromatin?
3. Which regions differ between pericentral and periportal hepatocytes, and are they near known zonated genes?

The practical consists of four notebooks:

1. [Part 1: quality control and preprocessing](part_1_qc_preprocessing.ipynb): fragment files, per-cell QC (fragments, TSS enrichment, fragment sizes), doublets, cell-type and zone annotations checked with gene activity, clustering on genome-wide bins without the annotations, pseudobulk peak calling per cluster and a consensus peak set (SnapATAC2, MACS3).
2. [Part 2: topic modelling with pycisTopic](part_2_pycistopic.ipynb): topic models, model selection, batch correction (Harmony), clustering and UMAP, cell-type and zonation topics.
3. [Part 3: dimensionality reduction and clustering with SnapATAC2](part_3_snapatac2.ipynb): spectral embedding, batch correction, clustering, and a comparison with pycisTopic.
4. [Part 4: differentially accessible regions](part_4_dars.ipynb): cell-type DARs, and pericentral versus periportal DARs with SnapATAC2 (per mouse) and pycisTopic (imputed accessibility), compared with each other and with known zonated genes, pseudobulk coverage tracks per cell type and zone (pycisTopic) to look at them in IGV; a comparison with the paper and the conclusion.

## Software

All four notebooks use one mamba environment, `atac-sc`, with SnapATAC2, pycisTopic, scanpy and MACS3; the kernel is "Python (atac-sc)". The commands to create it and register the kernel are given at the start of part 1. Part 2 also uses [Mallet](https://mimno.github.io/Mallet/) for topic modelling, but only if you train models yourself; the notebook then downloads it into `data`. It needs Java, which is included in the environment.

## Data

You work in the folder `data`, next to the notebooks: each notebook writes everything it makes to `data/part_N`. The input data (the **course data**) are already prepared for you on the teaching server, in the read-only folder `/data/teachers/data/scATAC/liver_multiome`; the notebooks find them automatically, so you don't need to copy or download anything. You need about 9 GB of disk space in `data` (most of it for the imported samples of part 1 and the pseudobulk tracks of part 4).

If you run the notebooks somewhere else, they look for the course data in a folder `course_data` next to the notebooks; see below.

### Preparing the course data yourself

To run the notebooks somewhere else than on the teaching server (for example on your own computer or on the VSC), build the course data yourself, in a folder `course_data` next to the notebooks. The scripts put the downloads (the original GEO files, Mallet, the annotation loom) in a second folder, `course_data_build`, which the notebooks don't need. Download the scripts below from the download menu at the top of this page, put them in a folder `scripts` next to the notebooks, and run them from the folder of the notebooks (each script has a Slurm header, and the usage in its first lines):

`01_download_sc.sh` needs `bgzip` and `tabix`, which are in the `atac-cli` environment of the [bulk ATAC-seq practical](../bulk-atac-seq/index.md) (or: `mamba install -n atac-sc -c bioconda htslib`); the other scripts run with `atac-sc`.

| Script | Does | Time |
|---|---|---|
| `01_download_sc.sh` | fragment files (recompressed with bgzip and indexed), the paper's ATAC count matrix, mm10 annotation files; Mallet (for `12_topic_models.sh`) | ~30 min |
| `13_annotations.py` | `course_data/annotations/cell_annotations.tsv` and `paper_cells.tsv` from the public Scope loom of the paper (download instructions in the script) | minutes |
| part 1 | clusters, consensus peaks, cell list and peak matrix; copy `liver_peaks.h5ad`, `consensus_peaks.bed`, `cells.tsv` and `pseudobulk_peaks/` from `data/part_1` to `course_data/data_generated` | ~1 h |
| `12_topic_models.sh` | pycisTopic models with 10 to 50 topics (`course_data/data_generated/topic_models.pkl`) | ~3 h (16 CPU) |
| parts 2 and 3 | copy `cistopic_obj.pkl` from `data/part_2` and `liver_peaks_snapatac2.h5ad` from `data/part_3` to `course_data/data_generated` | |

Building the cisTopic object in part 2 needs about 42 GB of memory. On the teaching server it is provided, so that part 3 and 4 can also run without it.
