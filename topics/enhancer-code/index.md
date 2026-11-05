# Enhancer code

In this practical we study the **regulatory code**: the transcription factor binding sites inside enhancers that decide in which cell type a gene is switched on. We use a deep learning model trained on the single-cell ATAC-seq atlas of:

> Zemke NR, Armand EJ, Wang W, et al. *Conserved and divergent gene regulatory programs of the mammalian neocortex.* Nature (2023) 624:390–402. https://doi.org/10.1038/s41586-023-06819-6

The authors profiled gene expression and chromatin accessibility in single nuclei of the **primary motor cortex** of four mammals (human, macaque, marmoset and mouse) and grouped them into cell types. We use the mouse data through the deep learning models and tools of CREsted (Kempynck N, De Winter S, et al. *CREsted: modeling genomic and synthetic cell-type-specific enhancers across tissues and species.* Nature Methods (2026). https://doi.org/10.1038/s41592-026-03057-2). The consensus peaks and the trained models are downloaded from the CREsted resources; the precomputed contribution scores come from the TF-MINDI tutorial dataset on [Zenodo](https://doi.org/10.5281/zenodo.18757793).

An astrocyte, an oligodendrocyte, a microglium and a neuron share one genome but use different enhancers. The question running through both notebooks is what it is about the *sequence* of an enhancer that makes it work in one of these cell types and not the others — first by reading that code out of a trained model, then by writing new enhancers with it.

The practical consists of two notebooks:

1. [Part 1: reading the regulatory code of brain cell types](part_1_regulatory_code.ipynb) (Python, TF-MINDI): starting from precomputed contribution scores for 5,700 enhancers across all 19 cell types, cutting out the stretches of DNA the model relies on (seqlets), matching them to a motif collection, assigning a transcription factor family and clustering them into motifs, building the cell-type-by-motif **code table**, naming candidate transcription factors by correlating motif use with gene expression from the same nuclei, and finally asking whether whole enhancers group by cell type when described only by the motifs they contain.
2. [Part 2: designing synthetic enhancers](part_2_enhancer_design.ipynb) (Python, CREsted): scoring the atlas's consensus peaks to get a reference set of real enhancers, evolving a microglia enhancer from random sequence and following it mutation by mutation, converting an enhancer of one cell type into another, building enhancers from the motifs discovered in part 1, and checking the designs against a second model.

The two notebooks use different tools and are deliberately separate: part 1 analyses the code, part 2 writes with it. Run part 1 first, because part 2 builds its designs from the motif clusters part 1 saves in `results/`. Training a sequence model and computing contribution scores is covered in the CREsted part of the course; here both are taken as given.

## Software

Both notebooks use the same mamba environment, separate from the `regulatory_genomics` environment of the other practicals because the deep learning libraries need their own package versions. The exact commands (including how to register the environment as a Jupyter kernel) are given at the start of each notebook.

    mamba create -n regulatory_genomics_enhancers -c conda-forge python=3.12 pip ipykernel \
        ipywidgets numpy pandas matplotlib seaborn "scanpy>=1.12.5" anndata
    pip install tensorflow "crested[motif]" "tfmindi>=2"

| Notebook | Environment | Kernel |
|---|---|---|
| Part 1 | `regulatory_genomics_enhancers` | Python (regulatory_genomics_enhancers) |
| Part 2 | `regulatory_genomics_enhancers` | Python (regulatory_genomics_enhancers) |

[TF-MINDI](https://tf-mindi.readthedocs.io/en/latest/) carries part 1: it groups the important sequence stretches (seqlets) into transcription factor motifs, annotates them and builds the enhancer code table. It must be version 2 or newer. [CREsted](https://crested.readthedocs.io/en/stable/) carries part 2: loading the enhancer models, predictions and enhancer design. The environment needs `scanpy>=1.12.5`: older versions crash when clustering more than a few thousand seqlets with recent SciPy.

Everything runs on a normal CPU; a GPU is not needed. Part 1 takes about 20 minutes end to end and needs roughly 8 GB of free memory, because comparing 51,000 seqlets with 4,800 motifs is done in one large matrix. In part 2 the two slowest steps (scoring a sample of peaks, evolving an enhancer) take about 10 and 30 minutes, and both are cached.

## Data

All data is stored in the folder `data`, next to the notebooks, and is downloaded by the notebooks themselves.

Part 1 needs one archive (1.85 GB), published with the TF-MINDI tutorial, which contains three things: precomputed contribution scores for the 2,000 most specific regions of each of the 19 cell types, the single-cell RNA-seq of the same nuclei, and a table of mouse transcription factors (AnimalTFDB4). Computing those contribution scores from scratch would take many hours on a CPU.

Part 2 needs the consensus ATAC-seq peaks of the atlas (26 MB), the DeepBICCN2 and DeepBICCN models (60 MB each) and the mouse reference genome mm10 (about 900 MB compressed, 2.7 GB unpacked).

Part 1 also downloads the v11 motif collection (64 MB) from the Aerts lab cistarget resources, which it uses to recognise motifs and assign transcription factor families.

Together this needs roughly 6 GB of disk space.
