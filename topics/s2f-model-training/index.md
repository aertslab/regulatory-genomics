# s2f model training

In this practical we learn how to train and use **sequence-to-function** models. Sequence-to-function models take the DNA sequence as input and predict some kind of 'function' of that sequence, such as scATAC-seq data. Training and investigating such models provides valuable insights into how our genomes work in different cell types.

> Zemke NR, Armand EJ, Wang W, et al. *Conserved and divergent gene regulatory programs of the mammalian neocortex.* Nature (2023) 624:390–402. https://doi.org/10.1038/s41586-023-06819-6

The authors profiled gene expression and chromatin accessibility in single nuclei of the **primary motor cortex** of four mammals (human, macaque, marmoset and mouse) and grouped them into cell types. We use the mouse data through the deep learning models and tools of CREsted (Kempynck N, De Winter S, et al. *CREsted: modeling genomic and synthetic cell-type-specific enhancers across tissues and species.* Nature Methods (2026). https://doi.org/10.1038/s41592-026-03057-2). The consensus peaks and the trained models are downloaded from the CREsted resources; the precomputed contribution scores come from the TF-MINDI tutorial dataset on [Zenodo](https://doi.org/10.5281/zenodo.18757793).

How does each of the cell types in these dataset come to be different in their morphology and function, even though they all originated from one single cell and they all have the same DNA? We'll train a model here to explore exactly that question.

The practical consists of two notebooks:

1. [sequence to function modelling with CREsted](sequence_to_function_model_with_crested.ipynb) (Python, CREsted): we will start with bigwig files derived from the analysis of a scATAC-seq dataset and train a s2f model on that using the CREsted package. We will also learn how to evaluate the performance of that model and how to use it for initial downstream tasks such as making predictions and looking at attributions. In the next topic (enhancer code), we will use what we learned here to dive deeper into the biology.
2. [chromatin accessibility at base-pair resolution with ChromBPNet](sequence_to_function_model_with_chrombpnet.ipynb) (Python, ChromBPNet): we train a ChromBPNet model on the microglia of the same dataset. ChromBPNet predicts where Tn5 inserts, base by base, and separates the sequence preference of Tn5 from the signal of transcription factors. We evaluate the model, find the motifs it uses with TF-MoDISco, look at TF footprints, and mutate the FIRE enhancer in silico.

## Software

This notebook uses its own mamba environment, separate from the `regulatory_genomics` environment of the other practicals (expect for enhancer code topic) because the deep learning libraries need their own package versions. If you have some problems, do check the CREsted repository as there are relevant information there too. Some packages we install here will only be used later on (in enhancer code topic) but for the sake of simplicity, you can create the following environment that you can use for both topics:

    mamba create -n regulatory_genomics_enhancers -c conda-forge python=3.12 pip ipykernel \
        ipywidgets numpy pandas matplotlib seaborn "scanpy>=1.12.5" anndata
    pip install tensorflow "crested[motif]" "tfmindi>=2"

ChromBPNet needs older package versions (TensorFlow 2.8) and gets its own environment, `chrombpnet`. The exact commands (including how to register the environment as a Jupyter kernel) are given at the start of the ChromBPNet notebook.

| Notebook | Environment | Kernel |
|---|---|---|
| S2F with CREsted | `regulatory_genomics_enhancers` | Python (regulatory_genomics_enhancers) |
| S2F with ChromBPNet | `chrombpnet` | Python (chrombpnet) |

To be able to train models, you will need a GPU as it just takes too much time with a CPU. Inference tasks can usually be done on a CPU, though for big amounts of data that could also be quite time consuming.

## Data

The commands you need to download the data and the trained model is included in the notebook. They are commented out since they take a while to download and untar, but in your first time you'll need to run them yourself and change the paths accordingly (necessary explanations are in the notebooks). 

You can adapt this notebook to your own data, necessary explanations are in the notebook but in principle you just need to set the paths to your own data!

The ChromBPNet notebook downloads everything it needs into the folder `data`, next to the notebook: chromosomes 9 and 18 of the mouse genome from UCSC, and one file (about 340 MB) with the trained models, peaks and the precomputed TF-MoDISco and footprint results. A second, optional file (about 500 MB) contains the full microglia fragment file, needed only to rerun the training (on a GPU; see Step 4 of the notebook). That file was made with the scripts in `scripts/chrombpnet/` (`00_download_refs.sh` to `09_bundle.sh`): the microglia fragments from the raw reads in the Sequence Read Archive, and the models on an A100 GPU, about 20 hours in total, most of it downloading, aligning and the interpretation steps of the ChromBPNet pipeline.
