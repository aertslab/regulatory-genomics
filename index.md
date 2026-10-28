# Regulatory genomics

This website contains the practicals of the course regulatory genomics. Each topic consists of one or more Jupyter notebooks in which we analyse data from a published study, starting from the raw sequencing data.

## Topics

- [Bulk RNA-seq](topics/bulk-rna-seq/index.md): can the transcription factors FOXA1 and HNF4A activate liver genes in a blood cell line? From sequencing reads to differential expression and biological interpretation.

- [Bulk ATAC-seq](topics/bulk-atac-seq/index.md): how does chromatin accessibility in the liver change in a mouse model of diabetes and high-fat diet? From sequencing reads to peaks, quality control, consensus peaks and differential accessibility with replicates (MACS3, DiffBind, DESeq2).

- [Single-nucleus RNA-seq](topics/single-nucleus-rna-seq/index.md): how can we identify the different cell populations of the mouse liver and investigate hepatocyte transcriptional heterogeneity? From gene-expression counts to quality control, clustering, cell-type annotation and liver zonation.

- [Single-cell ATAC-seq](topics/single-cell-atac-seq/index.md): can liver cell types and hepatocyte zonation be recovered from chromatin accessibility alone? Quality control, consensus peaks, topic modelling (pycisTopic), spectral embedding (SnapATAC2) and differentially accessible regions.

- [Spatial transcriptomics](topics/spatial-transcriptomics/index.md): how is gene expression organized within mouse liver tissue? Comparing Visium, MERFISH and Resolve Molecular Cartography to explore spatial expression patterns, image segmentation, liver zonation and spatial neighborhoods.

## How to use these notebooks

The notebooks are meant to be read and run in your own time. Each topic page explains which software you need and how to install it. Some steps (downloading the data, aligning reads) take a long time, the expected run time is mentioned in the notebooks.
