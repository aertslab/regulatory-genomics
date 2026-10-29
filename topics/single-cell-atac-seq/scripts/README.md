# Single-cell ATAC-seq: preparing the course data

Notes for teachers: how the data of the notebooks are built and installed on the teaching server. Students don't need this; the topic page (`../index.md`) tells them where the data are.

## Layout

Three folders next to the notebooks, all git-ignored:

- `course_data/`: the **course data**, everything the notebooks read and nothing else (about 6 GB). On the teaching server the notebooks read it from the shared, read-only folder `/data/teachers/data/scATAC/liver_multiome`; when that folder doesn't exist, they use `course_data/` next to the notebooks. The setup cell of each notebook chooses between the two (variable `SC`).
- `course_data_build/`: downloads and intermediate files, only needed to (re)build `course_data/`: the original GEO fragment files and RNA count matrix, the Scope looms, Mallet, logs, and `old_v1_peaks_from_rna_labels/` (results of an earlier version that called peaks per RNA cell type; not used).
- `data/`: the working folder of the notebooks; each notebook writes to `data/part_N` (one per student, writable).

## Installing on the teaching server

Copy `course_data/` as a whole:

```bash
rsync -a course_data/ /data/teachers/data/scATAC/liver_multiome/
```

What it contains:

| Path | Contents | Size |
|---|---|---|
| `fragments/` | fragment files of the two mice, bgzip-compressed and indexed | 3.4 GB |
| `genome/` | mm10 chromosome sizes, blacklist, GENCODE vM25 annotation, TSS of protein-coding genes | 29 MB |
| `annotations/cell_annotations.tsv`, `annotations/paper_cells.tsv` | RNA-based cell types and zones; the cells the paper kept | 1 MB |
| `raw/GSE218468_multiome_atac_counts.tsv.gz` | the paper's count matrix; part 1 downloads it from GEO into `data` and only uses this copy if GEO cannot be reached | 169 MB |
| `data_generated/` | pre-generated results, so that each notebook can also run on its own: `cells.tsv`, `consensus_peaks.bed`, `liver_peaks.h5ad`, `pseudobulk_peaks/` (part 1), `topic_models.pkl`, `cistopic_obj.pkl` (part 2), `liver_peaks_snapatac2.h5ad` (part 3) | 1.8 GB |

## Rebuilding the course data

The steps are on the topic page (`../index.md`, section "Preparing the course data yourself"), so that students can also follow them. The scripts write the files the notebooks use to `course_data/` and everything else to `course_data_build/`.
