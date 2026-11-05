"""Text blocks shared by both enhancer-code notebooks."""

INSTALL = """
### Selecting the kernel and installing the Python environment

This practical uses its own **Conda/Mamba environment**, `regulatory_genomics_enhancers`. It is separate from the `regulatory_genomics` environment of the single-nucleus and spatial practicals because the deep learning libraries (TensorFlow, CREsted, TF-MINDI) need their own, specific package versions. You only need to create it **once**, and the same environment is used for part 1 and part 2.

Run the following in a **terminal on the computer/VSC cluster where you use Jupyter**, not in a notebook code cell:

```bash
mamba create -n regulatory_genomics_enhancers -c conda-forge python=3.12 pip ipykernel ipywidgets numpy pandas matplotlib seaborn "scanpy>=1.12.5" anndata
conda activate regulatory_genomics_enhancers
pip install tensorflow "crested[motif]" "tfmindi>=2"
```

- **TensorFlow** is the deep learning library that runs the neural network.
- **CREsted** loads trained enhancer models, makes predictions, explains them and designs new sequences.
- **TF-MINDI** groups the important pieces of DNA found by the model into transcription factor motifs, and annotates each motif with the transcription factor family that binds it. Version 2 or newer is required.

The `scanpy>=1.12.5` is not optional: older versions crash when clustering more than a few thousand seqlets with recent SciPy.

The installation may take 10–15 minutes. Then **register a Jupyter kernel**:

```bash
conda activate regulatory_genomics_enhancers
python -m ipykernel install --user --name regulatory_genomics_enhancers --display-name "Python (regulatory_genomics_enhancers)"
jupyter kernelspec list
```

Open the notebook in VS Code or Jupyter and choose **Python (regulatory_genomics_enhancers)** via **Select Kernel** (VS Code) or **Kernel → Change Kernel** (Jupyter).

Everything in this practical runs on a normal laptop CPU; a GPU is not needed (though it will speed things up). The slowest steps take 10–20 minutes and are marked in the text.
"""

DEEPBICCN = """
### The DeepBICCN2 model: where it comes from and what it predicts

We do not train a model in this practical (training is covered in a separate session). Instead we use **DeepBICCN2**, a model published with CREsted and distributed with it.

Watch the digit. CREsted also ships a model called **DeepBICCN**, and the two are easy to confuse. They have the same architecture and the same 6.3 million parameters, and were trained on the same atlas — but on two different ways of summarising the same ATAC-seq data. `DeepBICCN` predicts the **mean coverage** over the middle of a peak; `DeepBICCN2` predicts **Tn5 cut-site counts**, the positions where the enzyme actually cut. We use `DeepBICCN2` here, and the precomputed contribution scores of part 1 were made with it too.

**What goes in.** A DNA sequence of exactly 2,114 base pairs.

**What comes out.** 19 numbers, one per cortex cell type:

| Group | Cell types |
|---|---|
| Non-neuronal | `Astro` (astrocytes), `Oligo` (oligodendrocytes), `OPC`, `Micro_PVM` (microglia), `Endo`, `VLMC` |
| Excitatory neurons | `L2_3IT`, `L5IT`, `L6IT`, `L5ET`, `L5_6NP`, `L6CT`, `L6b` |
| Inhibitory neurons | `Pvalb`, `Sst`, `SstChodl`, `Vip`, `Lamp5`, `Sncg` |

Each number is the **accessibility the model expects in that cell type**: how many Tn5 cut sites (the ATAC-seq signal) it predicts over the central 1,000 bp of the sequence. Unlike a probability these values are not capped at 1; a strong, specific enhancer can score 30 or more in its own cell type and near zero everywhere else.

**How it was trained (in short).** Training means showing the model many sequences together with the right answer and adjusting its 6.3 million internal parameters until its predictions match. DeepBICCN2 was trained in two stages:

1. First on **all** consensus peaks of a single-cell ATAC-seq atlas of the mouse primary motor cortex, learning the general relationship between sequence and accessibility.
2. Then **fine-tuned** on the subset of regions that are specific to one cell type, which sharpens exactly the distinction we care about here.

Crucially, the model only ever saw **DNA sequence**: no gene names, no positions in the genome, no information about which TFs exist. If it predicts cell-type-specific accessibility correctly, it must have learned which sequence features (TF binding sites) make a region accessible in a microglium and not in a neuron. Reading that knowledge back out is what this practical is about.

> Kempynck N, De Winter S, et al. *CREsted: modeling genomic and synthetic cell-type-specific enhancers across tissues and species.* Nature Methods (2026). The data come from Zemke NR, Armand EJ, Wang W, et al. *Conserved and divergent gene regulatory programs of the mammalian neocortex.* Nature (2023) 624:390–402, which profiled the primary motor cortex.
"""

PALETTE_TEXT = """
### A colour convention for the whole practical

Before we plot anything, we fix **one colour per cell type** and use it in every figure of both notebooks. This is worth doing deliberately rather than letting each plotting function pick its own defaults: once you have learned that orange-red means microglia and blue means oligodendrocytes, every later figure becomes readable at a glance.

We also fix the **order** here — non-neuronal, excitatory neurons, inhibitory interneurons — and reuse it in every table and heatmap, so that related cell types always sit next to each other.
"""

PALETTE = '''
# One colour per cell type, used in every figure of both notebooks.
# The order is the one we use everywhere: non-neuronal, then excitatory, then inhibitory.
CELL_TYPE_COLOURS = {
    # non-neuronal
    "Astro":     "#009E73",
    "Oligo":     "#0072B2",
    "OPC":       "#543005",
    "Micro_PVM": "#D55E00",
    "Endo":      "#332288",
    "VLMC":      "#66A61E",
    # excitatory neurons
    "L2_3IT":    "#CC79A7",
    "L5IT":      "#003C30",
    "L6IT":      "#DDCC77",
    "L5ET":      "#B2182B",
    "L5_6NP":    "#E7298A",
    "L6CT":      "#C49C94",
    "L6b":       "#EE8866",
    # inhibitory interneurons
    "Pvalb":     "#E69F00",
    "Sst":       "#882255",
    "SstChodl":  "#A6761D",
    "Vip":       "#56B4E9",
    "Lamp5":     "#BBCC33",
    "Sncg":      "#8C564B",
}
GREY = "#CCCCCC"  # anything without a colour of its own

CELL_TYPES = list(CELL_TYPE_COLOURS)          # all 19, in the order above
CELL_TYPE_GROUPS = {"non-neuronal": CELL_TYPES[:6],
                    "excitatory": CELL_TYPES[6:13],
                    "inhibitory": CELL_TYPES[13:]}

def cell_colour(cell_type):
    """The agreed colour of a cell type."""
    return CELL_TYPE_COLOURS.get(cell_type, GREY)
'''

CONTRIB_HELPER = '''
from matplotlib.ticker import MultipleLocator

def motif_span(scores, one_hot, class_index=0, sequence_index=0,
               pad=30, min_span=200, max_span=320):
    """The stretch of the sequence that actually carries contribution.

    Returns (start, end) positions in the sequence, so we can zoom a plot onto the
    motifs instead of showing 2,114 bp of mostly-empty track. If the important
    positions are spread over more than `max_span`, we keep the `max_span` window
    that carries the most contribution: a window wide enough to hold every motif is
    usually too wide to read the letters in.
    """
    length = scores.shape[2]
    track = np.abs((scores[sequence_index, class_index] * one_hot[sequence_index]).sum(axis=1))
    important = np.where(track > 0.1 * track.max())[0]
    start, end = int(important.min()) - pad, int(important.max()) + pad
    if end - start > max_span:
        densest = int(np.argmax(np.convolve(track, np.ones(max_span), mode="valid")))
        start, end = densest, densest + max_span
    if end - start < min_span:                      # never zoom in further than min_span
        centre = (start + end) // 2
        start, end = centre - min_span // 2, centre + min_span // 2
    return max(0, start), min(length, end)

def plot_contributions(scores, one_hot, labels, class_labels, span=None, **kwargs):
    """Draw contribution logos.

    `span` is a (start, end) pair in sequence coordinates; without it the whole
    sequence is drawn. The x axis always shows the real position in the sequence,
    so a zoomed plot starts at the position it zoomed to rather than at 0.
    """
    length = scores.shape[2]
    start, end = span if span is not None else (0, length)
    zoom = end - start
    x_shift = start - length // 2 + zoom // 2
    fig, axes = crested.pl.explain.contribution_scores(
        scores, one_hot,
        sequence_labels=labels, class_labels=class_labels,
        zoom_n_bases=zoom, x_shift=x_shift,
        coordinates=(0, length),   # label the x axis with positions in the sequence
        show=False, **kwargs,
    )
    tick = max(10, int(round(zoom / 8 / 10)) * 10)   # aim for about eight ticks
    for ax in np.atleast_1d(axes):
        ax.xaxis.set_major_locator(MultipleLocator(tick))
    plt.show()
'''
