"""Cell annotations for the scATAC notebooks, from the public Scope loom of Bravo Gonzalez-Blas et al. 2024.

Input:  course_data_build/annotations/Liver-cell_gene-all-2.loom (Scope: Bravo_et_al_Liver)

The loom has no plain download link; Scope's own read-only request fetches it (run from the topic folder):
    mkdir -p course_data_build/annotations
    curl -L -o course_data_build/annotations/Liver-cell_gene-all-2.loom -F "loomFilePath=Bravo_et_al_Liver/Liver-cell_gene-all-2.loom" \
        -F "UUID=Bravo_et_al_Liver" -F "file-type=Loom" https://scope.aertslab.org/upload/
Run with the atac-sc python: python scripts/13_annotations.py
Output: course_data/annotations/cell_annotations.tsv with barcode (as in the fragment files), sample, cell_type, cell_type_fine,
        cell_type_atac and zone; course_data/annotations/paper_cells.tsv with the barcode and sample of the same cells (the
        cells the paper kept, used in part 1 to compare cell selection).

Zone: the hepatocyte states Hep_PP_PC_1..5 are ordered along the periportal -> pericentral axis. Their mean RNA
expression of pericentral (Cyp2e1, Glul, Oat, Axin2, Cyp1a2, Slc1a2) minus periportal (Cyp2f2, Hal, Sds, Pck1, Arg1,
Hsd17b13) markers increases monotonically from state 1 to 5 (z-scored scores -2.16, -1.38, 0.14, 0.94, 2.46), so:
states 1-2 = PP, 3 = mid, 4-5 = PC.
"""
from pathlib import Path

import h5py
import pandas as pd

SC = Path(__file__).resolve().parents[1] / "course_data"
BUILD = Path(__file__).resolve().parents[1] / "course_data_build"
SAMPLES = {"TEW__043783__Mouse_4_Multiome_10x": "Mouse_4_10x", "TEW__ebb273__Mouse_5_Multiome_NST": "Mouse_5_NST"}
ZONE = {"Hep_PP_PC_1": "PP", "Hep_PP_PC_2": "PP", "Hep_PP_PC_3": "mid", "Hep_PP_PC_4": "PC", "Hep_PP_PC_5": "PC"}
BROAD = {"Kupffer_cycling": "Kupffer", "cDC": "DC", "pDC": "DC", "VEC_PP": "VEC", "VEC_PC": "VEC", "VEC_LY": "VEC"}


def decode(values):
    return [v.decode() if isinstance(v, bytes) else v for v in values]


with h5py.File(f"{BUILD}/annotations/Liver-cell_gene-all-2.loom", "r") as loom:
    attrs = loom["col_attrs"]
    cells = pd.DataFrame({key: decode(attrs[key][:]) for key in
                          ["CellID", "sample_id", "Refined_cell_type", "Refined_cell_type_ATAC"]})

cells = cells[cells["sample_id"].isin(SAMPLES)]
table = pd.DataFrame({
    "barcode": cells["CellID"].str.extract(r"^([ACGT]+-1)")[0],
    "sample": cells["sample_id"].map(SAMPLES),
    "cell_type": cells["Refined_cell_type"].map(lambda c: "Hepatocyte" if c.startswith("Hep_") else BROAD.get(c, c)),
    "cell_type_fine": cells["Refined_cell_type"],
    "cell_type_atac": cells["Refined_cell_type_ATAC"],
    "zone": cells["Refined_cell_type"].map(ZONE).fillna(""),
})
assert table["barcode"].notna().all()
(SC / "annotations").mkdir(parents=True, exist_ok=True)
table.to_csv(f"{SC}/annotations/cell_annotations.tsv", sep="\t", index=False)
table[["barcode", "sample"]].to_csv(f"{SC}/annotations/paper_cells.tsv", sep="\t", index=False)
print(table.groupby("sample").size())
print(table["cell_type"].value_counts())
print(table["zone"].value_counts())
