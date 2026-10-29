"""Collect the per-sample QC of the bulk pipeline (10_bulk_pipeline.sh, in course_data_build/pipeline) into
course_data/data_generated/qc_summary.tsv, and move the chr17 FASTQs of the two part-1 samples to course_data/fastq_chr17/.

Run with the atac-cli python after all 4 array tasks have finished:
    python scripts/11_bulk_summary.py
"""
import gzip
import json
import shutil
import subprocess
from pathlib import Path

TOPIC = Path(__file__).resolve().parents[1]
BULK = TOPIC / "course_data"           # used by the notebooks
BUILD = TOPIC / "course_data_build"    # pipeline outputs and sample sheet
PART1_SAMPLES = ["FOXA1HNF4A_U1", "FOXA1HNF4A_I1"]


def count(*args):
    return int(subprocess.run(["samtools", "view", "-c", *args], capture_output=True, text=True, check=True).stdout)


def qc(sample):
    d = BUILD / "pipeline" / sample
    fastp = json.load(open(d / f"{sample}.fastp.json"))["summary"]
    pairs = fastp["before_filtering"]["total_reads"] // 2
    kept = fastp["after_filtering"]["total_reads"] / fastp["before_filtering"]["total_reads"]

    aligned = (d / f"{sample}.bowtie2.log").read_text().strip().splitlines()[-1].split("%")[0]

    lines = (d / f"{sample}.markdup_metrics.txt").read_text().splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith("LIBRARY"))
    metrics = dict(zip(lines[i].split("\t"), lines[i + 1].split("\t")))
    dup = float(metrics["PERCENT_DUPLICATION"])

    mapped = {l.split("\t")[0]: int(l.split("\t")[2]) for l in (d / f"{sample}.markdup.idxstats").read_text().splitlines()}
    chrm = mapped.get("chrM", 0) / sum(mapped.values())

    ataqv = json.load(gzip.open(d / f"{sample}.ataqv.json.gz"))[0]["metrics"]

    bam = str(d / f"{sample}.filtered.bam")
    peaks = d / "macs3" / f"{sample}_peaks.narrowPeak"
    fragments = count("-f", "64", bam)
    in_peaks = count("-f", "64", "-L", str(peaks), bam)
    n_peaks = sum(1 for _ in open(peaks))

    return {
        "sample": sample,
        "read_pairs_M": f"{pairs / 1e6:.1f}",
        "kept_after_trim_%": f"{100 * kept:.1f}",
        "aligned_%": aligned,
        "duplicates_%": f"{100 * dup:.1f}",
        "chrM_%": f"{100 * chrm:.1f}",
        "filtered_fragments_M": f"{fragments / 1e6:.1f}",
        "TSS_enrichment": f"{ataqv['tss_enrichment']:.1f}",
        "peaks": n_peaks,
        "FRiP": f"{in_peaks / fragments:.2f}",
    }


samples = [l.split("\t")[0] for l in (BUILD / "samplesheet.tsv").read_text().splitlines()[1:]]
rows = [qc(s) for s in samples]
with open(BULK / "data_generated" / "qc_summary.tsv", "w") as out:
    out.write("\t".join(rows[0]) + "\n")
    for r in rows:
        out.write("\t".join(str(v) for v in r.values()) + "\n")
print((BULK / "data_generated" / "qc_summary.tsv").read_text())

fq = BULK / "fastq_chr17"
fq.mkdir(exist_ok=True)
for s in PART1_SAMPLES:
    for r in ("R1", "R2"):
        src = BUILD / "pipeline" / s / f"{s}_chr17_{r}.fastq.gz"
        if src.exists():
            shutil.move(src, fq / src.name)
print(sorted(p.name for p in fq.iterdir()))
