# BacTraq

SNP-based clustering with consistent naming across sequencing runs.

## Introduction

BacTraq is a Python package that clusters, tracks, and maintains consistent naming of bacterial clusters during SNP-based analysis. It compares clusters from the current sequencing run with those from previous analyses to determine whether clusters have remained the same or require new naming.

> [!IMPORTANT]
> The name structure follows a cumulative dot-notation rule across threshold levels (Finer threshold name caries information of its coarse threshold)
> A sample assigned cluster 1 at 20 SNPs, cluster 1 at 10 SNPs, and cluster 2 at 5 SNPs is reported as:
>
> | sample  | 20 SNPs | 10 SNPs | 5 SNPs |
> | ------- | ------- | ------- | ------ |
> | sample1 | 1       | 1.1     | 1.1.2  |
>
> If a sample is unclustered at a level, that field is left empty.
> E.g. sample1 clusters at 20 SNPs only → `sample1,1,,`

### Key Features

- Clusters sequences from SNP distance matrices (input as a tab-separated file)
- Compare and handle the changes in clusters' memberships and singleton status (split-to-singleton, dissolved, bridged reclassification) if a clustering history is provided
- Provides consistent cluster naming for easy tracking across routine analyses
- Writes a structured run log with per-cluster naming decisions


## Installation

Create a new environment:

```bash
conda create -n bactraq python=3.10
conda activate bactraq
```

**Install from source:**

```bash
git clone https://github.com/qhgenomics/BacTraq.git
cd BacTraq
pip install .
```

Verify installation:

```bash
bactraq --help
```

```
usage: bactraq [-h] -d FILE -o FILENAME [-t THRESHOLD] [--history HISTORY] [--nRef] [--summary]

SNP-distance clustering with consistent naming across runs.

options:
  -h, --help            show this help message and exit
  -d FILE, --distance-matrix FILE
                        Path to SNPdist matrix (tab-separated). (default: None)
  -o FILENAME, --output FILENAME
                        Output base name. Produces <output>.csv and <output>_history.parquet.gz. (default: None)
  -t THRESHOLD, --threshold THRESHOLD
                        Comma-separated SNP thresholds, largest first. E.g: --threshold 20,10,5 (default: 20,10,5)
  --history HISTORY     Path to history .parquet.gz from a previous run. Generate with `bactraq-history`. Omit for clustering only. (default: None)
  --nRef                No Reference sample in the SNP distance matrix. `Reference` is usually included when you run the workflow with Snippy and Snippy-core. (default: True)
  --summary             Also write <output>_summary.csv: a cluster-change table (one row per history cluster lineage, Status/New Name per threshold, plus singleton samples) comparing
                        this run against --history. No-op without --history. (default: False)
```


## Basic Usage

### First run (no history)

```bash
bactraq -d cc152_snpdists.tsv -o 20250110_cc152 -t 20,10,5
```

Outputs:

- `20250110_cc152.csv` — cluster assignments in dot-notation
- `20250110_cc152_history.parquet.gz` — history file for the next run
- `20250110_cc152_bactraq.log` — full run log

### Subsequent runs (with history)

```bash
bactraq -d cc152_snpdists.tsv -o 20250210_cc152 -t 20,10,5 \
        --history 20250110_cc152_history.parquet.gz
```

### Without a Reference sample

```bash
bactraq -d cc152_snpdists.tsv -o 20250110_cc152 -t 20,10,5 --nRef
```

### Input file

BacTraq takes the default tab-separated output of SNP-dists. The matrix must include a `Reference` row/column unless `--nRef` is passed.


## History File

If you do not have a `.parquet.gz` history file, generate one from a previous cluster CSV using `bactraq-history`:

**Input CSV format:**

```
sample,20 SNPs,10 SNPs,5 SNPs
SS24M01,10,10.1,10.1.1
SS24M02,10,10.2,10.2.1
SS24M03,8,8.1,8.1.1
SS24M04,8,8.1,8.1.1
SS24M05,8,8.1,
```

**Command:**

```bash
bactraq-history 20241108_cc152_cluster.csv -o 20241108_cc152_history -t 20,10,5
```

**Output:**

```
History file: 20241108_cc152_history.parquet.gz
```

## Run Log

Every run writes a `<output>_bactraq.log` file. Each cluster naming decision is logged with an event tag:

```
2026-07-02 14:05:06 [INFO] BacTraq run started — input: new_run_snpdist.tsv
2026-07-02 14:05:06 [INFO] Thresholds: [20, 10, 5]
2026-07-02 14:05:06 [INFO] History file: prev_history.parquet.gz
2026-07-02 14:05:06 [INFO] Clustering complete — 19 samples
2026-07-02 14:05:06 [INFO] [20 SNPs] MERGE  clusters [2, 3] -> new cluster 7; members: [s4, s5, s6, s7]
2026-07-02 14:05:06 [INFO] [20 SNPs] BRIDGED RECLASSIFICATION  [s16] (previously singleton) joined cluster 1; renamed: 1 -> 8; all members: [s1, s16, s2, s3]
2026-07-02 14:05:06 [INFO] [20 SNPs] SPLIT  cluster 4 (had [s10, s8, s9]) -> new cluster 9: [n1, s10]
```

### Note
Check out name rules and example [here](docs/example.md)