# BacTraq

SNP-based clustering with consistent naming.

## Intro:

BacTraq is a Python package designed to track and maintain consistent naming of clusters during SNP-based analysis. It compares clusters from the current sequencing run with those from previous analyses to determine whether clusters have remained the same or require new naming.

### Key Features

- Detects and handles overlapping clusters between sequencing runs.

- Provides consistent cluster naming for easy tracking across analyses.

- Identifies merge and split events for clusters.

### Cluster Naming Rules
BacTraq evaluates pairwise overlaps between current and previous clusters at each threshold level, applying the following rules:

1. Merge Case:

If a current cluster overlaps with multiple previous clusters, it is considered a merge.

A new name is assigned to the merged cluster.

[merge image](docs/merge_case.PNG)

2. Split Case:

If multiple current clusters overlap with a single previous cluster, it is considered a split.

New names are assigned to each split cluster.

[split image](docs/split.PNG)

3. Entirely New Cluster Case:

If a current cluster contains only new samples and no overlap with previous clusters, it is considered entirely new.

A new name is assigned to the cluster.

[entire new image](docs/entire_new.PNG)

4. Consistent Naming Case:

If a current cluster overlaps with only one previous cluster, including both shared and new samples, the previous name is retained.

[consistent naming image](docs/consitent_name.PNG)
