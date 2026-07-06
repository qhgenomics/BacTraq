## Cluster Naming Rules

BacTraq evaluates pairwise overlaps between current and previous clusters at each threshold level:

![how BacTraq traces and renames clusters across runs](cluster-trace-diagram.svg)

The diagram above walks through the worked example in [`example/`](example): `prev_clusters.csv`
(previous run) traced forward into `new_run_cluster.csv` (current run), matching the events in
`new_run_cluster_bactraq.log`. Each current-run box shows the name BacTraq assigned; the small gray
line underneath it is what a history-blind re-clustering pass would have called the same samples —
a fresh, unrelated number every run, with no way to tell it's the same cluster as before.

1. **Consistent / Expand**

A current cluster overlaps with exactly one previous cluster → previous name retained.

2. **Split**

Multiple current clusters overlap with a single previous cluster → new names assigned to each.

3. **Merge**

A current cluster overlaps with multiple previous clusters → new name assigned.

4. **Entirely New**

Current cluster has no overlap with any previous cluster → new name assigned.

5. **Split-to-Singleton** *(new in v1.1.0)*

Members of a previous cluster became singletons this run -> treated a split case and rename the cluster.

6. **Dissolved** *(new in v1.1.0)*

All members of a previous cluster became singletons → cluster logged as dissolved.

7. **Bridged Reclassification** *(new in v1.1.0)*

A previously-singleton sample joins an existing cluster → cluster renamed to reflect the new composition.