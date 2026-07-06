import pytest
import pandas as pd
import numpy as np


def make_sample_df(assignments: dict, threshold: list) -> pd.DataFrame:
    """Build a cluster-number DataFrame as cluster_from_distance_matrix would return."""
    cols = {f"cluster{t}": v for t, v in zip(threshold, assignments.values())}
    df = pd.DataFrame(cols)
    df.index.name = "sample"
    return df


def make_history_df(tuples_by_sample: dict, columns: list) -> pd.DataFrame:
    """Build a history DataFrame with tuple values, as assign_cluster_database expects."""
    df = pd.DataFrame(tuples_by_sample, index=pd.Index(list(tuples_by_sample.keys()), name="sample"))
    df.columns = columns
    return df
