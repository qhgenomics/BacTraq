import pytest
import pandas as pd
import numpy as np
import os

from BacTraq.bactraq import (
    get_name,
    get_parent,
    by_suffix,
    by_last,
    new_in_first_list,
    increment_tuple_cluster,
    match_db_cluster,
    intesection_betweet_cluster_dicts,
    cluster_number_generate,
    detect_singleton_transitions,
    _inject_split_if_affected,
    assign_cluster_database,
    pretty_name_save,
    cluster_from_distance_matrix,
)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

class TestGetName:
    def test_single_level(self):
        assert get_name((1,)) == "1"

    def test_two_levels(self):
        assert get_name((1, 2)) == "1.2"

    def test_three_levels(self):
        assert get_name((1, 2, 3)) == "1.2.3"


class TestGetParent:
    def test_root_level(self):
        assert get_parent((1,)) == "root"

    def test_second_level(self):
        assert get_parent((1, 2)) == "1"

    def test_third_level(self):
        assert get_parent((1, 2, 3)) == "1.2"


class TestBySuffixByLast:
    def test_by_suffix(self):
        assert by_suffix((1, 2, 3)) == (1, 2)

    def test_by_last(self):
        assert by_last((1, 2, 3)) == 3


class TestNewInFirstList:
    def test_basic_difference(self):
        result = new_in_first_list([1, 2, 3], [2, 3])
        assert result == [1]

    def test_no_difference(self):
        assert new_in_first_list([1, 2], [1, 2]) == []

    def test_empty_list1(self):
        assert new_in_first_list([], [1, 2]) == []

    def test_all_new(self):
        result = set(new_in_first_list([1, 2], []))
        assert result == {1, 2}

    def test_column_name_strings(self):
        result = new_in_first_list(["20 SNPs", "10 SNPs"], ["20 SNPs"])
        assert result == ["10 SNPs"]


class TestIncrementTupleCluster:
    def test_single_level(self):
        assert increment_tuple_cluster((1,)) == (2,)

    def test_two_levels(self):
        assert increment_tuple_cluster((1, 2)) == (1, 3)

    def test_three_levels(self):
        assert increment_tuple_cluster((1, 2, 3)) == (1, 2, 4)


# ---------------------------------------------------------------------------
# intesection_betweet_cluster_dicts
# ---------------------------------------------------------------------------

class TestIntersectionBetweenClusterDicts:
    def test_full_overlap(self):
        query = {(1,): ["s1", "s2"]}
        db = {(1,): ["s1", "s2"]}
        result = intesection_betweet_cluster_dicts(query, db)
        assert result.loc[(1,), (1,)] == True  # noqa: E712  (np.True_ != True with `is`)

    def test_partial_overlap(self):
        query = {(1,): ["s1", "s2"]}
        db = {(1,): ["s1", "s3"]}
        result = intesection_betweet_cluster_dicts(query, db)
        assert result.loc[(1,), (1,)] == True  # noqa: E712

    def test_no_overlap(self):
        query = {(1,): ["s1"]}
        db = {(2,): ["s2"]}
        result = intesection_betweet_cluster_dicts(query, db)
        assert result.loc[(1,), (2,)] == False  # noqa: E712

    def test_multiple_clusters(self):
        query = {(1,): ["s1", "s2"], (2,): ["s3"]}
        db = {(1,): ["s1"], (2,): ["s3"]}
        result = intesection_betweet_cluster_dicts(query, db)
        assert result.loc[(1,), (1,)] == True  # noqa: E712
        assert result.loc[(1,), (2,)] == False  # noqa: E712
        assert result.loc[(2,), (2,)] == True  # noqa: E712


# ---------------------------------------------------------------------------
# match_db_cluster
# ---------------------------------------------------------------------------

class TestMatchDbCluster:
    def _make_row(self, bool_values, db_clusters):
        return pd.Series(bool_values, index=db_clusters)

    def test_single_match(self):
        row = self._make_row([True], [(1,)])
        result = match_db_cluster(row, [(1,)])
        assert list(result) == [(1,)]

    def test_no_match(self):
        row = self._make_row([False], [(1,)])
        result = match_db_cluster(row, [(1,)])
        assert result is None

    def test_multiple_matches(self):
        row = self._make_row([True, True], [(1,), (2,)])
        result = match_db_cluster(row, [(1,), (2,)])
        assert set(result) == {(1,), (2,)}


# ---------------------------------------------------------------------------
# cluster_number_generate
# ---------------------------------------------------------------------------

class TestClusterNumberGenerate:
    def _make_df(self, cluster_cols: dict) -> pd.DataFrame:
        df = pd.DataFrame(cluster_cols)
        df.index.name = "sample"
        return df

    def test_single_threshold_no_singletons(self):
        df = self._make_df({"cluster20": {"s1": 1, "s2": 1, "s3": 2, "s4": 2}})
        result = cluster_number_generate(df, [20])
        assert result.loc["s1", "20 SNPs"] == (1,)
        assert result.loc["s2", "20 SNPs"] == (1,)
        assert result.loc["s3", "20 SNPs"] == (2,)
        assert result.loc["s4", "20 SNPs"] == (2,)

    def test_single_threshold_with_singleton(self):
        df = self._make_df({"cluster20": {"s1": 1, "s2": 1, "s3": 2}})
        result = cluster_number_generate(df, [20])
        assert result.loc["s1", "20 SNPs"] == (1,)
        assert result.loc["s2", "20 SNPs"] == (1,)
        assert result.loc["s3", "20 SNPs"] is None

    def test_all_singletons(self):
        df = self._make_df({"cluster20": {"s1": 1, "s2": 2, "s3": 3}})
        result = cluster_number_generate(df, [20])
        assert result.loc["s1", "20 SNPs"] is None
        assert result.loc["s2", "20 SNPs"] is None
        assert result.loc["s3", "20 SNPs"] is None

    def test_two_thresholds_hierarchical(self):
        # s1,s2,s3 in same 20-SNP cluster; s1,s2 in same 10-SNP sub-cluster; s3 is sub-singleton
        df = self._make_df({
            "cluster20": {"s1": 1, "s2": 1, "s3": 1, "s4": 2},
            "cluster10": {"s1": 1, "s2": 1, "s3": 2, "s4": 3},
        })
        result = cluster_number_generate(df, [20, 10])
        # s4: singleton at 20 SNP level → None at both levels
        assert result.loc["s4", "20 SNPs"] is None
        assert result.loc["s4", "10 SNPs"] is None
        # s1,s2 share sub-cluster → (1,) and (1,1)
        assert result.loc["s1", "20 SNPs"] == (1,)
        assert result.loc["s2", "20 SNPs"] == (1,)
        assert result.loc["s1", "10 SNPs"] == (1, 1)
        assert result.loc["s2", "10 SNPs"] == (1, 1)
        # s3: in 20-SNP cluster 1 but sub-singleton at 10 SNPs
        assert result.loc["s3", "20 SNPs"] == (1,)
        assert result.loc["s3", "10 SNPs"] is None

    def test_raises_on_none_threshold(self):
        df = self._make_df({"cluster20": {"s1": 1}})
        with pytest.raises(AttributeError):
            cluster_number_generate(df, None)

    def test_cluster_numbers_are_sequential(self):
        # 3 clusters, non-contiguous sklearn numbering → should be renumbered 1,2,3
        df = self._make_df({"cluster20": {"s1": 3, "s2": 3, "s3": 7, "s4": 7, "s5": 11, "s6": 11}})
        result = cluster_number_generate(df, [20])
        names = {result.loc[s, "20 SNPs"] for s in ["s1","s2","s3","s4","s5","s6"]}
        assert names == {(1,), (2,), (3,)}


# ---------------------------------------------------------------------------
# detect_singleton_transitions
# ---------------------------------------------------------------------------

class TestDetectSingletonTransitions:
    def _make_df(self, vals: dict) -> pd.DataFrame:
        df = pd.DataFrame({"20 SNPs": vals})
        df.index.name = "sample"
        return df

    def test_to_singleton(self):
        df = self._make_df({"s1": (1,), "s2": None})
        db = self._make_df({"s1": (1,), "s2": (1,)})
        to_s, from_s, affected = detect_singleton_transitions(df, db, "20 SNPs")
        assert to_s == {"s2": (1,)}
        assert from_s == {}
        assert affected == {(1,)}

    def test_from_singleton(self):
        df = self._make_df({"s1": (1,), "s2": (1,)})
        db = self._make_df({"s1": (1,), "s2": None})
        to_s, from_s, affected = detect_singleton_transitions(df, db, "20 SNPs")
        assert to_s == {}
        assert from_s == {"s2": (1,)}
        assert affected == set()

    def test_no_transitions(self):
        df = self._make_df({"s1": (1,), "s2": (1,)})
        db = self._make_df({"s1": (1,), "s2": (1,)})
        to_s, from_s, affected = detect_singleton_transitions(df, db, "20 SNPs")
        assert to_s == {}
        assert from_s == {}
        assert affected == set()

    def test_new_sample_ignored(self):
        # s3 is in current but not in history — should not appear in transitions
        df = self._make_df({"s1": (1,), "s2": (1,), "s3": (1,)})
        db = self._make_df({"s1": (1,), "s2": (1,)})
        to_s, from_s, affected = detect_singleton_transitions(df, db, "20 SNPs")
        assert "s3" not in to_s and "s3" not in from_s


# ---------------------------------------------------------------------------
# _inject_split_if_affected
# ---------------------------------------------------------------------------

def _make_match_array(*cluster_tuples):
    """Build the 1-D object numpy array that match_db_cluster actually returns."""
    arr = np.empty(len(cluster_tuples), dtype=object)
    for i, t in enumerate(cluster_tuples):
        arr[i] = t
    return arr


class TestInjectSplitIfAffected:
    def _row(self, split: bool, match):
        return pd.Series({"split": split, "20 SNPs_db_match": match})

    def test_already_split(self):
        row = self._row(True, None)
        assert _inject_split_if_affected(row, "20 SNPs_db_match", {(1,)}) is True

    def test_affected_cluster(self):
        row = self._row(False, _make_match_array((1,)))
        assert _inject_split_if_affected(row, "20 SNPs_db_match", {(1,)}) is True

    def test_not_affected(self):
        row = self._row(False, _make_match_array((2,)))
        assert _inject_split_if_affected(row, "20 SNPs_db_match", {(1,)}) is False

    def test_no_match_not_affected(self):
        row = self._row(False, None)
        assert _inject_split_if_affected(row, "20 SNPs_db_match", {(1,)}) is False


# ---------------------------------------------------------------------------
# assign_cluster_database — naming rules
# ---------------------------------------------------------------------------

def _tuple_df(data: dict, cols: list) -> pd.DataFrame:
    """Helper: build a cluster assignment DataFrame with tuple values."""
    df = pd.DataFrame(data).T
    df.columns = cols
    df.index.name = "sample"
    return df


class TestAssignClusterDatabaseNoHistory:
    def test_no_history_returns_data_unchanged(self):
        df = _tuple_df(
            {"s1": [(1,)], "s2": [(1,)]},
            ["20 SNPs"],
        )
        result, new_db = assign_cluster_database(df, None)
        pd.testing.assert_frame_equal(result, df)

    def test_no_history_db_equals_data(self):
        df = _tuple_df({"s1": [(1,)], "s2": [(2,)]}, ["20 SNPs"])
        result, new_db = assign_cluster_database(df, None)
        pd.testing.assert_frame_equal(new_db, df)


class TestAssignClusterDatabaseNewThreshold:
    def test_new_threshold_returns_data(self):
        df = _tuple_df({"s1": [(1,), (1, 1)], "s2": [(1,), (1, 1)]}, ["20 SNPs", "10 SNPs"])
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, ["20 SNPs"])
        result, new_db = assign_cluster_database(df, history)
        pd.testing.assert_frame_equal(result, df)


class TestAssignClusterDatabaseConsistent:
    def test_consistent_keeps_name(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        assert result.loc["s1", "20 SNPs"] == (1,)
        assert result.loc["s2", "20 SNPs"] == (1,)

    def test_consistent_new_member_joins(self):
        # s3 is new but joins existing cluster — history cluster is still consistent
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)], "s3": [(1,)]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        assert result.loc["s1", "20 SNPs"] == (1,)
        assert result.loc["s3", "20 SNPs"] == (1,)


class TestAssignClusterDatabaseSplit:
    def test_split_generates_new_names(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)], "s3": [(2,)]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)], "s3": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        # Both resulting clusters must be NEW names (not (1,))
        name_s1 = result.loc["s1", "20 SNPs"]
        name_s3 = result.loc["s3", "20 SNPs"]
        assert name_s1 != (1,)
        assert name_s3 != (1,)
        assert name_s1 != name_s3

    def test_split_names_are_tuples(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(2,)]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        assert isinstance(result.loc["s1", "20 SNPs"], tuple)
        assert isinstance(result.loc["s2", "20 SNPs"], tuple)


class TestAssignClusterDatabaseMerge:
    def test_merge_generates_new_name(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(2,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        # Merged cluster must get a new name (not 1 or 2)
        name = result.loc["s1", "20 SNPs"]
        assert name not in {(1,), (2,)}
        assert result.loc["s1", "20 SNPs"] == result.loc["s2", "20 SNPs"]


class TestAssignClusterDatabaseNew:
    def test_entirely_new_cluster_gets_name(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)], "s3": [(2,)]}, cols)
        # s3 is completely new, not in history
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        # s1,s2 consistent; s3 gets a new name
        assert result.loc["s1", "20 SNPs"] == (1,)
        new_name = result.loc["s3", "20 SNPs"]
        assert isinstance(new_name, tuple)
        assert new_name != (1,)

    def test_singleton_stays_none(self):
        cols = ["20 SNPs"]
        df = _tuple_df({"s1": [(1,)], "s2": [(1,)], "s3": [None]}, cols)
        history = _tuple_df({"s1": [(1,)], "s2": [(1,)]}, cols)
        result, _ = assign_cluster_database(df, history)
        assert result.loc["s3", "20 SNPs"] is None


# ---------------------------------------------------------------------------
# pretty_name_save
# ---------------------------------------------------------------------------

class TestPrettyNameSave:
    def test_writes_dot_notation(self, tmp_path):
        df = pd.DataFrame(
            {"20 SNPs": {"s1": (1,), "s2": (1,)}, "10 SNPs": {"s1": (1, 1), "s2": (1, 2)}},
            dtype=object,
        )
        df.index.name = "sample"
        out = str(tmp_path / "out.csv")
        pretty_name_save(df, "matrix.tsv", out)
        saved = pd.read_csv(out, index_col=0, dtype=str)
        assert saved.loc["s1", "20 SNPs"] == "1"
        assert saved.loc["s2", "10 SNPs"] == "1.2"

    def test_singleton_written_as_nan(self, tmp_path):
        # Singletons (None) are written as NaN, not omitted from the output
        df = pd.DataFrame({"20 SNPs": {"s1": (1,), "s2": None}}, dtype=object)
        df.index.name = "sample"
        out = str(tmp_path / "out.csv")
        pretty_name_save(df, "matrix.tsv", out)
        saved = pd.read_csv(out, index_col=0, dtype=str)
        assert saved.loc["s1", "20 SNPs"] == "1"
        assert saved.loc["s2", "20 SNPs"] != saved.loc["s2", "20 SNPs"]  # NaN != NaN


# ---------------------------------------------------------------------------
# cluster_from_distance_matrix (integration)
# ---------------------------------------------------------------------------

def _write_snpdist(path: str, samples: list, matrix: list) -> None:
    """Write a snp-dists-style TSV."""
    header = "\t".join([""] + samples)
    rows = [header]
    for i, s in enumerate(samples):
        row = "\t".join([s] + [str(matrix[i][j]) for j in range(len(samples))])
        rows.append(row)
    with open(path, "w") as f:
        f.write("\n".join(rows) + "\n")


class TestClusterFromDistanceMatrix:
    def test_basic_clustering_with_reference(self, tmp_path):
        samples = ["Reference", "s1", "s2", "s3"]
        # s1,s2 within 10 SNPs; s3 far
        matrix = [
            [0, 5, 8, 30],
            [5, 0, 6, 32],
            [8, 6, 0, 28],
            [30, 32, 28, 0],
        ]
        tsv = str(tmp_path / "dist.tsv")
        _write_snpdist(tsv, samples, matrix)
        result = cluster_from_distance_matrix(tsv, [20], ref_exist=True)
        # Reference should be dropped
        assert "Reference" not in result.index
        # s1 and s2 should be in the same cluster (within 20 SNPs)
        assert result.loc["s1", "20 SNPs"] == result.loc["s2", "20 SNPs"]
        # s3 should be singleton (far from others at 20 SNP threshold)
        assert result.loc["s3", "20 SNPs"] is None

    def test_clustering_without_reference(self, tmp_path):
        samples = ["s1", "s2", "s3", "s4"]
        matrix = [
            [0, 5, 25, 30],
            [5, 0, 22, 28],
            [25, 22, 0, 4],
            [30, 28, 4, 0],
        ]
        tsv = str(tmp_path / "dist.tsv")
        _write_snpdist(tsv, samples, matrix)
        result = cluster_from_distance_matrix(tsv, [20], ref_exist=False)
        assert set(result.index) == {"s1", "s2", "s3", "s4"}
        assert result.loc["s1", "20 SNPs"] == result.loc["s2", "20 SNPs"]
        assert result.loc["s3", "20 SNPs"] == result.loc["s4", "20 SNPs"]

    def test_missing_reference_exits(self, tmp_path):
        samples = ["s1", "s2"]
        matrix = [[0, 5], [5, 0]]
        tsv = str(tmp_path / "dist.tsv")
        _write_snpdist(tsv, samples, matrix)
        with pytest.raises(SystemExit):
            cluster_from_distance_matrix(tsv, [20], ref_exist=True)

    def test_returns_dataframe_with_correct_columns(self, tmp_path):
        samples = ["Reference", "s1", "s2"]
        matrix = [[0, 5, 8], [5, 0, 6], [8, 6, 0]]
        tsv = str(tmp_path / "dist.tsv")
        _write_snpdist(tsv, samples, matrix)
        result = cluster_from_distance_matrix(tsv, [20, 10], ref_exist=True)
        assert list(result.columns) == ["20 SNPs", "10 SNPs"]

    def test_threshold_sorted_descending(self, tmp_path):
        samples = ["Reference", "s1", "s2"]
        matrix = [[0, 5, 8], [5, 0, 6], [8, 6, 0]]
        tsv = str(tmp_path / "dist.tsv")
        _write_snpdist(tsv, samples, matrix)
        # pass thresholds in wrong order — should still produce correct column order
        result = cluster_from_distance_matrix(tsv, [5, 20, 10], ref_exist=True)
        assert list(result.columns) == ["20 SNPs", "10 SNPs", "5 SNPs"]
