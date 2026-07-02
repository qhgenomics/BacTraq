import pytest
import pandas as pd
import os

from BacTraq.makeHistory import read_cluster_table, row_str_to_tuple, history_generate

# ---------------------------------------------------------------------------
# row_str_to_tuple
# ---------------------------------------------------------------------------

class TestRowStrToTuple:
    def _row(self, values: list, cols: list) -> pd.Series:
        return pd.Series(values, index=cols, name="s1")

    def test_single_level(self):
        row = self._row(["1"], ["20 SNPs"])
        result = row_str_to_tuple(row)
        assert result["20 SNPs"] == (1,)

    def test_two_levels(self):
        row = self._row(["1", "2"], ["20 SNPs", "10 SNPs"])
        result = row_str_to_tuple(row)
        assert result["20 SNPs"] == (1,)
        assert result["10 SNPs"] == (2,)

    def test_dot_notation_two_levels(self):
        row = self._row(["1", "1.2"], ["20 SNPs", "10 SNPs"])
        result = row_str_to_tuple(row)
        assert result["20 SNPs"] == (1,)
        assert result["10 SNPs"] == (1, 2)

    def test_dot_notation_three_levels(self):
        row = self._row(["1", "1.2", "1.2.3"], ["20 SNPs", "10 SNPs", "5 SNPs"])
        result = row_str_to_tuple(row)
        assert result["5 SNPs"] == (1, 2, 3)

    def test_singleton_returns_none(self):
        # Singletons are encoded as "0" after fillna in read_cluster_table
        row = self._row(["0"], ["20 SNPs"])
        result = row_str_to_tuple(row)
        assert result["20 SNPs"] is None

    def test_mixed_singleton_and_named(self):
        row = self._row(["1", "0"], ["20 SNPs", "10 SNPs"])
        result = row_str_to_tuple(row)
        assert result["20 SNPs"] == (1,)
        assert result["10 SNPs"] is None


# ---------------------------------------------------------------------------
# read_cluster_table
# ---------------------------------------------------------------------------

class TestReadClusterTable:
    def _write_csv(self, path: str, content: str) -> None:
        with open(path, "w") as f:
            f.write(content)

    def test_reads_columns_in_sorted_threshold_order(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        self._write_csv(csv, "sample,col1,col2\ns1,1,1\ns2,1,2\n")
        result = read_cluster_table(csv, [20, 10])
        assert list(result.columns) == ["20 SNPs", "10 SNPs"]

    def test_sorted_descending_regardless_of_input_order(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        # passing thresholds in ascending order — columns should still be sorted descending
        self._write_csv(csv, "sample,col1,col2\ns1,1,1\n")
        result = read_cluster_table(csv, [5, 20])
        assert list(result.columns) == ["20 SNPs", "5 SNPs"]

    def test_index_is_sample_column(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        self._write_csv(csv, "sample,col1\ns1,1\ns2,1\n")
        result = read_cluster_table(csv, [20])
        assert "s1" in result.index
        assert "s2" in result.index

    def test_na_filled_with_zero_string(self, tmp_path):
        # Missing values become a string starting with "0" so row_str_to_tuple detects singletons
        csv = str(tmp_path / "clusters.csv")
        self._write_csv(csv, "sample,col1,col2\ns1,1,\ns2,1,2\n")
        result = read_cluster_table(csv, [20, 10])
        assert result.loc["s1", "10 SNPs"].startswith("0")


# ---------------------------------------------------------------------------
# history_generate (integration)
# ---------------------------------------------------------------------------

class TestHistoryGenerate:
    def _write_csv(self, path: str, content: str) -> None:
        with open(path, "w") as f:
            f.write(content)

    def test_generates_parquet_file(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        self._write_csv(csv, "sample,col1,col2\ns1,1,1\ns2,1,2\n")
        out_base = str(tmp_path / "history")
        history_generate(csv, thres_list=[20, 10], filename=out_base)
        assert os.path.exists(f"{out_base}.parquet.gz")

    def test_parquet_contains_tuples(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        self._write_csv(csv, "sample,col1,col2\ns1,1,1\ns2,1,2\n")
        out_base = str(tmp_path / "history")
        history_generate(csv, thres_list=[20, 10], filename=out_base)
        result = pd.read_parquet(f"{out_base}.parquet.gz")
        # parquet stores tuples as lists — check the first element is tuple-like
        val = result.loc["s1", "20 SNPs"]
        # pandas/pyarrow may deserialize tuple as list
        assert tuple(val) == (1,) or val == (1,)

    def test_singleton_encoded_as_none(self, tmp_path):
        csv = str(tmp_path / "clusters.csv")
        # s2 has no sub-cluster (singleton at 10 SNP level)
        self._write_csv(csv, "sample,col1,col2\ns1,1,1\ns2,1,\n")
        out_base = str(tmp_path / "history")
        history_generate(csv, thres_list=[20, 10], filename=out_base)
        result = pd.read_parquet(f"{out_base}.parquet.gz")
        assert result.loc["s2", "10 SNPs"] is None
