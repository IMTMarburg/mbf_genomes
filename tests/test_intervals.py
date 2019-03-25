import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from mbf_genomes import intervals


class TestIntervals:
    def test_merge_intervals(self):
        df = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1000, "key": "a"},
                {"chr": "1", "start": 850, "stop": 860, "key": "b"},
                {"chr": "1", "start": 900, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )
        merged = intervals.merge_intervals(df)
        should = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )
        assert_frame_equal(merged, should)

    def test_merge_intervals_with_callback(self):
        df = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1000, "key": "a"},
                {"chr": "1", "start": 850, "stop": 860, "key": "b"},
                {"chr": "1", "start": 900, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )

        def cb(sub_df):
            res = sub_df.iloc[0].to_dict()
            res["key"] = "".join(sub_df["key"])
            return res

        merged = intervals.merge_intervals_with_callback(df, cb)

        should = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1100, "key": "abc"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )
        assert_frame_equal(merged, should)

    def test_merge_intervals_with_callback_raises_on_non_dict(self):
        df = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1000, "key": "a"},
                {"chr": "1", "start": 850, "stop": 860, "key": "b"},
                {"chr": "1", "start": 900, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )

        def cb(sub_df):
            res = sub_df.iloc[0]
            res["key"] = "".join(sub_df["key"])
            return res

        with pytest.raises(ValueError):
            intervals.merge_intervals_with_callback(df, cb)

    def test_merge_intervals_with_callback_raises_on_column_removal(self):
        df = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1000, "key": "a"},
                {"chr": "1", "start": 850, "stop": 860, "key": "b"},
                {"chr": "1", "start": 900, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )

        def cb(sub_df):
            res = sub_df.iloc[0].to_dict()
            del res["key"]
            return res

        with pytest.raises(ValueError):
            intervals.merge_intervals_with_callback(df, cb)

    def test_merge_intervals_with_callback_raises_on_column_add(self):
        df = pd.DataFrame(
            [
                {"chr": "1", "start": 0, "stop": 1000, "key": "a"},
                {"chr": "1", "start": 850, "stop": 860, "key": "b"},
                {"chr": "1", "start": 900, "stop": 1100, "key": "c"},
                {"chr": "2", "start": 900, "stop": 1100, "key": "d"},
            ]
        )

        def cb(sub_df):
            res = sub_df.iloc[0].to_dict()
            res["key_merged"] = "".join(sub_df["key"])
            return res

        with pytest.raises(ValueError):
            intervals.merge_intervals_with_callback(df, cb)
