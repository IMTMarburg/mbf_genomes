import pytest
import pandas as pd

# from pathlib import Path
# from mbf_externals import PrebuildManager
# from mbf_genomes import EnsemblGenome
# import pypipegraph as ppg


@pytest.mark.usefixtures("new_pipegraph")
class TestBase:
    pass


def test_normalize_strand_wrong_calls():
    from mbf_genomes.base import normalize_strand

    with pytest.raises(ValueError):
        normalize_strand(pd.Series(["+", "-", "?"]))
    with pytest.raises(ValueError):
        normalize_strand(pd.Series([0, 1, -1]))
    with pytest.raises(ValueError):
        normalize_strand(["+", "-"])
    with pytest.raises(ValueError):
        normalize_strand("+")
    with pytest.raises(ValueError):
        normalize_strand("-")


def test_msgpack_unpacking_class_wrong_property_name():
    from mbf_genomes.base import msgpack_unpacking_class, MsgPackProperty

    with pytest.raises(NotImplementedError):

        @msgpack_unpacking_class
        class Shu:
            def _pepare_dfnolower():
                return pd.DataFrame()

            dfnolower = MsgPackProperty()
