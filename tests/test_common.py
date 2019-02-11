import pytest
from pathlib import Path
from mbf_genomes.common import (
    reverse_complement,
    open_file,
    chunkify,
    iter_fasta,
    EukaryoticCode,
)

data_path = Path(__file__).parent.absolute() / "sample_data"


def test_reverse_complement():
    assert reverse_complement("AGTC") == "GACT"
    assert reverse_complement("agtc") == "gact"
    assert reverse_complement("Agtc") == "gacT"
    assert reverse_complement("Ngtc") == "gacN"


def test_open_file():
    import gzip
    import tempfile
    import bz2

    tf = tempfile.TemporaryFile()
    assert open_file(tf) is tf

    tf2 = tempfile.NamedTemporaryFile(suffix=".gz", mode="w")
    g = gzip.GzipFile(tf2.name, "w")
    g.write(b"hello")
    g.close()
    tf2.flush()
    assert open_file(tf2.name).read() == b"hello"

    tf3 = tempfile.NamedTemporaryFile(suffix=".bz2", mode="w")
    b = bz2.BZ2File(tf3.name, "w")
    b.write(b"world")
    b.close()
    tf3.flush()
    assert open_file(tf3.name).read() == b"world"


def test_chunkify():
    import tempfile

    tf = tempfile.TemporaryFile("w+")
    tf.write("hello world")
    tf.flush()
    tf.seek(0, 0)
    c = list(chunkify(tf, " ", 2))
    assert c == ["hello", "world"]


def test_iter_fasta():
    fn = data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
    a = list(iter_fasta(fn))
    assert len(a[0][1]) == 159662
    b = list(iter_fasta(fn, block_size=10))
    assert a == b


def test_translate_raises_on_non_multiple_of_three():
    with pytest.raises(ValueError):
        EukaryoticCode.translate_dna("AA")


def test_translate_non_start():
    assert EukaryoticCode.translate_dna("TTG") == "L"


def test_translate_till_stop():
    assert EukaryoticCode.translate_dna_till_stop("ATGTTGTAATTG") == "ML*"


def test_translate_till_stop_non_start_start():
    assert EukaryoticCode.translate_dna_till_stop("TTGTTGTAATTG") == "LL*"


def test_translate_till_stop_not_three():
    with pytest.raises(ValueError):
        EukaryoticCode.translate_dna_till_stop("TTGTTGTA")
