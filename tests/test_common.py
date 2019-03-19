import pytest
from pathlib import Path
from mbf_genomes.common import reverse_complement, iter_fasta, EukaryoticCode

data_path = Path(__file__).parent.absolute() / "sample_data"


def test_reverse_complement():
    assert reverse_complement("AGTC") == "GACT"
    assert reverse_complement("agtc") == "gact"
    assert reverse_complement("Agtc") == "gacT"
    assert reverse_complement("Ngtc") == "gacN"


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
