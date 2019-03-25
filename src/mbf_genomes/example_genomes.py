from pathlib import Path
from . import FileBasedGenome

data_path = (
    Path(__file__).parent.parent.parent.absolute() / "tests" / "sample_data"
).absolute()


def get_Candidatus_carsonella_ruddii_pv(name=None, **kwargs):
    """A FilebasedGenome used by other libraries for their tests"""
    if name is None:  # pragma: no cover
        name = "Candidatus_carsonella"
    return FileBasedGenome(
        name,
        data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
        data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
        data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz",
        data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz",
        **kwargs,
    )
