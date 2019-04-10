from . import FileBasedGenome
from mbf_sampledata import get_sample_path


def get_Candidatus_carsonella_ruddii_pv(name=None, **kwargs):
    """A FilebasedGenome used by other libraries for their tests"""
    if name is None:  # pragma: no cover
        name = "Candidatus_carsonella"
    return FileBasedGenome(
        name,
        get_sample_path("mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"),
        get_sample_path("mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"),
        get_sample_path("mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"),
        get_sample_path("mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"),
        **kwargs,
    )
