import pytest
from pathlib import Path
import pypipegraph as ppg
from mbf_genomes import FileBasedGenome
from mbf_genomes.common import iter_fasta
from mbf_externals.util import UpstreamChangedError


data_path = Path(__file__).parent.absolute() / "sample_data"


@pytest.mark.usefixtures("new_pipegraph")
class TestFilebased:
    def test_indexing(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz",
        )
        g.download_genome()
        g.job_transcripts()
        g.job_genes()
        with pytest.raises(ValueError):
            g.df_transcripts
        with pytest.raises(ValueError):
            g.df_genes
        ppg.run_pipegraph()
        assert g.find_file("genome.fasta").exists()
        assert g.find_prebuild("genome.fasta") == g.genome_fasta_dependencies
        assert g.find_file("genome.fasta").with_suffix(".fasta.fai").exists()
        for should_file, actual_file in [
            (
                data_path
                / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
                g.find_file("genome.fasta"),
            ),
            (
                data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz",
                g.find_file("cdna.fasta"),
            ),
        ]:
            should = dict(iter_fasta(should_file))
            should = {k[: k.find(b" ")]: v for (k, v) in should.items()}
            actual = dict(iter_fasta(actual_file))
            if should != actual:
                assert len(should) == len(actual)
                assert set(should.keys()) == set(actual.keys())
                assert False == should_file  # noqa:E712
        tf = g.df_transcripts
        assert "BAF35033" in tf.index
        assert tf.loc["BAF35033"].exons == ((1313, 2816),)

        gf = g.df_genes
        assert len(gf) == 246
        assert set([len(x) for x in gf.transcript_stable_ids]) == set([1])

    def test_cdna_creation(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
            None,
        )
        g.download_genome()
        g.job_transcripts()
        ppg.run_pipegraph()
        assert g.find_file("genome.fasta").exists()
        assert g.find_file("genome.fasta").with_suffix(".fasta.fai").exists()
        tf = g.df_transcripts
        assert "BAF35033" in tf.index
        assert tf.loc["BAF35033"].exons == ((1313, 2816),)

        should = dict(
            iter_fasta(
                data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            )
        )
        should = {k[: k.find(b" ")]: v for (k, v) in should.items()}
        actual = dict(iter_fasta(g.find_file("cdna.fasta")))
        if actual != should:
            assert not set(should.keys()).difference(
                set(actual.keys())
            )  # they are all here, we just have more (tRNA...)
            for k in should:
                assert actual[k] == should[k]

    def test_empty_gtf_and_cdna(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
            None,
            None,
        )
        g.download_genome()
        assert g.gtf_filename is None
        assert g.cdna_fasta_filename is None
        g.job_transcripts()
        g.job_genes()
        ppg.run_pipegraph()
        assert len(g.df_transcripts) == 0
        assert len(g.get_gtf()) == 0
        assert len(g.df_genes) == 0

    def test_job_creating_fasta(self, new_pipegraph):
        new_pipegraph.quiet = False

        def gen_fasta():
            import shutil

            shutil.copy(
                data_path
                / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
                "shu.fasta.gz",
            )

        fasta_job = ppg.FileGeneratingJob("shu.fasta.gz", gen_fasta)
        g = FileBasedGenome(
            "Candidatus_carsonella",
            fasta_job,
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
            None,
        )
        g.download_genome()
        ppg.run_pipegraph()
        assert (
            g.get_cdna_sequence("BAF35032")
            == "ATGAATACTATATTTTCAAGAATAACACCATTAGGAAATGGTACGTTATGTGTTATAAGAAT"
            "TTCTGGAAAAAATGTAAAATTTTTAATACAAAAAATTGTAAAAAAAAATATAAAAGAAAAAATAG"
            "CTACTTTTTCTAAATTATTTTTAGATAAAGAATGTGTAGATTATGCAATGATTATTTTTTTTAAA"
            "AAACCAAATACGTTCACTGGAGAAGATATAATCGAATTTCATATTCACAATAATGAAACTATTGT"
            "AAAAAAAATAATTAATTATTTATTATTAAATAAAGCAAGATTTGCAAAAGCTGGCGAATTTTTAG"
            "AAAGACGATATTTAAATGGAAAAATTTCTTTAATAGAATGCGAATTAATAAATAATAAAATTTTA"
            "TATGATAATGAAAATATGTTTCAATTAACAAAAAATTCTGAAAAAAAAATATTTTTATGTATAAT"
            "TAAAAATTTAAAATTTAAAATAAATTCTTTAATAATTTGTATTGAAATCGCAAATTTTAATTTTA"
            "GTTTTTTTTTTTTTAATGATTTTTTATTTATAAAATATACATTTAAAAAACTATTAAAACTTTTA"
            "AAAATATTAATTGATAAAATAACTGTTATAAATTATTTAAAAAAGAATTTCACAATAATGATATT"
            "AGGTAGAAGAAATGTAGGAAAGTCTACTTTATTTAATAAAATATGTGCACAATATGACTCGATTG"
            "TAACTAATATTCCTGGTACTACAAAAAATATTATATCAAAAAAAATAAAAATTTTATCTAAAAAA"
            "ATAAAAATGATGGATACAGCAGGATTAAAAATTAGAACTAAAAATTTAATTGAAAAAATTGGAAT"
            "TATTAAAAATATAAATAAAATTTATCAAGGAAATTTAATTTTGTATATGATTGATAAATTTAATA"
            "TTAAAAATATATTTTTTAACATTCCAATAGATTTTATTGATAAAATTAAATTAAATGAATTAATA"
            "ATTTTAGTTAACAAATCAGATATTTTAGGAAAAGAAGAAGGAGTTTTTAAAATAAAAAATATATT"
            "AATAATTTTAATTTCTTCTAAAAATGGAACTTTTATAAAAAATTTAAAATGTTTTATTAATAAAA"
            "TCGTTGATAATAAAGATTTTTCTAAAAATAATTATTCTGATGTTAAAATTCTATTTAATAAATTT"
            "TCTTTTTTTTATAAAGAATTTTCATGTAACTATGATTTAGTGTTATCAAAATTAATTGATTTTCA"
            "AAAAAATATATTTAAATTAACAGGAAATTTTACTAATAAAAAAATAATAAATTCTTGTTTTAGAA"
            "ATTTTTGTATTGGTAAATGA"
        )

    def test_multiple_fasta_files(self, new_pipegraph):
        import tempfile

        tf = tempfile.NamedTemporaryFile(suffix=".fasta")
        tf.write(b">Extra\nAGTC")
        tf.flush()
        g = FileBasedGenome(
            "Candidatus_carsonella",
            [
                data_path
                / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
                tf.name,
            ],
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
            None,
        )
        g.download_genome()
        ppg.run_pipegraph()
        assert g.get_genome_sequence("Extra", 0, 4) == "AGTC"
        assert g.get_chromosome_lengths() == {"Extra": 4, "Chromosome": 159662}

        # test that changing the fasta leads to an explosion
        new_pipegraph.new_pipegraph()
        tf.seek(0, 0)
        tf.write(b">Extra\nAGTCA")
        tf.flush()
        g = FileBasedGenome(
            "Candidatus_carsonella",
            [
                data_path
                / "Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz",
                tf.name,
            ],
            data_path / "Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz",
            None,
        )
        g.download_genome()
        with pytest.raises(UpstreamChangedError):
            ppg.run_pipegraph()
