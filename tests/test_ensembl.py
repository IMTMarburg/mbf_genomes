import pytest
from mbf_genomes import EnsemblGenome
from mbf_externals import PrebuildManager
from mbf_externals.aligners.subread import Subread
from pathlib import Path
import pypipegraph as ppg
from pypipegraph.util import checksum_file


@pytest.mark.usefixtures("new_pipegraph")
class TestEnsembl:
    def test_download(self, new_pipegraph, mock_download):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        species = (
            "Ashbya_gossypii"
        )  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome(species, "41", prebuild_manager=pb)

        def shorten_genome_fasta(output_path):
            with open(g.find_file("genome.fasta")) as op:
                head = op.read(1024 * 100)
            (output_path / "test.fasta").write_text(head)

        test_fasta_job = g.prebuild_manager.prebuild(
            f"ensembl/{g.species}_{g.revision}/test_fasta",
            "1",
            [],
            ["test.fasta"],
            shorten_genome_fasta,
        )
        test_fasta_job.depends_on(g.download_genome())
        g._prebuilds.append(test_fasta_job)

        subread = Subread(version="1.6.3")
        index = g.build_index(subread, "test.fasta")
        subread_old = Subread(version="1.4.3-p1")
        index_old = g.build_index(subread_old, "test.fasta")

        new_pipegraph.run()
        # note that these are not the checksums from CHECKSUMS files (those are fore
        # the gziped variants, we keep them ungziped and let the filesystem handle
        # the gzip, since we can't rely on the downstream reading gzip...
        assert (
            checksum_file(g.find_file("genome.fasta"))
            == "584a734589964a654c7c1dc23b0167ab"
        )
        assert (
            checksum_file(g.find_file("cdna.fasta"))
            == "3fc1f19ab829573169cb2488abe39211"
        )
        assert (
            checksum_file(g.find_file("genes.gtf"))
            == "8bdeec9b3db5278668dbff8b34e9d93b"
        )
        assert (
            checksum_file(g.find_file("genes.gtf"))
            == "8bdeec9b3db5278668dbff8b34e9d93b"
        )
        assert (
            checksum_file(g.find_file("pep.fasta"))
            == "9580fd44832d419c38469d657f6e2484"
        )
        with pytest.raises(OSError):
            g.find_file("no such file")
        assert index.name_file("subread_index.reads").exists()
        assert index.name_file("subread_index.files").exists()
        assert index.name_file("subread_index.00.b.array").exists()
        assert index_old.name_file("subread_index.reads").exists()
        assert index_old.name_file("subread_index.files").exists()
        assert index_old.name_file("subread_index.00.b.array").exists()
        assert index.name_file("subread_index.reads") != index_old.name_file(
            "subread_index.reads"
        )
        assert g.find_file("test.fasta.md5sum").exists()
        with pytest.raises(OSError):
            assert g.find_file("test.fasta.md5sum.nosuchfile").exists()
        assert g.find_prebuild("test.fasta") is test_fasta_job
        with pytest.raises(OSError):
            assert g.find_prebuild("test.fasta.md5sum.nosuchfile").exists()
        assert g.find_file("genome.fasta.fai").exists()
        assert g.find_file("cdna.fasta.fai").exists()

        new_pipegraph.new_pipegraph()
        pb = PrebuildManager(p)
        g = EnsemblGenome(species, "41", prebuild_manager=pb)
        g.download_genome()
        test_fasta_job = g.prebuild_manager.prebuild(
            f"ensembl/{g.species}_{g.revision}/test_fasta",
            "1",
            [],
            ["test.fasta"],
            shorten_genome_fasta,
        )
        g._prebuilds.append(test_fasta_job)

        subread_intermediate = Subread(version="1.5.0")
        index_intermediate = g.build_index(subread_intermediate, "test.fasta")
        assert index_intermediate.name_file(
            "subread_index.reads"
        ) == index_old.name_file("subread_index.reads")
        index_genome = g.build_index(subread_intermediate)
        assert "/genome/" in str(index_genome.filenames[0])

        assert g.get_chromosome_lengths() == {
            "IV": 1_467_287,
            "MT": 23564,
            "V": 1_519_140,
            "III": 907_494,
            "II": 870_771,
            "VII": 1_800_949,
            "I": 693_414,
            "VI": 1_836_693,
        }

        assert g.get_genome_sequence("VI", 20, 30) == "ACCGCTGAGA"
        assert (
            g.get_cdna_sequence("EFAGOT00000000349")
            == "GCTCGCGTGGCGTAATGGCAACGCGTCTGACTTCTAATCAGAAGATTGTGGGTTCGACCC"
            "CCACCGTGAGTG"
        )
        assert (
            g.get_protein_sequence("AAS53315")
            == "MFSTRICSLLARPFMVPIVPRFGSALLQKPLNGVVVPQFTRGFKVRTSVKKFCAHCYIVR"
            "RKGRVYVYCKSNNKHKQRQG"
        )

    def test_species_formating(self,):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        species = (
            "ashbya_gossypii"
        )  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        with pytest.raises(ValueError):
            EnsemblGenome(species, "41", prebuild_manager=pb)

    def test_unknown_species_raises(self, mock_download):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        species = (
            "Unknown_unknown"
        )  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        EnsemblGenome(species, "41", prebuild_manager=pb).download_genome()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()

    def test_all_genes(self, mock_download):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        g = EnsemblGenome("Ustilago_maydis", 33, pb)
        g.download_genome()
        g.job_genes()
        ppg.run_pipegraph()

        df = g.df_genes
        assert len(df) == 6910
        assert df.loc["UMAG_12015"].strand == 1
        assert df.loc["UMAG_12015"].tss == 4370
        assert df.loc["UMAG_12015"].tes == 5366
        assert df.loc["UMAG_12015"].transcript_stable_ids == ("KIS66832",)

        assert df.loc["UMAG_00663"].strand == -1
        assert df.loc["UMAG_00663"].tss == 1_947_590
        assert df.loc["UMAG_00663"].tes == 1_945_040
        assert df.loc["UMAG_00663"].transcript_stable_ids == ("KIS72250",)

        assert df.loc["UMAG_03168"].transcript_stable_ids == ("KIS68597", "KIS68596")
        assert df["chr"].dtype.name == "category"
        assert df["biotype"].dtype.name == "category"

    def test_all_transcripts(self, mock_download):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        g = EnsemblGenome("Ustilago_maydis", 33, pb)
        g.download_genome()
        g.job_transcripts()
        ppg.run_pipegraph()

        df = g.df_transcripts
        assert "gene_stable_id" in df.columns
        assert len(df) == 6928
        assert df["chr"].dtype.name == "category"
        assert df["biotype"].dtype.name == "category"
        assert df.loc["KIS71021"].chr == "2"
        assert df.loc["KIS71021"].strand == 1
        assert df.loc["KIS71021"].start == 354_742
        assert df.loc["KIS71021"].stop == 356_690
        assert df.loc["KIS71021"].gene_stable_id == "UMAG_12118"
        assert df.loc["KIS71021"].biotype == "protein_coding"
        assert df.loc["KIS71021"].exons == ((354_742, 354_936), (355_222, 356_690))
        assert df.loc["KIS71021"].exon_stable_ids == ("KIS71021-1", "KIS71021-2")

    def test_all_proteins(self, mock_download):
        p = Path("prebuild")
        p.mkdir()
        pb = PrebuildManager(p)
        g = EnsemblGenome("Ustilago_maydis", 33, pb)
        g.download_genome()
        g.job_proteins()
        ppg.run_pipegraph()

        df = g.df_proteins
        assert df.strand.isin([1, -1]).all()
