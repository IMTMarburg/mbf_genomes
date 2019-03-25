import pytest
import numpy as np
from pathlib import Path
from mbf_genomes import HardCodedGenome
import pandas as pd


data_path = Path(__file__).parent.absolute() / "sample_data"

default_chr_lengths = {
    "1": 100_000,
    "2": 200_000,
    "3": 300_000,
    "4": 400_000,
    "5": 500_000,
}


def DummyGenome(df_genes, df_transcripts=None):

    df_genes = df_genes.rename(columns={"stable_id": "gene_stable_id"})
    if not "start" in df_genes.columns:
        starts = []
        stops = []
        for idx, row in df_genes.iterrows():
            if row["strand"] == 1:
                starts.append(row["tss"])
                stops.append(row["tes"])
            else:
                starts.append(row["tes"])
                stops.append(row["tss"])
        df_genes = df_genes.assign(start=starts, stop=stops)
    if not "biotype" in df_genes.columns:
        df_genes = df_genes.assign(biotype="protein_coding")
    df_genes = df_genes.sort_values(["chr", "start"])
    df_genes = df_genes.set_index("gene_stable_id")
    if df_transcripts is not None:
        if not "biotype" in df_transcripts.columns:
            df_transcripts = df_transcripts.assign(biotype="protein_coding")
        if "exons" in df_transcripts.columns:
            if len(df_transcripts["exons"].iloc[0]) == 3:
                df_transcripts = df_transcripts.assign(
                    exons=[(x[0], x[1]) for x in df_transcripts["exons"]]
                )
        df_transcripts = df_transcripts.set_index("transcript_stable_id")
    return HardCodedGenome("dummy", default_chr_lengths, df_genes, df_transcripts, None)


def test_transcript_get_introns():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4900, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3750, 4000)],
                    [(4900, 5000), (5100, 5400)],
                    [(4900, 5000), (5100, 5200), (5222, 5400)],
                ],
            }
        ),
    )
    assert genome.transcript("trans1a").introns == [(3000, 3100)]
    assert genome.transcript("trans1b").introns == [(3500, 3750), (4000, 4900)]
    assert genome.transcript("trans2").introns == [(5000, 5100)]
    assert genome.transcript("trans3").introns == [(5000, 5100), (5200, 5222)]


def test_gene_get_introns_merging():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4900, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4800)],
                    [(3000, 3500), (3750, 4000)],
                    [(4900, 5000), (5100, 5400)],
                    [(4900, 5000), (5100, 5200), (5222, 5400)],
                ],
            }
        ),
    )
    print(genome.gene("fake1").introns)
    assert genome.gene("fake1").introns == [(4800, 4900)]
    # assert genome.transcript("trans1b").introns == [(3500, 3750), (4000, 4900)]
    assert genome.gene("fake2").introns == [(5000, 5100)]
    assert genome.gene("fake3").introns == [(5000, 5100), (5200, 5222)]


def test_intronify_more_complex():
    transcript = {
        "chr": "2R",
        "exons": [
            (14_243_005, 14_244_766),
            (14_177_040, 14_177_355),
            (14_176_065, 14_176_232),
            (14_175_632, 14_175_893),
            (14_172_742, 14_175_244),
            (14_172_109, 14_172_226),
            (14_170_836, 14_172_015),
            (14_169_750, 14_170_749),
            (14_169_470, 14_169_683),
            (14_169_134, 14_169_402),
            (14_167_751, 14_169_018),
            (14_166_570, 14_167_681),
        ],
        "gene_stable_id": "FBgn0010575",
        "start": 14_166_570,
        "stop": 14_244_766,
        "strand": -1,
        "transcript_stable_id": "FBtr0301547",
    }
    gene = {
        "biotype": "protein_coding",
        "chr": "2R",
        "description": "CG5580 [Source:FlyBase;GeneId:FBgn0010575]",
        "name": "sbb",
        "stable_id": "FBgn0010575",
        "strand": -1,
        "tes": 14_166_570,
        "tss": 14_244_766,
    }
    genome = DummyGenome(pd.DataFrame([gene]), pd.DataFrame([transcript]))
    g = genome.transcript("FBtr0301547")
    assert g.gene_id == "FBgn0010575"
    introns = g.introns
    assert (
        np.array(introns)
        == [
            (14_167_681, 14_167_751),
            (14_169_018, 14_169_134),
            (14_169_402, 14_169_470),
            (14_169_683, 14_169_750),
            (14_170_749, 14_170_836),
            (14_172_015, 14_172_109),
            (14_172_226, 14_172_742),
            (14_175_244, 14_175_632),
            (14_175_893, 14_176_065),
            (14_176_232, 14_177_040),
            (14_177_355, 14_243_005),
        ]
    ).all()
    assert genome.get_chromosome_lengths() == default_chr_lengths
    assert genome.job_genes() is None


def test_intron_intervals_raises_on_inverted():
    transcript = {
        "chr": "2R",
        "exons": [
            (14_243_005, 14_244_766),
            (14_177_040, 14_177_355),
            (14_176_065, 14_176_232),
            (14_175_632, 14_175_893),
            (14_172_742, 14_175_244),
            (14_172_109, 14_172_226),
            (14_172_015, 14_170_836),  # inverted
            (14_169_750, 14_170_749),
            (14_169_470, 14_169_683),
            (14_169_134, 14_169_402),
            (14_167_751, 14_169_018),
            (14_166_570, 14_167_681),
        ],
        "gene_stable_id": "FBgn0010575",
        "start": 14_166_570,
        "stop": 14_244_766,
        "strand": -1,
        "transcript_stable_id": "FBtr0301547",
    }
    gene = {
        "biotype": "protein_coding",
        "chr": "2R",
        "description": "CG5580 [Source:FlyBase;GeneId:FBgn0010575]",
        "name": "sbb",
        "stable_id": "FBgn0010575",
        "strand": -1,
        "tes": 14_166_570,
        "tss": 14_244_766,
    }
    genome = DummyGenome(pd.DataFrame([gene]), pd.DataFrame([transcript]))
    g = genome.transcript("FBtr0301547")
    with pytest.raises(ValueError):
        g.introns


def test_get_gene_introns():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5500,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    one = genome.gene("fake1").introns
    assert len(one) == 0

    two = genome.gene("fake2").introns
    assert two == [(4900, 4910), (5000, 5100), (5400, 5500)]


def test_get_gene_exons():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    g = genome.gene("fake2")
    two = g.exons_merged
    assert (two["start"] == [4910, 5100]).all()
    assert (two["stop"] == [5000, 5400]).all()
    assert "chr" in two.columns
    assert "strand" in two.columns
    four = genome.gene("fake4").exons_merged
    assert len(four) == 0


def test_get_gene_exons_protein_coding():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                    "biotype": "lincRNA",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "biotype": ["protein_coding", "whatever", "protein_coding", "lincRNA"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [3200, 4000, 5400, 5400],
                "exons": [
                    [(3100, 3200)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    g = genome.gene("fake1")
    one = g.exons_protein_coding_merged
    print(one)
    assert (one["start"] == [3100]).all()
    assert (one["stop"] == [3200]).all()

    g = genome.gene("fake2")
    two = g.exons_protein_coding_merged
    assert (two["start"] == [4910, 5100]).all()
    assert (two["stop"] == [5000, 5400]).all()

    g = genome.gene("fake3")
    three = g.exons_protein_coding_merged
    assert (three["start"] == [4900]).all()
    assert (three["stop"] == [5400]).all()


def test_gene_exons_overlapping():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    one = genome.gene("fake1").exons_overlapping
    assert (one["start"] == [3000, 3100, 3300, 3750]).all()
    assert (one["stop"] == [3500, 4900, 3330, 4000]).all()
    assert "chr" in one.columns
    assert "strand" in one.columns
    assert len(genome.gene("fake4").exons_overlapping) == 0
