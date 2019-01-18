from pathlib import Path
import gtfparse
from dppd import dppd
import pandas as pd
import pysam
from mbf_externals.util import lazy_property

dp, X = dppd()


def normalize_strand(strand):
    return strand.replace({"+": 1, "-": -1, ".": 0})


def include_in_downloads(func):
    """A decorator to collect the download funcs"""
    func._include_in_downloads = True
    return func


def class_with_downloads(cls):
    cls._download_methods = []
    for f in cls.__dict__.items():
        if hasattr(f[1], "_include_in_downloads"):
            cls._download_methods.append(f[1])
    return cls


class GenomeBase:
    def __init__(self):
        self._filename_lookups = []
        self._prebuilds = []
        self._download_jobs = []

    def download_genome(self):
        """All the jobs needed to download the genome and prepare it for usage"""
        result = []
        for method in self.__class__._download_methods:
            j = method(self)
            if isinstance(j, list):
                if j:
                    result.extend(j)
            else:
                result.append(j)
            # for j in result:
            # if isinstance(j, list):
            # raise ValueError(method)
        for j in self._download_jobs:
            result.append(j)
        for j in result:
            if not j in self._prebuilds:
                self._prebuilds.append(j)
        return result

    def find_file(self, name):
        if name in self._filename_lookups:
            return self._filename_lookups[name]
        for job in self._prebuilds:
            if hasattr(job, "find_file"):
                try:
                    return job.find_file(name)
                except KeyError:
                    pass
            else:
                for f in job.filenames:
                    if Path(f).name == name:
                        return Path(f)
        # now search for undeclared, but created files
        # mostly for aligners, where we only track the sentinels
        for job in self._prebuilds:
            if hasattr(job, "name_file"):
                if job.name_file(name).exists():
                    return job.name_file(name)
        raise OSError(f"File not found: {name}")

    def find_prebuild(self, name):
        """Find which prebuild created the file named @name.
        Must be in the list of job.filenames"""

        for job in self._prebuilds:
            if hasattr(job, "find_file"):
                try:
                    job.find_file(name)
                    return job
                except KeyError:
                    pass
            else:
                for f in job.filenames:
                    if Path(f).name == name:
                        return job
        raise OSError(f"File not found: {name}")

    def build_index(self, aligner, fasta_to_use=None, gtf_to_use=None):
        if fasta_to_use is None:
            fasta_to_use = "genome.fasta"
        if gtf_to_use is None:
            gtf_to_use = "genes.gtf"
        name = Path(fasta_to_use).stem

        def do_align(output_path):
            aligner.build_index(
                [self.find_file(fasta_to_use)],
                self.find_file(gtf_to_use) if gtf_to_use is not None else None,
                output_path,
            )

        min_ver, max_ver = aligner.get_index_version_range()

        job = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/indices/{name}/{aligner.name}",
            aligner.version,
            [],
            ["sentinel.txt", "stdout.txt", "stderr.txt", "cmd.txt"],
            do_align,
            minimum_acceptable_version=min_ver,
            maximum_acceptable_version=max_ver,
        )
        job.depends_on(self.find_prebuild(fasta_to_use))
        job.depends_on(self.find_prebuild(gtf_to_use))
        return job

    def get_chromosome_lengths(self):
        f = pysam.FastaFile(str(self.find_file("genome.fasta")))
        return dict(zip(f.references, f.lengths))

    def get_genome_sequence(self, chr, start, stop):
        f = pysam.FastaFile(str(self.find_file("genome.fasta")))
        return f.fetch(chr, start, stop)

    def get_cdna_sequence(self, transcript_stable_id):
        with pysam.FastaFile(str(self.find_file("cdna.fasta"))) as f:
            return f.fetch(transcript_stable_id)

    def get_gtf(self):
        gtf_filename = self.find_file("genes.gtf")
        if gtf_filename is None:
            return pd.DataFrame({})
        return gtfparse.read_gtf(gtf_filename)

    def _prepare_df_genes(self):
        df = self.get_gtf()
        if len(df) == 0:  # a genome without gene information
            return pd.DataFrame(
                {
                    "gene_stable_id": [],
                    "name": [],
                    "chr": [],
                    "start": [],
                    "stop": [],
                    "strand": [],
                    "tss": [],
                    "tes": [],
                    "biotype": [],
                }
            )
        genes = df[df.feature == "gene"]
        transcripts = (
            df[df.feature == "transcript"]
            .set_index("gene_id")
            .sort_values(["start", "end"])
        )
        genes = (
            dp(genes)
            .transassign(
                gene_stable_id=X.gene_id,
                name=X.gene_name,
                chr=pd.Categorical(X.seqname),
                start=X.start - 1,
                stop=X.end,
                strand=normalize_strand(X.strand),
                tss=(genes.start - 1).where(genes.strand == "+", genes.end),
                tes=(genes.end).where(genes.strand == "+", genes.start - 1),
                biotype=pd.Categorical(X.gene_biotype),
            )
            .set_index("gene_stable_id")
            .pd
        )
        tr = {}
        for gene_stable_id, transcript_stable_id in transcripts[
            "transcript_id"
        ].items():
            if not gene_stable_id in tr:
                tr[gene_stable_id] = []
            tr[gene_stable_id].append(transcript_stable_id)
        genes = genes.assign(
            transcript_stable_ids=pd.Series(list(tr.values()), index=list(tr.keys()))
        )
        return genes

    @lazy_property
    def df_genes(self):
        fn = self.find_file("df_genes.msgpack")
        return pd.read_msgpack(fn)

    def _prepare_df_transcripts(self):
        df = self.get_gtf()
        if len(df) == 0:
            df = pd.DataFrame(
                {
                    "transcript_stable_id": [],
                    "gene_stable_id": [],
                    "name": [],
                    "chr": [],
                    "start": [],
                    "stop": [],
                    "strand": [],
                    "biotype": [],
                    "exons": [],
                    "exon_stable_ids": [],
                    "translation_start": [],
                    "translation_start_exon": [],
                    "translation_stop": [],
                    "translation_stop_exon": [],
                    "protein_id": [],
                }
            ).set_index("transcript_stable_id")
            return df
        transcripts = df[df.feature == "transcript"]
        all_exons = (
            df[df.feature == "exon"].set_index("transcript_id").sort_values("start")
        )

        result = (
            dp(transcripts)
            .transassign(
                transcript_stable_id=X.transcript_id,
                gene_stable_id=X.gene_id,
                name=X.transcript_name,
                chr=pd.Categorical(X.seqname),
                start=X.start - 1,
                stop=X.end,
                strand=normalize_strand(X.strand),
                biotype=pd.Categorical(X.transcript_biotype),
            )
            .set_index("transcript_stable_id")
            .pd
        )
        result_exons = {}
        result_exon_ids = {}
        for transcript_stable_id, exon in all_exons.iterrows():
            estart = exon["start"] - 1
            estop = exon["end"]
            if not transcript_stable_id in result_exons:
                result_exons[transcript_stable_id] = []
                result_exon_ids[transcript_stable_id] = []
            result_exons[transcript_stable_id].append((estart, estop))
            result_exon_ids[transcript_stable_id].append(exon["exon_id"])
        result_exons = pd.Series(
            list(result_exons.values()), index=list(result_exons.keys())
        )
        result_exon_ids = pd.Series(
            list(result_exon_ids.values()), index=list(result_exon_ids.keys())
        )
        result = result.assign(exons=result_exons, exon_stable_ids=result_exon_ids)
        return result

    @lazy_property
    def df_transcripts(self):
        """Get a DataFrame with all the transcript information
            transcript_stable_id (index),
            gene_stable_id,
            name,
            chr, start, stop, strand,
            biotype,
            exons  - list (start, stop, phase)
        """
        fn = self.find_file("df_transcripts.msgpack")
        if not fn.exists():
            raise ValueError(
                "df_transcripts called before job_transcripts() job has run"
            )
        return pd.read_msgpack(fn)
