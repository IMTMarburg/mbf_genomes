from pathlib import Path
from abc import ABC, abstractmethod
import pandas as pd
import gtfparse
from dppd import dppd
import pysam
from .common import reverse_complement

dp, X = dppd()


def normalize_strand(strand):
    if not isinstance(strand, pd.Series):
        raise ValueError("normalize_strand expects a pd.series")
    if not strand.isin(["+", "-", "."]).all():
        raise ValueError("normalize_strand called on something that was not +-.")
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


class MsgPackProperty:
    """
        a message pack property is a property x_y that get's
        calculated by a method _prepare_x_y
        and automatically stored/loaded by a caching job
        as msgpack file.
        the actual job used depends on the GenomeBase subclass

        The dependency_callback get's called with the GenomeBase subclass
        instance and can return dependencys for the generated job

        The object has three members afterwards:
            x_y -> get the value returned by _prepare_x_y (lazy load)
            _prepare_x_y -> that's the one you need to implement,
                        it's docstring is copied to this propery
            job_y -> the job that caches _prepare_x_y() results

        """

    def __init__(self, dependency_callback=None):
        self.dependency_callback = dependency_callback


def msgpack_unpacking_class(cls):
    for d in list(cls.__dict__):
        v = cls.__dict__[d]
        if isinstance(v, MsgPackProperty):
            if not "_" in d:
                raise NotImplementedError(
                    "Do not know how to create job name for msg_pack_properties that do not containt  _"
                )
            job_name = "job_" + d[d.find("_") + 1 :]
            filename = d + ".msgpack"
            calc_func = getattr(cls, f"_prepare_{d}")

            def load(self, d=d, filename=filename):
                if not hasattr(self, "__" + d):
                    fn = self.find_file(filename)
                    if not fn.exists():
                        raise ValueError(
                            f"{d} accessed before the respecting {job_name} call"
                        )
                    df = pd.read_msgpack(fn)
                    setattr(self, "__" + d, df)
                return getattr(self, "__" + d)

            p = property(load)
            p.__doc__ == calc_func.__doc__
            setattr(cls, d, p)
            if not hasattr(cls, job_name):

                def gen_job(
                    self,
                    d=d,
                    filename=filename,
                    calc_func=calc_func,
                    dependency_callback=v.dependency_callback,
                ):
                    j = self._msg_pack_job(d, filename, calc_func)
                    j.depends_on(dependency_callback(self))
                    return j

                setattr(cls, job_name, gen_job)
            else:  # pragma: no cover
                pass

    return cls


@msgpack_unpacking_class
class GenomeBase(ABC):
    def __init__(self):
        self._filename_lookups = []
        self._prebuilds = []
        self._download_jobs = []

    @abstractmethod
    def _msg_pack_job(self, property_name, filename, callback_function):
        raise NotImplementedError()

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
            if not j in self._prebuilds:  # pragma: no branch
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
        # mostly for aligners, where we only track the sentinels, not the index
        # files
        for job in self._prebuilds:
            if hasattr(job, "name_file"):  # pragma: no branch
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
        if fasta_to_use is None:  # pragma: no cover
            fasta_to_use = "genome.fasta"
        if gtf_to_use is None:  # pragma: no cover
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
        try:
            job.depends_on(self.find_prebuild(fasta_to_use))
            job.depends_on(self.find_prebuild(gtf_to_use))
        except OSError:
            self.download_genome()  # so that the jobs are there
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

    def get_cds_sequence(self, protein_id, protein_info=None):
        """Get the coding sequence (rna) of a protein"""
        if protein_info is None:
            protein_info = self.df_proteins.loc[protein_id]
        elif protein_info.name != protein_id:
            raise ValueError("protein_id != protein_info['protein_id']")
        cdna = ""
        chr = protein_info["chr"]
        for start, stop in protein_info["cds"]:
            cdna += self.get_genome_sequence(chr, start, stop)
        if protein_info["strand"] not in (1, -1):  # pragma: no cover
            raise ValueError(f'{protein_info["strand"]} was not 1/-1')
        if protein_info["strand"] == -1:
            cdna = reverse_complement(cdna)
        return cdna

    def get_protein_sequence(self, protein_id):
        """Get the AA sequence of a protein"""
        with pysam.FastaFile(str(self.find_file("pep.fasta"))) as f:
            return f.fetch(protein_id)

    def get_gtf(self):
        filenames = [self.find_file("genes.gtf")]
        if hasattr(self, "get_additional_gene_gtfs"):
            filenames.extend(self.get_additional_gene_gtfs())
        dfs = []
        for gtf_filename in filenames:
            if gtf_filename is None:
                dfs.append(pd.DataFrame({}))
            else:
                dfs.append(gtfparse.read_gtf(gtf_filename))
        return pd.concat(dfs)

    def _prepare_df_genes(self):
        """Return a DataFrame with  gene information:
            gene_stable_id
            name
            chr
            start
            stop
            strand
            tss
            tes
            biotype
        """
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

    def _prepare_df_transcripts(self):
        """Get a DataFrame with all the transcript information
            transcript_stable_id (index),
            gene_stable_id,
            name,
            chr, start, stop, strand,
            biotype,
            exons  - list (start, stop, phase)
        """
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

    def _prepare_df_proteins(self):
        """Get a DataFrame with protein information
            protein_stable_id (index)
            transcript_stable_id,
            gene_stable_id,
            chr,
            strand
            cds - [(start, stop)]  # in genomic coordinates
        """
        df = self.get_gtf()
        if len(df) == 0:
            df = pd.DataFrame(
                {
                    "protein_stable_id": [],
                    "transcript_stable_id": [],
                    "gene_stable_id": [],
                    "chr": [],
                    "strand": [],
                    "cds": [],
                }
            ).set_index("protein_stable_id")
            return df
        cds = df[df.feature == "CDS"]
        result = {
            "protein_stable_id": [],
            "transcript_stable_id": [],
            "gene_stable_id": [],
            "chr": [],
            "strand": [],
            "cds": [],
        }

        for protein_stable_id, sub_df in cds.groupby("protein_id"):
            transcript_stable_id = sub_df.transcript_id.iloc[0]
            gene_stable_id = sub_df.gene_id.iloc[0]
            chr = sub_df.seqname.iloc[0]
            strand = sub_df.strand.iloc[0]
            cds = list(zip(sub_df.start - 1, sub_df.end))
            result["protein_stable_id"].append(protein_stable_id)
            result["gene_stable_id"].append(gene_stable_id)
            result["transcript_stable_id"].append(transcript_stable_id)
            result["chr"].append(chr)
            result["strand"].append(strand)
            result["cds"].append(cds)
        result = pd.DataFrame(result).set_index("protein_stable_id")
        result = result.assign(strand=normalize_strand(result.strand))
        return result

    df_genes = MsgPackProperty(lambda self: self.gtf_dependencies)
    df_transcripts = MsgPackProperty(lambda self: self.gtf_dependencies)
    df_proteins = MsgPackProperty(lambda self: self.gtf_dependencies)
