import pypipegraph as ppg
from pathlib import Path
from .base import GenomeBase, include_in_downloads, class_with_downloads
from .common import reverse_complement, iter_fasta, wrappedIterator
from mbf_externals.prebuild import PrebuildFileInvariantsExploding


@class_with_downloads
class FileBasedGenome(GenomeBase):
    def __init__(self, name, genome_fasta_file, gtf_file=None, cdna_fasta_file=None):
        """
        Parameters
        ----------
            name: str
            genome_fasta_file: job / Path / str, or list of such
            gtf_file: job / Path / str
            cdna_fasta_file: job / path / str / list of such, None
                generated from gtf + genome if not set
        """
        super().__init__()
        self.name = name
        ppg.assert_uniqueness_of_object(self)
        self.cache_dir = (
            Path(ppg.util.global_pipegraph.cache_folder) / "FileBasedGenome" / self.name
        )
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.genome_fasta_filename = self.cache_dir / "dna" / "genome.fasta"
        self.genome_fasta_dependencies = self.prep_fasta(
            genome_fasta_file, self.genome_fasta_filename
        )
        self.gtf_filename, self.gtf_dependencies = ppg.util.job_or_filename(
            gtf_file,
            lambda x: PrebuildFileInvariantsExploding(self.name + "_gtf_file", [x]),
        )
        self.cdna_fasta_filename = self.cache_dir / "cdna" / "cdna.fasta"
        if cdna_fasta_file:
            self.cdna_fasta_dependencies = self.prep_fasta(
                cdna_fasta_file, self.cdna_fasta_filename
            )
        elif gtf_file:
            self.cdna_fasta_dependencies = self.create_cdna_from_genome_and_gtf()
        else:
            self.cdna_fasta_filename = None
            self.cdna_fasta_dependencies = []

        self._filename_lookups = {
            "genome.fasta": self.genome_fasta_filename,
            "cdna.fasta": self.cdna_fasta_filename,
            "genes.gtf": self.gtf_filename,
        }

    def prep_fasta(self, input_filenames, output_filename):
        if isinstance(input_filenames, ppg.Job):
            filenames = input_filenames.filenames
            deps = input_filenames
        elif isinstance(input_filenames, (str, Path)):
            filenames = [str(input_filenames)]
            deps = PrebuildFileInvariantsExploding(self.name + "_prep_fasta", filenames)
        else:
            filenames = [str(x) for x in input_filenames]
            deps = PrebuildFileInvariantsExploding(self.name + "_prep_fasta", filenames)

        def prep(output_filename):
            import pysam

            with open(output_filename, "wb") as op:
                for fn in filenames:
                    for key, seq in iter_fasta(
                        fn, lambda x: x[: x.find(b" ")] if b" " in x else x
                    ):
                        op.write(
                            b">%s\n%s\n" % (key, b"\n".join(wrappedIterator(80)(seq)))
                        )
            pysam.faidx(output_filename)

        Path(output_filename).parent.mkdir(exist_ok=True)
        job = ppg.FileGeneratingJob(output_filename, prep)
        job.depends_on(deps)
        self._download_jobs.append(job)
        return job

    def create_cdna_from_genome_and_gtf(self):
        Path(self.cdna_fasta_filename).parent.mkdir(exist_ok=True)

        def create(output_filename):
            with open(output_filename, "w") as op:
                for (
                    transcript_stable_id,
                    transcript_info,
                ) in self.df_transcripts.iterrows():
                    seq = ""
                    for start, stop in transcript_info["exons"]:
                        seq += self.get_genome_sequence(
                            transcript_info["chr"], start, stop
                        )
                    if transcript_info["strand"] == -1:
                        seq = reverse_complement(seq)
                    seq = "".join(wrappedIterator(80)(seq))
                    op.write(f">{transcript_stable_id}\n{seq}\n")

        job = ppg.FileGeneratingJob(self.cdna_fasta_filename, create).depends_on(
            self.job_transcripts()
        )
        self._download_jobs.append(job)
        return job

    @include_in_downloads
    def prepare_gtf(self):
        return self.gtf_dependencies

    def job_genes(self):
        out_dir = self.cache_dir / "lookup"
        out_dir.mkdir(exist_ok=True)

        def dump(output_filename):
            df = self._prepare_df_genes()
            df.to_msgpack(output_filename)

        j = ppg.FileGeneratingJob(out_dir / "df_genes.msgpack", dump).depends_on(
            ppg.FunctionInvariant(
                out_dir / "df_genes.msgpack" / "func", self.__class__._prepare_df_genes
            )
        )
        j.depends_on(self.gtf_dependencies)
        self._prebuilds.append(j)

    def job_transcripts(self):
        out_dir = self.cache_dir / "lookup"
        out_dir.mkdir(exist_ok=True)

        def dump(output_filename):
            df = self._prepare_df_transcripts()
            df.to_msgpack(output_filename)

        j = ppg.FileGeneratingJob(out_dir / "df_transcripts.msgpack", dump).depends_on(
            ppg.FunctionInvariant(
                out_dir / "df_transcripts.msgpack" / "func",
                self.__class__._prepare_df_genes,
            )
        )
        j.depends_on(self.gtf_dependencies)
        self._prebuilds.append(j)
        return j
