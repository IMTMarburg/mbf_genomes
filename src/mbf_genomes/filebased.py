import pypipegraph as ppg
from pathlib import Path
from .base import GenomeBase, include_in_downloads, class_with_downloads
from .common import reverse_complement, iter_fasta, wrappedIterator, EukaryoticCode
from mbf_externals.prebuild import PrebuildFileInvariantsExploding


@class_with_downloads
class FileBasedGenome(GenomeBase):
    def __init__(
        self,
        name,
        genome_fasta_file,
        gtf_file=None,
        cdna_fasta_file=None,
        protein_fasta_file=None,
        genetic_code=EukaryoticCode,
        cache_dir=None,
    ):
        """
        Parameters
        ----------
            name: str
            genome_fasta_file: job / Path / str, or list of such
            gtf_file: job / Path / str
            cdna_fasta_file: job / path / str / list of such, None
                generated from gtf + genome if not set
            protein_fasta_file: job / path / str / list of such, None
                generated from gtf + genome if not set
            genetic_code: a common.GeneticCode object defining the translation tables and start codons
        """
        super().__init__()
        self.name = name
        ppg.assert_uniqueness_of_object(self)
        self.genetic_code = genetic_code
        if cache_dir is None:
            self.cache_dir = (
                Path(ppg.util.global_pipegraph.cache_folder)
                / "FileBasedGenome"
                / self.name
            )
        else:  # pragma: no cover
            self.cache_dir = Path(cache_dir)
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

        self.protein_fasta_filename = self.cache_dir / "protein" / "pep.fasta"
        if protein_fasta_file:
            self.protein_fasta_dependencies = self.prep_fasta(
                protein_fasta_file, self.protein_fasta_filename
            )
        elif gtf_file:
            self.protein_fasta_dependencies = self.create_protein_from_genome_and_gtf()
        else:
            self.protein_fasta_filename = None
            self.protein_fasta_dependencies = []

        self._filename_lookups = {
            "genome.fasta": self.genome_fasta_filename,
            "cdna.fasta": self.cdna_fasta_filename,
            "genes.gtf": self.gtf_filename,
            "proteins.fasta": self.protein_fasta_filename,
        }

    def prep_fasta(self, input_filenames, output_filename):
        if isinstance(input_filenames, ppg.Job):
            filenames = input_filenames.filenames
            deps = input_filenames
        elif isinstance(input_filenames, (str, Path)):
            filenames = [str(input_filenames)]
            deps = PrebuildFileInvariantsExploding(
                str(output_filename) + "_prep_fasta", filenames
            )
        else:
            filenames = [str(x) for x in input_filenames]
            deps = PrebuildFileInvariantsExploding(
                str(output_filename) + "_prep_fasta", filenames
            )

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
            self.job_transcripts(), self.genome_fasta_dependencies
        )
        self._download_jobs.append(job)
        return job

    def create_protein_from_genome_and_gtf(self):
        Path(self.protein_fasta_filename).parent.mkdir(exist_ok=True)

        def create(output_filename):
            proteins = self.df_proteins
            with open(output_filename, "w") as op:
                for protein_stable_id, protein_info in proteins.iterrows():
                    cdna = self.get_cds_sequence(protein_stable_id, protein_info)
                    seq = self.genetic_code.translate_dna(
                        cdna, raise_on_non_multiple_of_three=False
                    )
                    seq = "".join(wrappedIterator(80)(seq))
                    op.write(f">{protein_stable_id}\n{seq}\n")

        job = ppg.FileGeneratingJob(self.protein_fasta_filename, create).depends_on(
            self.job_proteins(), self.genome_fasta_dependencies
        )
        self._download_jobs.append(job)
        return job

    @include_in_downloads
    def prepare_gtf(self):
        return self.gtf_dependencies

    def _msg_pack_job(self, property_name, filename, callback_function):
        out_dir = self.cache_dir / "lookup"
        out_dir.mkdir(exist_ok=True)

        def dump(output_filename):
            df = callback_function(self)
            df.to_msgpack(output_filename)

        j = ppg.FileGeneratingJob(out_dir / filename, dump).depends_on(
            ppg.FunctionInvariant(out_dir / filename / property_name, callback_function)
        )
        self._prebuilds.append(j)
        return j


@class_with_downloads
class InteractiveFileBasedGenome(GenomeBase):
    def __init__(
        self,
        name,
        genome_fasta_filename,
        cdna_fasta_filename,
        protein_fasta_filename,
        gtf_filename,
        cache_dir,
    ):
        """
        A FileBasedGenome used for interactive work,
        uses files that a FileBasedGenome has created in a previous ppg run.
        """
        super().__init__()
        self.name = name
        self.cache_dir = Path(cache_dir)

        self.genome_fasta_filename = genome_fasta_filename
        self.cdna_fasta_filename = cdna_fasta_filename
        self.protein_fasta_filename = protein_fasta_filename
        self.gtf_filename = gtf_filename

        self._filename_lookups = {
            "genome.fasta": self.genome_fasta_filename,
            "cdna.fasta": self.cdna_fasta_filename,
            "protein.fasta": self.protein_fasta_filename,
            "genes.gtf": self.gtf_filename,
            "df_genes.msgpack": self.cache_dir / "lookup" / "df_genes.msgpack",
            "df_transcripts.msgpack": self.cache_dir
            / "lookup"
            / "df_transcripts.msgpack",
        }

        if ppg.util.inside_ppg():
            self.gtf_dependencies = ppg.FileInvariant(self.gtf_filename)
        else:
            self.gtf_dependencies = []

    def _msg_pack_job(self, property_name, filename, callback_function):
        out_dir = self.cache_dir / "lookup"
        out_dir.mkdir(exist_ok=True)

        if not ppg.util.inside_ppg():
            if not Path(filename).exists():  # pragma: no branch
                df = callback_function(self)
                df.to_msgpack(out_dir / filename)
        else:

            def dump(output_filename):
                df = callback_function(self)
                df.to_msgpack(output_filename)

            j = ppg.FileGeneratingJob(out_dir / filename, dump).depends_on(
                ppg.FunctionInvariant(
                    out_dir / filename / property_name, callback_function
                )
            )
            self._prebuilds.append(j)
            return j

        return
