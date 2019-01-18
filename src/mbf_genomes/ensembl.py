import re
from pathlib import Path
import mbf_externals
from mbf_externals.util import download_file_and_gunzip, lazy_property, get_page
from .base import GenomeBase, include_in_downloads, class_with_downloads
import pypipegraph as ppg


@class_with_downloads
class EnsemblGenome(GenomeBase):
    def __init__(self, species, revision, prebuild_manager=None):
        super().__init__()
        if prebuild_manager is None:  # pragma: no cover
            prebuild_manager = mbf_externals.prebuild.global_manager
        self.prebuild_manager = prebuild_manager

        self.species = species
        if not re.match(r"^[A-Z][a-z]+_[a-z]+$", species):
            raise ValueError("Species must be capitalized like 'Homo_sapiens")
        self.revision = str(int(revision))
        self.name = f"{self.species}_{self.revision}"
        ppg.util.assert_uniqueness_of_object(self)

    @include_in_downloads
    def _pb_find_server(self):
        ensembl_urls = [
            "http://ftp.ensembl.org/pub/release-%i/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/fungi/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/metazoa/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/plants/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/protists/",
            # "http://ftp.ensemblgenomes.org/pub/release-%i/bacteria/", # bacteria are complicated / subdivided?
        ]

        def find_ensembl_server(output_path):
            for proto_url in ensembl_urls:
                url = proto_url % (int(self.revision),) + "fasta/"
                r = get_page(url)
                if self.species.lower() in r:
                    (output_path / "url.txt").write_text(
                        proto_url % (int(self.revision),)
                    )
                    return
            raise ValueError("Could not find this species on any ensembl server")

        server_job = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/server",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "1",
            [],
            ["url.txt"],
            find_ensembl_server,
            minimum_acceptable_version="1",
            maximum_acceptable_version="1",
        )
        self.server_job = server_job
        return server_job

    @include_in_downloads
    def _pb_download_gtf(self):
        return self._pb_download(
            "gtf",
            "gtf/" + self.species.lower() + "/",
            fr"{self.species}\..+\.{self.revision}.gtf.gz",
            "genes.gtf",
        )

    @include_in_downloads
    def _pb_download_genome_fasta(self):
        return self._pb_download(
            "dna",
            "fasta/" + self.species.lower() + "/dna/",
            fr"{self.species}\..+\.dna.toplevel.fa.gz",
            "genome.fasta",
        )

    @include_in_downloads
    def _pb_download_cdna_fasta(self):
        return self._pb_download(
            "cdna",
            "fasta/" + self.species.lower() + "/cdna/",
            fr"{self.species}\..+\.cdna.all.fa.gz",
            "cdna.fasta",
        )

    def _pb_download(self, pb_name, url, regexps, output_filename):
        def do_download(output_path):
            real_url = self.base_url + url
            raw = get_page(real_url)
            if not raw:  # pragma: no cover
                raise ValueError("Retrieving url failed: %s" % real_url)
            matches = re.findall(regexps, raw)
            if len(matches) != 1:  # pragma: no cover
                raise ValueError(
                    "Found wrong number of matches for regexps: %s. \nRaw was %s"
                    % (matches, raw)
                )
            download_file_and_gunzip(
                real_url + matches[0], output_path / output_filename
            )
            if Path(output_filename).suffix == ".fasta":
                import pysam

                pysam.faidx(str((output_path / output_filename).absolute()))

        job = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/{pb_name}",
            "1",
            [],
            [output_filename],
            do_download,
        )
        job.depends_on(self.server_job)
        return job

    @lazy_property
    def base_url(self):
        return self.server_job.find_file("url.txt").read_text()

    def job_genes(self):
        def dump(output_filename):
            df = self._prepare_df_genes()
            df.to_msgpack(output_filename / "df_genes.msgpack")

        j = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/genes_msgpack",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "1",
            [],
            ["df_genes.msgpack"],
            dump,
        )
        j.depends_on(
            ppg.FunctionInvariant(
                Path("EnsemblGenome") / self.name / "df_genes.msgpack" / "func",
                self.__class__._prepare_df_genes,
            ),
            self._pb_download_gtf()
        )
        self._prebuilds.append(j)
        return j

    def job_transcripts(self):
        def dump(output_filename):
            df = self._prepare_df_transcripts()
            df.to_msgpack(output_filename / "df_transcripts.msgpack")

        j = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/transcripts_msgpack",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "1",
            [],
            ["df_transcripts.msgpack"],
            dump,
        )
        j.depends_on(
            ppg.FunctionInvariant(
                Path("EnsemblGenome") / self.name / "df_transcripts.msgpack" / "func",
                self.__class__._prepare_df_transcripts,
            ),
            self._pb_download_gtf()
        )
        self._prebuilds.append(j)
        return j
