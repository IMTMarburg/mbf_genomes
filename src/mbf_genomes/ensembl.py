import re
from pathlib import Path
import mbf_externals
from mbf_externals.util import (
    download_file_and_gunzip,
    lazy_property,
    get_page,
    lazy_method,
)
from .base import GenomeBase, include_in_downloads, class_with_downloads
import pypipegraph as ppg
from .common import EukaryoticCode


@class_with_downloads
class EnsemblGenome(GenomeBase):
    def __init__(self, species, revision, prebuild_manager=None):
        super().__init__()
        if prebuild_manager is None:  # pragma: no cover
            prebuild_manager = mbf_externals.get_global_manager()
        self.prebuild_manager = prebuild_manager

        self.species = species
        if not re.match(r"^[A-Z][a-z]+_[a-z]+$", species):
            raise ValueError("Species must be capitalized like 'Homo_sapiens")
        self.revision = str(int(revision))
        self.name = f"{self.species}_{self.revision}"
        if ppg.inside_ppg():
            ppg.util.assert_uniqueness_of_object(self)
        self.genetic_code = EukaryoticCode
        self.download_genome()

    @include_in_downloads
    def _pb_find_server(self):
        ensembl_urls = [
            "ftp://ftp.ensembl.org/pub/release-%i/",
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
            (fr"{self.species}\..+\.{self.revision}.gtf.gz",),
            "genes.gtf",
        )  # can't have this unziped star wants it unziped

    def get_additional_gene_gtfs(self):
        if self.species == "Homo_sapiens":
            if int(self.revision) <= 74:
                return [
                    Path(__file__).parent.parent.parent
                    / "data"
                    / "ribosomal_genes_grch37.gtf.gz.full.gtf.gz"
                ]
            else:
                return [
                    Path(__file__).parent.parent.parent
                    / "data"
                    / "ribosomal_genes_grch38.gtf.gz.full.gtf.gz"
                ]
        elif self.species == "Mus_musculus":
            if int(self.revision) > 67:
                return [
                    Path(__file__).parent.parent.parent
                    / "data"
                    / "ribosomal_genes_mm10.gtf.gz.full.gtf.gz"
                ]
        return []

    @property
    def gene_gtf_dependencies(self):
        return [self._pb_download_gtf()]

    @include_in_downloads
    def _pb_download_genome_fasta(self):
        return self._pb_download(
            "dna",
            "fasta/" + self.species.lower() + "/dna/",
            (
                fr"{self.species}\..+\.dna.primary_assembly.fa.gz",
                fr"{self.species}\..+\.dna.toplevel.fa.gz",
            ),
            "genome.fasta",
        )

    @include_in_downloads
    def _pb_extract_keys_from_genome(self):
        output_filename = "references.txt"

        def extract(output_path):
            from .common import iter_fasta

            fn = self.find_file("genome.fasta")
            keys = []
            for key, seq in iter_fasta(fn):
                keys.append(key)
            (output_path / output_filename).write_bytes(b"\n".join(keys))

        job = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/chromosomes_and_contigs",
            "1",
            [],
            [output_filename],
            extract,
        ).depends_on(self._pb_download_genome_fasta())
        return job

    @include_in_downloads
    def _pb_download_cdna_fasta(self):
        return self._pb_download(
            "cdna",
            "fasta/" + self.species.lower() + "/cdna/",
            (fr"{self.species}\..+\.cdna.all.fa.gz",),
            "cdna.fasta",
        )

    @include_in_downloads
    def _pb_download_proptein_fasta(self):
        return self._pb_download(
            "pep",
            f"fasta/{self.species.lower()}/pep/",
            (fr"{self.species}\..+\.pep.all.fa.gz",),
            "pep.fasta",
        )

    def _pb_download(self, pb_name, url, regexps, output_filename):
        """regexps may be multiple - then the first one matching is used"""

        def do_download(output_path):
            real_url = self.base_url + url
            raw = get_page(real_url)
            if not raw:  # pragma: no cover
                raise ValueError("Retrieving url failed: %s" % real_url)
            for aregexps in regexps:
                matches = re.findall(aregexps, raw)
                if len(matches) == 1:
                    Path(output_filename + ".url").write_text((real_url + matches[0]))
                    download_file_and_gunzip(
                        real_url + matches[0], output_path / output_filename
                    )
                    break
            else:
                raise ValueError(  # pragma: no cover - defensive
                    "Found either too few or too many for every regexps. \nRaw was %s"
                    % (raw,)
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
        job.depends_on(self._pb_find_server())
        return job

    @lazy_property
    def base_url(self):
        return self.server_job.find_file("url.txt").read_text()

    def _msg_pack_job(self, property_name, filename, callback_function):
        def dump(output_filename):
            df = callback_function(self)
            df.to_msgpack(output_filename / filename)

        j = self.prebuild_manager.prebuild(
            f"ensembl/{self.species}_{self.revision}/{property_name}",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "1",
            [],
            [filename],
            dump,
        )
        j.depends_on_func(property_name, callback_function)
        self._prebuilds.append(j)
        return j

    @lazy_method
    def get_true_chromosomes(self):
        """Get the names of 'true' chromosomes, ie. no scaffolds/contigs
        in genomes that have chromosomes, otherwise all"""
        fn = self.find_file("references.txt")
        keys = Path(fn).read_text().split("\n")
        chroms = [x for x in keys if "chromosome:" in x]
        if not chroms:
            chroms = keys
        return [x[: x.find(" ")] for x in chroms]
