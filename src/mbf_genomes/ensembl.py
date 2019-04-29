import re
from pathlib import Path
import pandas as pd
import mbf_externals
from mbf_externals.util import (
    download_file_and_gunzip,
    download_file_and_gzip,
    download_file,
    lazy_property,
    get_page,
    lazy_method,
)
from .base import (
    GenomeBase,
    include_in_downloads,
    class_with_downloads,
    MsgPackProperty,
    msgpack_unpacking_class,
)
import pypipegraph as ppg
from .common import EukaryoticCode


@msgpack_unpacking_class
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
        return self._pb_download_and_gunzip(
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
        return self._pb_download_and_gunzip(
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
        return self._pb_download_and_gunzip(
            "cdna",
            "fasta/" + self.species.lower() + "/cdna/",
            (fr"{self.species}\..+\.cdna.all.fa.gz",),
            "cdna.fasta",
        )

    @include_in_downloads
    def _pb_download_proptein_fasta(self):
        return self._pb_download_and_gunzip(
            "pep",
            f"fasta/{self.species.lower()}/pep/",
            (fr"{self.species}\..+\.pep.all.fa.gz",),
            "pep.fasta",
        )

    @include_in_downloads
    def _pb_download_sql_table_definitions(self):
        return self._pb_download_straight(
            "sql/core/sql_def",
            "mysql/",
            (fr"{self.species.lower()}_core.+",),
            "core.sql.gz",
            lambda match: f"{match.strip()}/{match.strip()}.sql.gz",
        )

    def _pb_download_sql_table(self, table_name):
        """Helper to download sql tables as mysql dumps"""
        return self._pb_download_straight(
            f"sql/core/{table_name}",
            "mysql/",
            (fr"{self.species.lower()}_core.+",),
            f"{table_name}.txt.gz",
            lambda match: f"{match.strip()}/{table_name.strip()}.txt.gz",
        ).depends_on(self._pb_download_sql_table_definitions())

    @include_in_downloads
    def _pb_download_sql_table_gene(self):
        return self._pb_download_sql_table("gene")

    def _pb_download(
        self,
        pb_name,
        url,
        regexps,
        output_filename,
        download_func,
        match_transformer=lambda x: x,
    ):
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
                    download_func(
                        real_url + match_transformer(matches[0]),
                        output_path / output_filename,
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

        return self._pb_download_job(pb_name, output_filename, do_download)

    def _pb_download_straight(
        self, pb_name, url, regexps, output_filename, match_transformer=lambda x: x
    ):
        def df(url, filename):
            with open(filename, "wb") as op:
                download_file(url, op)

        return self._pb_download(
            pb_name, url, regexps, output_filename, df, match_transformer
        )

    def _pb_download_and_gunzip(self, pb_name, url, regexps, output_filename):
        return self._pb_download(
            pb_name, url, regexps, output_filename, download_file_and_gunzip
        )

    def _pb_download_and_gzip(self, pb_name, url, regexps, output_filename):
        return self._pb_download(
            pb_name, url, regexps, output_filename, download_file_and_gzip
        )

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

    def _prepare_df_genes_meta(self):
        """Meta data for genes.
        Currently contains:
            'description'
        """
        columns = self._get_sql_table_column_names("gene")
        if not "stable_id" in columns:
            raise ValueError(
                "No stable_id column found - "
                "old ensembl, split into seperate table, add support code?"
            )
        df = pd.read_csv(
            self.find_file("gene.txt.gz"),
            sep="\t",
            header=None,
            names=columns,
            usecols=["stable_id", "description"],
        )
        res = df.set_index("stable_id")
        res.index.name = "gene_stable_id"
        return res

    df_genes_meta = MsgPackProperty(lambda self: [self._pb_download_sql_table_gene()])

    def _get_sql_table_column_names(self, sql_table_name):
        """Read the sql definition and extract column names"""
        import gzip

        with gzip.GzipFile(self.find_file("core.sql.gz")) as op:
            raw = op.read().decode("utf-8")
            for p in raw.split("\n-- Table structure")[1:]:
                p = p[p.find("CREATE TABLE") :]
                if "  PRIMARY" in p:
                    p = p[: p.find("  PRIMARY")]
                elif "  KEY" in p:
                    p = p[: p.find("  KEY")]
                elif "  UNIQUE" in p:
                    p = p[: p.find("  UNIQUE")]

                else:
                    raise ValueError(p)
                names = re.findall("`([^`]+)`", p)
                table_name, *columns = names
                if table_name == sql_table_name:
                    return columns
        raise KeyError(f"{sql_table_name} not in core.sql.gz")
