import re
from pathlib import Path
import pandas as pd
import numpy as np
import mbf_externals
from mbf_externals.util import (
    download_file_and_gunzip,
    # download_file_and_gzip,
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
        self._seq_region_is_canonical = {}

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
    def _pb_download_protein_fasta(self):
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
        job = self._pb_download_straight(
            f"sql/core/{table_name}",
            "mysql/",
            (fr"{self.species.lower()}_core.+",),
            f"{table_name}.txt.gz",
            lambda match: f"{match.strip()}/{table_name.strip()}.txt.gz",
        ).depends_on(self._pb_download_sql_table_definitions())
        job.table_name = table_name
        return job

    @include_in_downloads
    @lazy_method
    def _pb_download_sql_tables(self):
        tables = [
            ("gene"),  # for description
            ("stable_id_event"),  # for stable_id changes
            ("external_db"),  # for external name lookup
            ("object_xref"),  # for external name lookup
            ("xref"),  # for external name lookup
            ("alt_allele"),  # for finding 'canonical' ids
            ("seq_region"),  # for finding 'canonical' ids
            ("seq_region_attrib"),  # for finding 'canonical' ids
            ("attrib_type"),  # for finding 'canonical' ids
        ]
        return [self._pb_download_sql_table(x) for x in tables]

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
                    Path(str(output_path / output_filename) + ".url").write_text(
                        (real_url + matches[0])
                    )
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

    def _pb_download_straight(
        self,
        pb_name,
        url,
        regexps,
        output_filename,
        match_transformer=lambda x: x,  # pragma: no cover
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

    # def _pb_download_and_gzip(self, pb_name, url, regexps, output_filename):
    # return self._pb_download(
    # pb_name, url, regexps, output_filename, download_file_and_gzip
    # )

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

    def _load_from_sql(
        self, table_name, columns=None, check_for_columns=None, **kwargs
    ):
        table_columns = self._get_sql_table_column_names(table_name)
        df = pd.read_csv(
            self.find_file(f"{table_name}.txt.gz"),
            sep="\t",
            header=None,
            names=table_columns,
            usecols=columns,
            na_values="\\N",
            lineterminator="\n",
            escapechar="\\",
            **kwargs,
        )
        if check_for_columns:
            for c in check_for_columns:
                if not c in table_columns:  # pragma: no cover
                    raise KeyError(c)
        return df

    def _prepare_df_genes_meta(self):
        """Meta data for genes.
        Currently contains:
            'description'
        """
        try:
            df = self._load_from_sql(
                "gene", ["stable_id", "description"], ["stable_id"]
            )
        except KeyError:  # pragma: no cover
            raise ValueError(
                "No stable_id column found - "
                "old ensembl, split into seperate table, add support code?"
            )
        res = df.set_index("stable_id")
        res.index.name = "gene_stable_id"
        return res

    df_genes_meta = MsgPackProperty(
        lambda self: [
            x for x in self._pb_download_sql_tables() if x.table_name == "gene"
        ]
    )

    def _get_sql_table_column_names(self, sql_table_name):
        """Read the sql definition and extract column names"""
        import gzip

        with gzip.GzipFile(self.find_file("core.sql.gz")) as op:
            raw = op.read().decode("utf-8")
            for p in raw.split("\n-- Table structure")[1:]:
                p = p[p.find("CREATE TABLE") :]
                if "  PRIMARY" in p:
                    p = p[: p.find("  PRIMARY")]
                elif "  UNIQUE" in p:
                    p = p[: p.find("  UNIQUE")]
                elif "  KEY" in p:
                    p = p[: p.find("  KEY")]

                else:  # pragma: no cover
                    raise ValueError(p)
                names = re.findall("`([^`]+)`", p)
                table_name, *columns = names
                if table_name == sql_table_name:
                    return columns
        raise KeyError(f"{sql_table_name} not in core.sql.gz")  # pragma: no cover

    def _prepare_lookup_stable_id_events(self):
        """Lookup old_stable_id -> new_stable_id"""
        df = self._load_from_sql("stable_id_event", ["old_stable_id", "new_stable_id"])
        lookup = {}
        olds = [str(x) for x in df["old_stable_id"].values]
        news = [str(x) for x in df["new_stable_id"].values]
        for old in olds:
            lookup[old] = set()
        for old, new in zip(olds, news):
            lookup[old].add(new)
        return pd.DataFrame(
            {"old": list(lookup.keys()), "new": [list(x) for x in lookup.values()]}
        ).set_index("old")

    lookup_stable_id_events = MsgPackProperty(
        lambda self: [
            x
            for x in self._pb_download_sql_tables()
            if x.table_name == "stable_id_event"
        ]
    )

    def newest_stable_ids_for(self, stable_id):
        """Get the most up to date and current stable_ids for genes, transcripts, proteins).
        Plurarl for gene might have split, or have been deleted.
        returns a set of new ids.
 """
        try:
            valid_ids = set(self.df_genes.index)
            valid_ids.update(self.df_transcripts.index)
            valid_ids.update(self.df_proteins.index)
            res = set(self.lookup_stable_id_events.loc[stable_id]["new"])
            res = set(
                [x for x in res if x in valid_ids]
            )  # filter those that are no longer in the database - no matter that they were m             apped somewhere else in between
            return res
        except KeyError as e:
            # see if it's a current id where we're simply lacking the stable_id_event for some reason
            if stable_id in valid_ids:
                return set([stable_id])
            else:
                raise e

    def get_external_dbs(self):
        """Return the names of all external dbs that actually have xrefs"""
        df_external_db = self._load_from_sql(
            "external_db", ["external_db_id", "db_name"]
        )
        with_data = set(
            self._load_from_sql("xref", ["external_db_id"])["external_db_id"].unique()
        )
        return sorted(
            df_external_db["db_name"][df_external_db["external_db_id"].isin(with_data)]
        )

    def get_external_db_to_gene_id_mapping(self, external_db_name):
        """Return a dict external id -> set(stable_id, ...)
        for a given external db - e.g. EntrezGene, list
        with get_external_dbs()
        """
        df_external_db = self._load_from_sql(
            "external_db", ["external_db_id", "db_name"]
        ).set_index("db_name")
        external_db_id = df_external_db.at[external_db_name, "external_db_id"]
        xref = self._load_from_sql(
            "xref", ["dbprimary_acc", "external_db_id", "xref_id"]
        ).set_index("xref_id")
        xref = xref[xref.external_db_id == external_db_id]
        object_xref = self._load_from_sql(
            "object_xref", ["ensembl_object_type", "xref_id", "ensembl_id"]
        )
        object_xref = object_xref[object_xref["ensembl_object_type"] == "Gene"]
        object_xref = object_xref[object_xref.xref_id.isin(set(xref.index))]
        result = {}
        genes = self._load_from_sql("gene", ["gene_id", "stable_id"]).set_index(
            "gene_id"
        )
        for internal_id, xref_id in zip(
            object_xref["ensembl_id"].values, object_xref["xref_id"].values
        ):
            gene_stable_id = genes.at[internal_id, "stable_id"]
            db_primary = xref.at[xref_id, "dbprimary_acc"]
            if not db_primary in result:
                result[db_primary] = set()
            result[db_primary].add(gene_stable_id)
        return result

    @lazy_property
    def allele_groups(self):
        df = self._load_from_sql("alt_allele", ["alt_allele_group_id", "gene_id"])
        gene_df = self._load_from_sql("gene", ["gene_id", "stable_id"]).rename(
            columns={"gene_stable_id": "stable_id"}
        )
        df = df.join(gene_df.set_index("gene_id"), "gene_id")
        return df

    def name_to_canonical_id(self, name):
        """Given a gene name, lookup up it's stable ids, and return the
        one that's on the primary assembly from the allele group"""
        name_candidates = set(
            [x for x in self.name_to_gene_ids(name) if not x.startswith("LRG")]
        )
        if not name_candidates:  # pragma: no cover
            raise KeyError("No gene named %s" % name)
        ag = self.allele_groups
        ag_ids = ag.alt_allele_group_id[ag.stable_id.isin(name_candidates)].unique()
        ag_candidates = set(ag.stable_id[ag.alt_allele_group_id.isin(ag_ids)])
        if len(ag_ids) == 1 and name_candidates.issubset(ag_candidates):
            # the easy case, everything matches
            on_primary = [
                x
                for x in ag_candidates
                if x
                in self.df_genes.index  # for there is no entry in genes.gtf if it's not on a not 'non_ref' chromosome.
            ]
            if len(on_primary) == 1:
                return on_primary[0]
            elif len(on_primary) == 0:
                #if self.species == "Homo_sapiens" and name == "HLA-DRB3": # HLA-DRB3 is not in genes.gtf!
                    # known issue - return basically any of the candidates on alternate regions, but be consistent.
                    #return sorted(ag_candidates)[0]
                raise ValueError("No primary gene found for %s" % name)
            else:
                raise ValueError(
                    "Multiple gene on primary assemblies found for %s" % name
                )
        elif len(ag_ids) == 0 and len(name_candidates) == 1:
            # another easy case, there are no alternatives
            return list(name_candidates)[0]
        else:
            raise ValueError(
                "Could not determine canonical gene for '%s'. use name_to_gene_ids()"
                " and have a look yourself (don't forget the allele groups).\n"
                "Name candidates: %s\n"
                "AG candidates: %s\n"
                "AG ids: %s" % (name, name_candidates, ag_candidates, ag_ids)
            )
