import pandas as pd
from .intervals import merge_intervals
import attr


def _intron_intervals_from_exons(exons, gene_start, gene_stop, merge=False):
    intron_start = gene_start
    res = []
    exon_no = 0
    if merge:
        # print 'merging', exons
        zip_exons = list(zip(*exons))
        # print zip_exons
        temp_df = pd.DataFrame(
            {"chr": ["X"] * len(exons), "start": zip_exons[0], "stop": zip_exons[1]}
        )
        # print temp_df
        temp_df = merge_intervals(temp_df)
        # print temp_df
        exons = zip(*[temp_df["start"], temp_df["stop"]])
        # print 'exons', exons
    for exon_start, exon_stop in exons:  # from left to right, please.
        if exon_stop < exon_start:
            raise ValueError("inverted exon")
        if intron_start != exon_start:
            if intron_start > exon_start:  # the overlapping exons case.
                raise ValueError(
                    "_intron_intervals_from_exons saw exons "
                    "that need merging by setting merge=True,"
                    " but merge was False - should not happen?"
                )  # pragma: no cover
                # arguably we could just recall with merge = True
                # but it's an upstream caller bug, it should know
                # whether to expect overlapping exons (=genes)
                # or not (transcripts) and this is defensive.

            res.append((intron_start, exon_start))
        exon_no += 1
        intron_start = exon_stop
    if intron_start < gene_stop:
        res.append((intron_start, gene_stop))
    return res


@attr.s(slots=True)
class Gene:
    gene_stable_id = attr.ib()
    name = attr.ib()
    chr = attr.ib()
    start = attr.ib()
    stop = attr.ib()
    strand = attr.ib()
    biotype = attr.ib()
    transcripts = attr.ib(default=None)

    @property
    def tss(self):
        return self.start if self.strand == 1 else self.stop

    @property
    def tes(self):
        return self.start if self.strand != 1 else self.stop

    @property
    def introns(self):
        """Get truly intronic regions - ie. not covered by any exon for this gene
        result is  [(start, stop),...]

        """
        gene_start = self.start
        gene_stop = self.stop
        introns = {"start": [], "stop": []}
        exons = []
        for tr in self.transcripts:
            exons.extend(tr.exons)

        exons.sort()
        transcript_introns = _intron_intervals_from_exons(
            exons, gene_start, gene_stop, True
        )
        for start, stop in transcript_introns:
            introns["start"].append(start)
            introns["stop"].append(stop)
        introns["chr"] = self.chr
        introns = pd.DataFrame(introns)
        introns = merge_intervals(introns)
        return list(zip(introns["start"], introns["stop"]))

    @property
    def _exons(self):
        """Common code to exons_merged and exons_overlapping"""
        exons = {"start": [], "stop": []}
        for tr in self.transcripts:
            for exon_start, exon_stop in tr.exons:
                exons["start"].append(exon_start)
                exons["stop"].append(exon_stop)
        if exons["start"]:
            exons["chr"] = self.chr
            exons["strand"] = self.strand
        else:
            exons["chr"] = []
            exons["strand"] = []
        return pd.DataFrame(exons)

    @property
    def exons_merged(self):
        """Get the merged exon regions for a gene given by gene_stable_id
        result is a DataFrame{chr, start, stop}
        """
        exons = self._exons
        exons = merge_intervals(exons)
        return exons

    @property
    def exons_overlapping(self):
        """Get the overlapping exon regions for a gene given by gene_stable_id
        result is a DataFrame{chr, strand, start, stop}
        """
        exons = self._exons
        if len(exons) > 1:
            # exons = merge_intervals(exons)
            exons = exons.sort_values(["start", "stop"])
        return exons

    @property
    def exons_protein_coding_merged(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a DataFrame{chr, strand, start, stop}
        """
        res = self.exons_protein_coding_overlapping
        return merge_intervals(res)

    @property
    def exons_protein_coding_overlapping(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a DataFrame{chr, strand, start, stop}
        """
        # New rule: if a gene has a single protein coding transcript, we only take protein coding transcripts
        # earlier on, we required the gene to be annotated as protein_coding.
        # but that's not strictly correct
        # for example polymorphismic_pseudogenes can have protein coding variants.
        exons = {"start": [], "stop": []}
        found = False
        for tr in self.transcripts:
            if tr.biotype == "protein_coding":
                found = True
                for exon_start, exon_stop in tr.exons:
                    exons["start"].append(exon_start)
                    exons["stop"].append(exon_stop)
        if not found:
            return self.exons_overlapping
        # no need to handle the empty exon case her, it's done in exon_overlapping
        exons["chr"] = self.chr
        exons["strand"] = self.strand
        exons = pd.DataFrame(exons)
        return exons


@attr.s(slots=True)
class Transcript:
    transcript_stable_id = attr.ib()
    gene_stable_id = attr.ib()
    name = attr.ib()
    chr = attr.ib()
    start = attr.ib()
    stop = attr.ib()
    strand = attr.ib()
    biotype = attr.ib()
    exons = attr.ib()
    exon_stable_ids = attr.ib()
    gene = attr.ib()

    @property
    def exons_tuples(self):
        return [(start, stop) for (start, stop) in self.exons]

    @property
    def introns(self):
        """Return [(start, stop),...] for all introns in the transcript
        Order is in genomic order.
        Intron is defined as everything inside tss..tes that is not an exon,
        so if a gene, by any reason would extend beyond it's exons,
        that region would also be covered.
        """
        gene_start = self.gene.start
        gene_stop = self.gene.stop
        exons = sorted(self.exons_tuples)
        return _intron_intervals_from_exons(exons, gene_start, gene_stop)
