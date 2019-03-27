import pandas as pd
from .intervals import merge_intervals
import attr
import numpy as np


def merge_exons(exons):
    """take a list of exons as [(start, stop), ...] tuples
    (ie. they already have to be on the same chromosome and strand)
    and return two lists [starts], [stops]
    which have the merged)

    This is much faster than the pandas based variants in intervals I'm afraid.

    Converting this to cython shaves off about half a second (reduction to 87%)
    from merging all exons of all genes of the human genome -> not worth it.
    """
    exons.sort()
    starts = np.array([x[0] for x in exons])
    stops = np.array([x[1] for x in exons])
    ii = 0
    lendf = len(starts)
    keep = np.zeros((lendf,), dtype=np.bool)
    last_stop = 0
    last_row = None
    while ii < lendf:
        if starts[ii] < last_stop:
            starts[ii] = starts[last_row]
            stops[ii] = max(stops[ii], last_stop)
        else:
            if last_row is not None:
                keep[last_row] = True
                # new_rows.append(df.get_row(last_row))
        if stops[ii] > last_stop:
            last_stop = stops[ii]
        last_row = ii
        ii += 1
    if last_row is not None:
        keep[last_row] = True
    return starts[keep], stops[keep]


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
        exons = []
        for tr in self.transcripts:
            exons.extend(tr.exons)
        return exons

    @property
    def exons_merged(self):
        """Get the merged exon regions for a gene given by gene_stable_id
        result is a a tuple of np arrays, (starts, stops)
        """
        return merge_exons(self._exons)

    @property
    def exons_overlapping(self):
        """Get the overlapping exon regions for a gene given by gene_stable_id
        result is a a tuple of np arrays, (starts, stops)
        not sorted
        """
        return self._reformat_exons(self._exons)

    def _reformat_exons(self, exons):
        """Turn exons [(start, stop), ...] into [[start, ...], [stop, ...]
        """
        exons.sort()
        return np.array([x[0] for x in exons]), np.array([x[1] for x in exons])

    @property
    def _exons_protein_coding(self):
        """common code for the exons_protein_coding_* propertys"""
        exons = []
        for tr in self.transcripts:
            if tr.biotype == "protein_coding":
                exons.extend(tr.exons)
        return exons

    @property
    def exons_protein_coding_merged(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a a tuple of np arrays, (starts, stops)
        """
        return merge_exons(self._exons_protein_coding)

    @property
    def exons_protein_coding_overlapping(self):
        """Get the overlapping exon regions for a gene, only for protein coding transcripts.
        Empty result on non protein coding genes

        Result is a DataFrame{chr, strand, start, stop}

        We test biotype on transcripts, not on genes,
        because for example polymorphismic_pseudogenes can have protein coding variants.
        """
        return self._reformat_exons(self._exons_protein_coding)


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
