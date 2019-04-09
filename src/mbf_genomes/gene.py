import attr
import numpy as np
from mbf_nested_intervals import IntervalSet


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
        result is a a tuple of np arrays, (starts, stops)
        """
        gene_start = self.start
        gene_stop = self.stop
        exons = []
        for tr in self.transcripts:
            exons.extend(tr.exons)
        return IntervalSet.from_tuples(exons).invert(gene_start, gene_stop).to_numpy()

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
        return IntervalSet.from_tuples(self._exons).merge_hull().to_numpy()

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
        return (
            IntervalSet.from_tuples(self._exons_protein_coding).merge_hull().to_numpy()
        )

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
        return IntervalSet.from_tuples(exons).invert(gene_start, gene_stop).to_tuples()
