import pandas as pd
from .intervals import merge_intervals


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


class Gene:
    def __init__(self, genome, gene_stable_id):
        self.genome = genome
        self.gene_stable_id = gene_stable_id

    @property
    def data(self):
        return self.genome.df_genes.loc[self.gene_stable_id]

    @property
    def transcript_ids(self):
        return self.genome.df_transcripts[
            self.genome.df_transcripts.gene_stable_id == self.gene_stable_id
        ].index

    @property
    def introns(self):
        """Get truly intronic regions - ie. not covered by any exon for this gene
        result is  [(start, stop),...]

        """
        gene = self.data
        gene_start = min(gene["tss"], gene["tes"])
        gene_stop = max(gene["tss"], gene["tes"])
        introns = {"start": [], "stop": []}
        exons = []
        for transcript_stable_id in self.transcript_ids:
            tr = self.genome.transcript(transcript_stable_id)
            exons.extend(tr.exons)

        exons.sort()
        transcript_introns = _intron_intervals_from_exons(
            exons, gene_start, gene_stop, True
        )
        for start, stop in transcript_introns:
            introns["start"].append(start)
            introns["stop"].append(stop)
        introns["chr"] = tr.data["chr"]
        introns = pd.DataFrame(introns)
        introns = merge_intervals(introns)
        return list(zip(introns["start"], introns["stop"]))

    @property
    def exons_merged(self):
        """Get the merged exon regions for a gene given by gene_stable_id
        result is a DataFrame{chr, start, stop}

        """
        exons = {"start": [], "stop": []}
        for transcript_stable_id in self.transcript_ids:
            tr = self.genome.transcript(transcript_stable_id)
            for exon_start, exon_stop in tr.exons:
                exons["start"].append(exon_start)
                exons["stop"].append(exon_stop)
        if exons['start']:
            exons["chr"] = tr.data["chr"]
            exons["strand"] = tr.data["strand"]
        else:
            exons['chr'] = []
            exons['strand'] = []
        exons = pd.DataFrame(exons)
        exons = merge_intervals(exons)
        return exons

    @property
    def exons_overlapping(self):
        """Get the overlapping exon regions for a gene given by gene_stable_id
        result is a DataFrame{chr, start, stop}

        """
        exons = {"start": [], "stop": []}
        for transcript_stable_id in self.transcript_ids:
            tr = self.genome.transcript(transcript_stable_id)
            for exon_start, exon_stop in tr.exons:
                exons["start"].append(exon_start)
                exons["stop"].append(exon_stop)
        if exons['start']:
            exons["chr"] = tr.data["chr"]
            exons["strand"] = tr.data["strand"]
        else:
            exons['chr'] = []
            exons['strand'] = []
        exons = pd.DataFrame(exons)
        if len(exons) > 1:
            # exons = merge_intervals(exons)
            exons = exons.sort_values(["start", "stop"])
        return exons

    @property
    def exons_protein_coding_merged(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a DataFrame{chr, start, stop}
        """
        res = self.exons_protein_coding_overlapping
        return merge_intervals(res)

    @property
    def exons_protein_coding_overlapping(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a DataFrame{chr, start, stop}
        """
        # New rule: if a gene has a single protein coding transcript, we only take protein coding transcripts
        # earlier on, we required the gene to be annotated as protein_coding.
        # but that's not strictly correct
        # for example polymorphismic_pseudogenes can have protein coding variants.
        biotypes = self.genome.df_transcripts[
            self.genome.df_transcripts.gene_stable_id == self.gene_stable_id
        ]["biotype"]
        if not (biotypes == "protein_coding").any():
            return self.exons_overlapping
        exons = {"start": [], "stop": []}
        for transcript_stable_id in self.transcript_ids:
            tr = self.genome.transcript(transcript_stable_id)
            if tr.data["biotype"] == "protein_coding":
                for exon_start, exon_stop in tr.exons:
                    exons["start"].append(exon_start)
                    exons["stop"].append(exon_stop)
        exons["chr"] = tr.data["chr"]
        exons["strand"] = tr.data["strand"]
        exons = pd.DataFrame(exons)
        # if len(exons) > 1:
        # exons = merge_intervals(exons)
        return exons


class Transcript:
    def __init__(self, genome, transcript_stable_id):
        self.genome = genome
        self.transcript_stable_id = transcript_stable_id

    @property
    def data(self):
        return self.genome.df_transcripts.loc[self.transcript_stable_id]

    @property
    def gene_id(self):
        return self.data["gene_stable_id"]

    @property
    def exons(self):
        return [(start, stop) for (start, stop) in self.data["exons"]]

    @property
    def introns(self):
        """Return [(start, stop),...] for all introns in the transcript
        Order is in genomic order.
        Intron is defined as everything inside tss..tes that is not an exon,
        so if a gene, by any reason would extend beyond it's exons,
        that region would also be covered.
        """
        print(self.genome.df_transcripts)
        transcript_info = self.genome.df_transcripts.loc[self.transcript_stable_id]
        gene_info = self.genome.df_genes.loc[transcript_info["gene_stable_id"]]
        gene_start = min(gene_info["tss"], gene_info["tes"])
        gene_stop = max(gene_info["tss"], gene_info["tes"])
        exons = sorted(transcript_info["exons"])
        return _intron_intervals_from_exons(exons, gene_start, gene_stop)
