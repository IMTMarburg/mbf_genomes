from pathlib import Path
from abc import ABC


def open_file(fileNameOrHandle, mode="rb"):
    """Transparently open compressed or uncompressed files"""
    if hasattr(fileNameOrHandle, "read"):
        return fileNameOrHandle
    elif isinstance(fileNameOrHandle, Path):
        fileNameOrHandle = str(fileNameOrHandle)
    if fileNameOrHandle.endswith(".gz"):
        import gzip

        return gzip.GzipFile(fileNameOrHandle, mode)
    elif fileNameOrHandle.endswith(".bz2"):
        import bz2

        return bz2.BZ2File(fileNameOrHandle, mode)
    else:
        return open(fileNameOrHandle, mode)


def chunkify(handle, separator, block_size=None):
    """take a file handle and split it at separator, reading in efficently in 50 mb blocks or so"""
    if block_size is None:
        block_size = 50 * 1024 * 1024
    chunk = handle.read(block_size)
    chunk = chunk.split(separator)
    while True:
        for k in chunk[:-1]:
            yield k
        next = handle.read(block_size)
        if next:
            chunk = chunk[-1] + next
            chunk = chunk.split(separator)
        else:
            yield chunk[-1]
            break


def iter_fasta(filenameOrFileLikeObject, keyFunc=None, block_size=None):
    """An iterator over a fasta file (raw or gzipped).
    Yields tupples of key, sequence  (bytes!) on each iteration
    """
    o = open_file(filenameOrFileLikeObject)
    key = ""
    for chunk in chunkify(o, b"\n>", block_size=block_size):
        if chunk:
            key = chunk[: chunk.find(b"\n")].strip()
            if key.startswith(b">"):
                key = key[1:]
            if keyFunc:
                key = keyFunc(key)
            if chunk.find(b"\n") != -1:
                seq = (
                    chunk[chunk.find(b"\n") + 1 :]
                    .replace(b"\r", b"")
                    .replace(b"\n", b"")
                )
            else:
                raise ValueError("Should not be reached")  # pragma: no cover
                # seq = b""
            yield (key, seq)
    return


def wrappedIterator(width):
    def inner(text):
        i = 0
        length = len(text)
        while i < length:
            yield text[i : i + width]
            i += width

    return inner


rc_table = str.maketrans("agctAGCT", "tcgaTCGA")


def reverse_complement(s):
    """return complementary, reversed sequence to x (keeping case)"""
    return s.translate(rc_table)[::-1]


universal_genenetic_code = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}


class GeneticCode(ABC):
    def translate_dna(self, sequence, raise_on_non_multiple_of_three=True):
        if raise_on_non_multiple_of_three and len(sequence) % 3 != 0:
            raise ValueError("len(sequence) was not a multiple of 3")
        genetic_code = self.genetic_code
        proteinseq = ""
        sequence = sequence.upper()
        if sequence[:3] in self.start_codons:
            proteinseq += "M"
        else:
            proteinseq += genetic_code[sequence[:3]]
        for n in range(3, len(sequence), 3):
            proteinseq += genetic_code[sequence[n : n + 3]]
        return proteinseq

    def translate_dna_till_stop(self, sequence, genetic_code=None):
        genetic_code = self.genetic_code
        proteinseq = ""
        sequence = sequence.upper()
        sequence = sequence.upper()
        if sequence[:3] in self.start_codons:
            proteinseq += "M"
        else:
            proteinseq += genetic_code[sequence[:3]]
        for n in range(3, len(sequence), 3):
            try:
                x = genetic_code[sequence[n : n + 3]]
                proteinseq += x
                if x == "*":
                    break
            except KeyError:
                if n + 3 < len(sequence):
                    raise
                else:
                    break
        return proteinseq


class EukaryoticCode(GeneticCode):
    """Genetic code for eukaryotes"""

    genetic_code = universal_genenetic_code
    start_codons = ["ATG"]


class ProkaryoticCode(GeneticCode):
    """Genetic code for prokaryotes - e.g. E. coli"""

    genetic_code = universal_genenetic_code
    # for e coli, from wikipedia) - 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG[7] and one or two others (e.g., an AUU and possibly a CUG
    start_codons = ["ATG", "GTG", "TTG", " ATT", "CTG"]
