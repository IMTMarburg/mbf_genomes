from pathlib import Path


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
