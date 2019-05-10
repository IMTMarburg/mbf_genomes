from pathlib import Path
from .base import GenomeBase, HardCodedGenome
from .ensembl import EnsemblGenome
from .filebased import FileBasedGenome, InteractiveFileBasedGenome

data_path = Path(__file__).parent.parent.parent / "data"
__version__ = "0.1"

__all__ = [
    "GenomeBase",
    "EnsemblGenome",
    "FileBasedGenome",
    "InteractiveFileBasedGenome",
    "data_path",
    "HardCodedGenome",
]
