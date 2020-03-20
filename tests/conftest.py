#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import sys
import pytest
from pathlib import Path
from pypipegraph.testing.fixtures import (  # noqa: F401
    new_pipegraph,  # noqa: F401
    pytest_runtest_makereport,  # noqa: F401
)  # noqa: F401
from mbf_externals.testing.fixtures import local_store  # noqa:F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))


def _mock_get_page(url):
    import hashlib
    import mbf_externals

    p = (
        Path(__file__).parent
        / ".testing_download_cache"
        / hashlib.md5(url.encode("utf-8")).hexdigest()
    )
    p.parent.mkdir(exist_ok=True)
    if not p.exists():
        p.write_text(mbf_externals.util.get_page(url))
    return p.read_text()


def _mock_download_file_and_gunzip(url, filename):
    import shutil
    import hashlib
    import mbf_externals

    p = (
        Path(__file__).parent
        / ".testing_download_cache"
        / hashlib.md5(url.encode("utf-8")).hexdigest()
    )
    p.parent.mkdir(exist_ok=True)
    if not p.exists():
        mbf_externals.util.download_file_and_gunzip(url, p)
    return shutil.copyfile(p, filename)


@pytest.fixture
def mock_download():
    import mbf_genomes

    org_get_page = mbf_genomes.ensembl.get_page
    org_download_file_and_gunzip = mbf_genomes.ensembl.download_file_and_gunzip
    mbf_genomes.ensembl.get_page = _mock_get_page
    mbf_genomes.ensembl.download_file_and_gunzip = _mock_download_file_and_gunzip
    yield
    mbf_genomes.ensembl.get_page = org_get_page
    mbf_genomes.ensembl.download_file_and_gunzip = org_download_file_and_gunzip


first_shared_prebuild = True


@pytest.fixture()
def shared_prebuild():
    global first_shared_prebuild
    p = Path("../prebuild")
    if first_shared_prebuild:
        if p.exists():
            import shutil

            shutil.rmtree(p)
        p.mkdir()
        first_shared_prebuild = False
    from mbf_externals import PrebuildManager

    return PrebuildManager(p)
