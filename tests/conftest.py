#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import sys

# import pytest
from pathlib import Path
from pypipegraph.testing.fixtures import (  # noqa: F401
    new_pipegraph,  # noqa: F401
    pytest_runtest_makereport,  # noqa: F401
)  # noqa: F401
from mbf_externals.testing.fixtures import local_store  # noqa:F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))

local_store_path = root / "tests" / "run" / "local_store"
local_store_path.mkdir(exist_ok=True, parents=True)
local_store = local_store(local_store_path)

from mbf_genomes.testing.fixtures import mock_download, shared_prebuild  # noqa: F401
