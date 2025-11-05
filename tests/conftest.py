# This file is part of the gmxtop project.
#
# The gmxtop project is based on or includes code from:
#    kimmdy (https://github.com/graeter-group/kimmdy/tree/main)
#    Copyright (C) graeter-group
#    Licensed under the GNU General Public License v3.0 (GPLv3).
#
# Modifications and additional code:
#    Copyright (C) 2025 graeter-group
#    Licensed under the GNU General Public License v3.0 (GPLv3).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""File for pytest configuration in python and fixture definition.

The name 'conftest.py' is recognized by pytest to execute it before tests.
"""

import pytest
import shutil
import os
from pathlib import Path


# fixtures for setup and teardown ##
@pytest.fixture
def arranged_tmp_path(tmp_path: Path, request: pytest.FixtureRequest):
    """Arrange temporary directory for tests.

    With files for the test and a symlink to forcefield.
    """

    # discover_plugins()

    # if fixture was parameterized, use this for directory with input files
    if hasattr(request, "param"):
        file_dir = Path(__file__).parent / "test_files" / request.param
    # else use stem of requesting file to find directory with input files
    else:
        file_dir = Path(__file__).parent / "test_files" / request.path.stem
    # arrange tmp_path
    shutil.copytree(file_dir, tmp_path, dirs_exist_ok=True)
    assetsdir = Path(__file__).parent / "test_files" / "assets"
    if not (tmp_path / "amber99sb-star-ildnp.ff").exists():
        Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
            assetsdir / "amber99sb-star-ildnp.ff",
            target_is_directory=True,
        )
    # change cwd to tmp_path
    os.chdir(tmp_path.resolve())
    return tmp_path
