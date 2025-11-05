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

"""
Utilities for shell convenience functions and GROMACS related functions
"""

import logging
import subprocess as sp
from pathlib import Path
from typing import Optional


logger = logging.getLogger(__name__)

TopologyAtomAddress = str
"""Address to an atom in the topology.

Corresponds to the 1-based id in the topology and coordinates file.
Note, that gromacs ids in the atomnr column of the gro file
can overflow due to the fixed width file format.
The line number - 2 (for the title and the number of atoms) is always
the correct the atom id.
"""


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


### GROMACS related functions ###


def get_gmx_dir(
    gromacs_alias: str = "gmx", grompp_prefix: Optional[str] = None
) -> Optional[Path]:
    """Returns the path to the gromacs installation"""

    # get the stder from calling `gmx` to search for the `Data prefix:`
    # line which contains the path to the gromacs installation

    # Add prefix if necesarry
    cmd = [gromacs_alias]
    if grompp_prefix:
        cmd.insert(0, grompp_prefix)
    try:
        r = sp.run(cmd, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        logger.warning("GROMACS not found.")
        return None

    from_source = False
    gmx_prefix = None
    for l in r.stderr.splitlines():
        if l.startswith("Data prefix:"):
            gmx_prefix = Path(l.split()[2])
            if "(source tree)" in l:
                from_source = True
            break

    if gmx_prefix is None:
        logger.warning("GROMACS data directory not found in gromacs message.")
        return None

    if from_source:
        return Path(gmx_prefix) / "share"
    else:
        return Path(gmx_prefix) / "share" / "gromacs"
