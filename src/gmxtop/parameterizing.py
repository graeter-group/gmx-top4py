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

"""Base classes and basic instances for parameterizing topologies."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from gmxtop.topology.topology import Topology


class Parameterizer(ABC):
    def __init__(self, **kwargs):
        self.type_scheme = dict()

    @abstractmethod
    def parameterize_topology(
        self, current_topology: Topology, focus_nrs: Optional[set[str]]
    ) -> Topology:
        pass


class BasicParameterizer(Parameterizer):
    """reconstruct base force field state"""

    def parameterize_topology(
        self, current_topology: Topology, focus_nrs: Optional[set[str]] = None
    ) -> Topology:
        """Do nothing,
        all necessary actions should already have happened in bind_bond and break_bond of Topology
        """
        return current_topology
