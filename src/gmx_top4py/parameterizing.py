"""Base classes and basic instances for parameterizing topologies."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from gmx_top4py.topology.topology import Topology


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
