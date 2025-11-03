# %%
# This file is part of the gmx-top4py project.
#
# The gmx-top4py project is based on or includes code from:
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
from pathlib import Path
from copy import deepcopy

from gmx_top4py.topology.topology import Topology
from gmx_top4py.parsing import read_top, write_top
from gmx_top4py.topology.utils import get_is_selected_moleculetype_f
# %% [markdown]
# # 1. Reading a topology file
# Internally, the `Topology.from_path` class method uses the `read_top` function to parse the topology file into a dictionary, and then creates an instance of the `Topology` object from that dictionary. <br>
# So the following two approaches are equivalent:
# %% [markdown] 
# a) Create the topology object directly
# %%
top_path = Path("urea.top")
top = Topology.from_path(top_path)
original_top = deepcopy(top)  # Keep a copy of the original topology for later use
print(top)
# %% [markdown]
# b) Create the topology object from a dictionary
# %%
top_path = Path("urea.top")
top_dict = read_top(top_path) 
top = Topology(top_dict)
print(top)
# %% [markdown]
# ## 1.1. Specifying a force field directory
# Per default, the force field will be loaded from the GROMACS installation path. However, you can specify a different force field directory when reading the topology file.
# Here, we read a topology file that requires the 'amber99sb-star-ildnp.ff' force field, which is not defautly shipped with GROMACS, so we need to specify the path to the force field directory explicitly.
# %% [markdown] 
# a) Create the topology object directly
# %%
top_path = Path("urea_amber99sb_star_ildnp.top")
top = Topology.from_path(top_path, ffdir=Path("amber99sb-star-ildnp.ff")) 
# %% [markdown]
# b) Create the topology object from a dictionary
# %%
top_dict = read_top(top_path, ffdir=Path("amber99sb-star-ildnp.ff"))
top = Topology(top_dict)
# %% [markdown]
# ## 1.2. Selecting specific molecules
# Per default, all molecules whose molecule type is not a solvent or a inorganic ion are selected when reading the topology file. In the case of the 'urea.top' example, the 'Urea' molecule type is selected.
# A selected molecule and its properties (i.e, atoms, bonds, angles, im-/proper_dihedrals, exclusions, settles, pairs, position restraints, dihedral restraints, radicals, and nrexcl) are directly linked to the instance of the topology object.
# %% [markdown]
# The molecule corresponding to the selected molecule type can be dircetly accessed via the 'selected_molecule' attribute of the Topology object:
# %%
print(original_top.selected_molecule)
# %% [markdown]
# All properties of the selected molecule type are also directly accessible via the Topology object. Thus, in the 'urea.top' example, the atoms of the urea molecule can be accessed in three different ways.
# %% [markdown]
# a) Via the `Topology` instance
# %%
print(original_top.atoms)
# %% [markdown]
# b) Via the `selected_molecule` attribute
# %%
print(original_top.selected_molecule.atoms)
# %% [markdown]
# c) Via the `moleculetypes` attribute
# %%
print(original_top.moleculetypes["Urea(1x)"].atoms)
# %% [markdown]
# More importantly, the multiplicity of the selected molecule is made explicit: 
# If there is more than one molecule of the selected molecule type, then the properties (e.g., atoms, bonds, angles, and dihedrals) of the selected molecule are duplicated for each copy of the molecule.
# %% [markdown]
# In case of the 'urea.top' example, there is only one molecule of the 'Urea' molecule type.
# %%
print("Single urea molecule in 'urea.top':")
print(original_top.selected_molecule)
# %% [markdown]
# In contrast, there are two urea molecules in the 'urea-times-2.top' file.
# %%
print("Two urea molecules in 'urea-times-2.top':")
top_path = Path("urea-times-2.top")
top = Topology.from_path(top_path)
print(top.selected_molecule)
# %% [markdown]
# If multiple different molecule types are selected, then these molecule types are merged into a single one.
# %%
top_path = Path("naa50-naa15.top")
top = Topology.from_path(top_path)
print(top.selected_molecule)
# %% [markdown]
# Please note, that selecting a molecule changes the name of its molecule type. For example: <br>
# 'Urea'  -->  'Urea(1x)'  (for one molecule as in 'urea.top') <br>
# 'Urea'  -->  'Urea(2x)'  (for two molecules as in 'urea-times-2.top') <br>
# 'NAA50' and 'NAA15'  -->  'NAA50(1x)-NAA15(1x)'  (for one molecule of each type as in 'naa50-naa15.top') <br>

# %% [markdown]
# It is also possible to manually select/deselect one or more molecule types by providing the selection function `is_selected_moleculetype_f`.
# In case of the 'naa50-naa15.top' example, we can select only the 'NAA15' molecule type by deselecting the 'NAA50' molecule type:
# %% [markdown]
# a) When creating the topology object directly
# %%
top_path = Path("naa50-naa15.top")
is_selected_moleculetype_f=get_is_selected_moleculetype_f(selected=[], deselected=["NAA50"])

top = Topology.from_path(top_path, 
is_selected_moleculetype_f=is_selected_moleculetype_f) 
print(top.selected_molecule)
# %% [markdown]
# b) When creating the topology object from a dictionary
# %%
top_dict = read_top(top_path) 
top = Topology(top_dict, is_selected_moleculetype_f=is_selected_moleculetype_f)
print(top.selected_molecule)
# %% [markdown]
# # 2. Writing a topology file
# Internally, the `Topology.to_path` method converts the `Topology` object to a dictionary using the `Topology.to_dict` method, and then writes the dictionary to a file using the `write_top` function. 
# Thus, the following two approaches are equivalent:
# %% [markdown]
# a) Write the topology to a file directly
# %%
out_path = Path("out.top")
top.to_path(out_path)
# %% [markdown]
# b) Write the topology to a file via explicit dictionary conversion first
# %%
out_path = Path("out.top")
out_dict = top.to_dict()
write_top(out_dict, out_path)
# %%