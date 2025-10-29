# %%
from pathlib import Path

from gmx_top4py.topology.topology import Topology
# %% [markdown]
# # 1. Reading a topology file
# Let's read the GROMACS topology file 'urea.top' into a `Topology` object.
# %%
top_path = Path("urea.top") # path to the topology file
top = Topology.from_path(top_path) # read the topology file into a Topology object
print(top)

# %% [markdown]
# The first paragraph of the printout shows which molecule types and how many molecules of each type are present in the topology file. <br>

# <br>
# Topology with the following molecules: <br>
# Urea(1x): 1 <br>
# SOL: 1000 <br>
# <br>

# In case of the 'urea.top' example, there are one urea molecule and 1000 solvent molecules.
# The second paragraph of the printout provides an overview of the topology properties of the selected molecule type. <br>

# <br>
# With Moleculetype Urea(1x) with: <br>
# 8 atoms, <br>
# 7 bonds, <br>
# 0 angles, <br>
# 0 pairs, <br>
# 0 settles, <br>
# 0 exclusions, <br>
# 8 proper dihedrals <br>
# 3 improper dihedrals <br>
# 3 position restraints <br>
# 2 dihedral restraints <br>
# <br>

# Here the selected molecule type is 'Urea(1x)'. This molecule type has 8 atoms, 7 bonds, 8 proper dihedrals, 3 improper dihedrals, 3 position restraints, and 2 dihedral restraints.
# Last but not least, the third paragraph gives information about the force field that was used to parameterize the topology. <br>

# <br>
# ForceField parameters with <br>
# 65 atomtypes, <br>
# 96 bondtypes, <br>
# 230 angletypes, <br>
# 125 proper dihedraltypes <br>
# 49 improper dihedraltypes <br>
# 125 residuetypes <br>
# <br>

# In the 'urea.top' file, the amber99 force field is used.
# %% [markdown]
# ## 1.1. The Topology object
# An instance of the `Topology` class has the following main attributes:
# - moleculetypes: A dictionary containing the instances of the `MoleculeType` class for each molecule type present in the topology file.
# - ff: An instance of the `FF` class containing all parameters of the specified force field.
# %% [markdown]
# ### 1.1.1. The moleculetypes attribute
# The `moleculetypes` attribute is a dictionary where the keys are the names of the molecule types and the values are instances of the `MoleculeType` class.
# %% [markdown]
# In the case of the 'urea.top' example, there are two molecule types: 'SOL' and 'Urea(1x)':
# %%
print(top.moleculetypes)
# %% [markdown]
# For example, we can access the topology information of the 'Urea(1x)' molecule type as follows:
# %%
print(top.moleculetypes["Urea(1x)"])
# %% [markdown]
# An instance of the `MoleculeType` class has the following main attributes:
# - atoms: A dictionary containing for each atom the atom index as key and the corresponding instance of the `Atom` class as value.
# - bonds: A dictionary containing for each bond the tuple of atom indices as key and the corresponding instance of the `Bond` class as value.
# - angles: A dictionary containing for each angle the tuple of atom indices as key and the corresponding instance of the `Angle` class as value.
# - proper_dihedrals: A dictionary containing for each proper dihedral the tuple of atom indices as key and the corresponding instance of the `MultipleDihedrals` class as value.
# - improper_dihedrals: A dictionary containing for each improper dihedral the tuple of atom indices as key and the corresponding instance of the `MultipleDihedrals` class as value.
# - position_restraints: A dictionary containing for each position restraint the atom index as key and the corresponding instance of the `PositionRestraint` class as value.
# - dihedral_restraints: A dictionary containing for each dihedral restraint the tuple of atom indices as key and the corresponding instance of the `DihedralRestraint` class as value.
# %% [markdown]
# For example, we can access the atoms and bonds of the urea molecule as follows:
# %%
print(top.moleculetypes["Urea(1x)"].atoms) 
# %%
print(top.moleculetypes["Urea(1x)"].bonds)
# %% [markdown]
# ### 1.1.2. The ff attribute
# The ff attribute provides access to the force field parameters used in the topology.
# An instance of the `FF` class has the following main attributes:
# - atomtypes: A dictionary containing for each atom type the atom type name as key and the corresponding instance of the `AtomType` class as value.
# - bondtypes: A dictionary containing for each bond type the tuple of atom type names as key and the corresponding instance of the `BondType` class as value.
# - angletypes: A dictionary containing for each angle type the tuple of atom type names as key and the corresponding instance of the `AngleType` class as value
# - proper_dihedraltypes: A dictionary containing for each proper dihedral type the tuple of atom type names and the periodicity as key and the corresponding instance of the `DihedralType` class as value.
# - improper_dihedraltypes: A dictionary containing for each improper dihedral type the tuple of atom type names as key and the corresponding instance of the `DihedralType` class as value.
# - residuetypes: A dictionary containing for each residue type the residue name as key and the corresponding instance of the `ResidueType` class as value.
# %% [markdown]
# For example, we can access the atom types and bond types of the force field as follows
print(top.ff.atomtypes)
# %%
print(top.ff.bondtypes)
# %% [markdown]
# # 3. Writing a topology file
# The 'Topology.to_path' method can be used to write the topology to a file.
# %%
out_path = Path("urea_out.top")
top.to_path(out_path)
# %%
