# %%
from pathlib import Path
from copy import deepcopy

from gmx_top4py.topology.topology import Topology
from gmx_top4py.parsing import read_top, write_top
from gmx_top4py.parameterizing import Parameterizer
from gmx_top4py.topology.atomic import Dihedral

# Read topology file and instantiate Topology object
top_dict = read_top(Path("hexala.top"), ffdir=Path("amber99sb-star-ildnp.ff"))
original_top = Topology(top_dict)

# Write topology to file
original_top.to_path("hexala_out.top")
# %%
# Manipulate partial charges of individual atoms
# Example: Change the charge of the CB atom of the first Ala residue and its attached H atoms
top = deepcopy(original_top)
print("Atomic charges before modification:")
for atom_id in ["11", "12", "13", "14"]:
    print(top.atoms[atom_id])

top.atoms["11"].charge = -0.3  # CB atom of the first ala residue
top.atoms["12"].charge = 0.1  # HB1 atom of the first ala residue
top.atoms["13"].charge = 0.1  # HB2 atom
top.atoms["14"].charge = 0.1  # HB3 atom

print("\nAtomic charges after modification:")
for atom_id in ["11", "12", "13", "14"]:
    print(top.atoms[atom_id])
# %%
# Simillarly the parameters can be modified of a single bond, angle, dihedral, etc.
# Example: Change the bond length and force constant of the bond between the CB atom and the CA atom of the first Ala residue
top = deepcopy(original_top)
print("Bond information before modification:")
print(
    top.bonds[("9", "11")]
)  # Bond between the CA atom (= atom 9) and the CB atom (= atom 9) of the first Ala residue

print(
    "\nAs the equilibrium bond length (c0) and the force constant (c1) are None, the bonded parameters for this CA-CB bond are taken from the bond type in the force field:"
)
print(top.ff.bondtypes[("CA", "CB")])

top.bonds[("9", "11")].c0 = 0.15  # New bond length in nm
top.bonds[("9", "11")].c1 = 500000.0  # New force constant in kJ mol-1 nm-2

print("\nBond information after modification:")
print(top.bonds[("9", "11")])

# %%
# We can also modify the parameters of all bonds of a certain type
# by overriding the parameters of the bond type in the force field
# Example: Change the bond length and force constant of all CA-CB bonds
top = deepcopy(original_top)
print("CA-CB bond parameters from the force field:")
print(top.ff.bondtypes[("CA", "CB")])

top.ff.bondtypes[("CA", "CB")].c0 = 0.15  # New bond length in nm
top.ff.bondtypes[("CA", "CB")].c1 = 500000.0  # New force constant in kJ mol-1 nm-2

print("\nCA-CB bond parameters after modification:")
print(top.ff.bondtypes[("CA", "CB")])


# %%
# For more complex modifications of the topology, one can implement a custom Parameterizer class
# Example: Change the N-CA-CB-HBX proper dihedral parameters only for Ala residues
class MyParameterizer(Parameterizer):
    def parameterize_topology(
        self, current_topology: Topology, focus_nrs: set[str] | None = None
    ) -> Topology:
        # iterate over all proper dihedrals
        for proper in current_topology.proper_dihedrals.values():
            # Only modify proper dihedrals in Ala residues
            if current_topology.atoms[proper.ai].residue == "ALA":
                # Get the atom types of the four atoms in the proper dihedral
                atom_names = (
                    current_topology.atoms[proper.ai].atom,
                    current_topology.atoms[proper.aj].atom,
                    current_topology.atoms[proper.ak].atom,
                    current_topology.atoms[proper.al].atom,
                )
                if atom_names in [
                    ("N", "CA", "CB", "HB1"),
                    ("N", "CA", "CB", "HB2"),
                    ("N", "CA", "CB", "HB3"),
                ]:
                    # Change the parameters of the N-CA-CB-HBX proper dihedral
                    proper.dihedrals[""] = Dihedral(
                        proper.ai,
                        proper.aj,
                        proper.ak,
                        proper.al,
                        funct="9",
                        c0="180.0",
                        c1="2.0",
                    )
        return current_topology


top = deepcopy(original_top)
# Assign the new parameterizer to the topology
top.parametrizer = MyParameterizer()
# Please note that the parameters are only updated when needs_parameterization is set to True
top.needs_parameterization = True
# Now we can update the parameters
# a) Directly in the topology instance by calling update_parameters() explicitly
top.update_parameters()
# b) when writing the topology to dict
top_dict = top.to_dict()
# c) or when writing the topology to a file
top.to_path("hexala_modified_proper.top")
# When executing to_dict() or to_path(), the parameters are automatically updated by calling internally update_parameters() if needs_parameterization is True
