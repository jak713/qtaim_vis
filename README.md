This requires a CPprop.txt file from multiwfn. To obtain this starting with Orca, you can follow these steps:
1. Single point calculation with desired level of theory
2. From calculation.gbw file, use the command below to get wavefunction input for multiwfn:
```bash
orca_2mkl calculation -molden
```
3. Run multiwfn with the new calculation.molden.input file and select these options:
```bash
multiwfn calculation.molden.input
2 Topology analysis
2 Search CPs from nuclear positions
3 Search CPs from midpoint of atomic pairs
8 Generating the paths connecting (3,-3) and (3,-1) CPs
9 Generating the paths connecting (3,+1) and (3,+3) CPs
7 Show real space function values at specific CP or all CPs
0
```
4. You have a CPprop.txt file ready after some time.


Example usage
```python
from qtaim import QTAIM

qtaim = QTAIM('CPprop.txt')

# Visualize using an XYZ file (same as the one used for wavefunction for multiwfn) for atomic coordinates
qtaim.visualise("your_structure.xyz",
                    xyz_outline=True, # shows xyz file as stick model (useful)
                    covalent=False, # hides (3,-3)-type CPs (nuclear)and non-positive BCPs
                    show_rho=False, 
                    show_pos_lap=True,  # generally what one might want for non-covalent interactions
                    hide_ring_cage=True, # hides BCPs of type (3,+3) and (3,+1)
                    show_atom_labels = False, # atom type andn number
                    show_only_same=False, # same connecting atoms like 1C2C C -- C
                    show_only_different=False, # different connecting atoms like 201N59H for N -- H
                    show_bond_lengths=False, 
                    print_parameters=True, # print a table with all you want to know for only the displayed (as per these options) BCPs
                    legend=True, # might not need this
                    connect_atoms_A_B=False, # when wanting to study particular interactions e.g. O and H for hydrogen bonding
                    A="C",
                    B="H")
```
