Simple Python module for parsing and visualising QTAIM critical points with py3dmol (best used in a jupyter notebook)

Extracts data from CPprop.txt multiwfn files and visualises points based on xyz structure.

Example usage
```python
from qtaim import QTAIM

# Initialize with a QTAIM output file
qtaim = QTAIM('CPprop.txt')

# Visualize using an XYZ file for atomic coordinates
qtaim.visualise(
    xyz_file='your_structure.xyz',
    show_lap=True,           # Show Laplacian at CPs
    show_bond_lengths=True,  # Display bond lengths at CPs
    print_parameters=True    # Print extracted parameter table
)
```
