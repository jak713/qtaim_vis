- `xyz_file: str`. This is required. Filename of the xyz file containing the structure. This has to be same structure as the one used for obtaining CPProp.txt otherwise the program may crash or the results will be nonsensical.

- `view=None`. This can be ignored. It becomes important when combining nbo and qtaim visualisation together.

- `display=True`. This can also be ignored. It becomes important when combining nbo and qtaim visualisation together.

- `show_cp=False`. Show critical point numbers (from CPprop).

- `show_rho=False`. Show electron density values.

- `show_lap=False`. Show Laplacian values.

- `show_pos_lap=False`. Show only BCPs with a positive Laplacian value (indicative of a non-covalent interaction).

- `hide_bcp=False`. Hide bond critical points.

- `show_bond_lengths=False`. Show bond lengths in Angstrom between identified bond paths.

- `aboveXangstrom=False`. Show paths with bond length above X (see below).

- `belowXangstrom=False`. Show paths with bond length below X (see below).

- `X=None`. Set value of X angstrom (see above).

- `hide_ring_cage=False`. Hide BCPs of type (3,+3) and (3,+1)

- `show_only_same=False`. Show only BCP identified between the same atoms.

- `show_only_different=False`. Show only BCP identified between differing atoms.

- `show_atom_labels=False`. Show element symbol at atom location.

- `connect_atoms_A_B=False`. Show only paths between specfic atoms A and B (see below)

- `A=None`. Set atom A (see above)

- `B=None`.Set atom B (see above)

- `covalent=True`. Show (3,-3)-type CPs (nuclear) and non-positive BCPs.

- `xyz_outline=False`. Show the structural outline of the molecule underneath the CPs.

- `print_parameters=False`. Print the table of identified connections inclusive of all available data.

- `legend=True`. Print a legend for the different components of the visualisation.

- `print_latex=False`. Print the table of identified connections with LaTeX formatting, ready to copy-paste.

- `make_proportional=False`. Set the sphere radius of BCPs proportional to its relative rho value.