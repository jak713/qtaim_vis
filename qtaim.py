import py3Dmol
import re
import warnings
import numpy as np
from ase.io import read
from ase import Atoms


class CP():
    """Class to represent a Critical Point (CP) in QTAIM analysis."""
    def __init__(self, cp_block):
        self.cp_block = cp_block
        self.properties = self.extract_properties()

    def __repr__(self):
        """String representation of the CP object."""
        return f"CP(CP_no={self.cp_no}, type={self.cp_type}, connected_atoms={self.connected_atoms}, rho={self.rho}, laplacian={self.laplacian}, energy_density={self.energy_density}, potential_energy_density={self.potential_energy_density}, lagrangian_kinetic_energy={self.lagrangian_kinetic_energy}, coordinates={self.coordinates})"

    def extract_properties(self) -> dict:
        """Extracts properties from the CP block.

        Returns:
            dict: A dictionary containing the properties of the CP.
        """

        properties = {

            "CP_no": None,
            "type": None,
            "connected_atoms": None,
            "rho": None,
            "laplacian": None,
            "energy_density": None,
            "potential_energy_density": None,
            "lagrangian_kinetic_energy": None,
            "coordinates": None
        }
        # Define regex patterns for extracting properties
        header_pattern = re.compile(r'\-*\s*CP\s*(\d+),\s*Type\s*\(3,\s*([\+\-]\d+)\)\s*\-*\n')
        connected_atoms_pattern = re.compile(r'Connected atoms\:\s*(\d+)\((\w{1,2})\s{0,1}\)\s*\-{2}\s*(\d+)\((\w{1,2})\s{0,1}\)')
        position_pattern = re.compile(r'Position\s*\(Angstrom\)\:\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)')
        density_pattern = re.compile(r'Density of all electrons\:\s*(-?\d+\.\d+)')
        laplacian_pattern = re.compile(r'Laplacian of electron density\:\s*(-?\d+\.\d+)')
        energy_density_pattern = re.compile(r'Energy density E\(r\) or H\(r\)\:\s*(-?\d+\.\d+)')
        potential_energy_density_pattern = re.compile(r'Potential energy density V\(r\)\:\s*(-?\d+\.\d+)')
        lagrangian_kinetic_energy_pattern = re.compile(r'Lagrangian kinetic energy G\(r\)\:\s*(-?\d+\.\d+)')

        # Extract CP number and type

        header_match = re.search(header_pattern, self.cp_block)
        cp_num = int(header_match.group(1))
        cp_type = int(header_match.group(2))

        if cp_type == -1:
            connected_atoms = re.findall(connected_atoms_pattern, self.cp_block)
            connected_atom_id1 = int(connected_atoms[0][0])
            connected_atom_type1 = connected_atoms[0][1]
            connected_atom_id2 = int(connected_atoms[0][2])
            connected_atom_type2 = connected_atoms[0][3]
            connected_atoms = (connected_atom_id1, connected_atom_type1, connected_atom_id2, connected_atom_type2)
        else:
            connected_atoms = None
        # Extract coordinates
        position_match = re.search(position_pattern, self.cp_block)
        coordinates = np.array([float(position_match.group(1)), float(position_match.group(2)), float(position_match.group(3))])
        # Extract density
        density_match = re.search(density_pattern, self.cp_block)
        try:
            density = float(density_match.group(1))
            properties["rho"] = density
        except AttributeError:
            warnings.warn("Density not found in CP block.", UserWarning)
            density = None
        try:
            # Extract laplacian
            laplacian_match = re.search(laplacian_pattern, self.cp_block)
            laplacian = float(laplacian_match.group(1))
        except AttributeError:
            warnings.warn("Laplacian not found in CP block.", UserWarning)
            laplacian = None
        # Extract energy density
        energy_density_match = re.search(energy_density_pattern, self.cp_block)
        try:
            energy_density = float(energy_density_match.group(1))
        except AttributeError:
            warnings.warn("Energy density not found in CP block.", UserWarning)
            energy_density = None
        # Extract potential energy density
        potential_energy_density_match = re.search(potential_energy_density_pattern, self.cp_block)
        try:
            potential_energy_density = float(potential_energy_density_match.group(1))
        except AttributeError:
            warnings.warn("Potential energy density not found in CP block.", UserWarning)
            potential_energy_density = None
        # Extract lagrangian kinetic energy
        lagrangian_kinetic_energy_match = re.search(lagrangian_kinetic_energy_pattern, self.cp_block)
        try:
            lagrangian_kinetic_energy = float(lagrangian_kinetic_energy_match.group(1))
        except AttributeError:
            warnings.warn("Lagrangian kinetic energy not found in CP block.", UserWarning)
            lagrangian_kinetic_energy = None
        # Store properties in the dictionary
        properties["CP_no"] = cp_num
        properties["type"] = cp_type
        properties["connected_atoms"] = connected_atoms
        properties["rho"] = density
        properties["laplacian"] = laplacian
        properties["energy_density"] = energy_density
        properties["potential_energy_density"] = potential_energy_density
        properties["lagrangian_kinetic_energy"] = lagrangian_kinetic_energy
        properties["coordinates"] = coordinates
        return properties

    @property
    def coordinates(self):
        """Returns the coordinates of the CP."""
        return self.properties["coordinates"]

    @property
    def cp_no(self):
        """Returns the CP number."""
        return self.properties["CP_no"]

    @property
    def cp_type(self):
        """Returns the CP type."""
        return self.properties['type']

    @property
    def connected_atoms(self):
        """Returns the connected atom."""
        return self.properties["connected_atoms"]

    @property
    def rho(self):
        """Returns the density."""
        return self.properties["rho"]

    @property
    def laplacian(self):
        """Returns the laplacian."""
        return self.properties["laplacian"]

    @property
    def energy_density(self):
        """Returns the energy density."""
        return self.properties["energy_density"]

    @property
    def potential_energy_density(self):
        """Returns the potential energy density."""
        return self.properties["potential_energy_density"]

    @property
    def lagrangian_kinetic_energy(self):
        """Returns the lagrangian kinetic energy."""
        return self.properties["lagrangian_kinetic_energy"]


class QTAIM:
    __version__ = "1.0.4"

    def __init__(self, file):
        self.file = file
        self.cps = self.get_cps()

    def __repr__(self):
        """String representation of the QTAIM object."""
        return f"QTAIM(file={self.file}, number of cps={len(self.cps)},"


    def get_cps(self) -> list[CP]:
        """Get the CPs from the file.

        Returns:
            list[CP]: A list of CP objects.
        """
        with open(self.file, 'r') as f:
            text = f.read()
        # split text into cp blocks
        cp_block_pattern = re.compile(r'(\-*\s*CP\s*\d+,\s*Type\s*\(3,\s*[\+\-]\d+\)\s*\-*\n)')
        # split the text block
        cp_blocks = re.split(cp_block_pattern, text)[1:]
        cp_blocks = [cp_blocks[i] + cp_blocks[i + 1] for i in range(0, len(cp_blocks), 2)]
        # create a list of CP objects
        cp_objects = []
        for cp_block in cp_blocks:
            cp = CP(cp_block)
            cp_objects.append(cp)
        return cp_objects

    def is_same_atom_type(self, cp_no) -> bool:
        """Check if the two connected atoms are of the same type.
        Args:
            cp_no (int): The index of the CP to check.
        Returns:
            bool: True if the two connected atoms are of the same type, False otherwise.
        """
        # checks the two atom (ignoring numbering) types in connected_atoms for the given index (cp_no)
        # returns True if they are the same atom type, False otherwise
        atom_types = self.cps[cp_no - 1].connected_atoms
        if atom_types is None:
            return False
        if atom_types[1] == atom_types[3]:
            return True
        else:
            return False

    def get_connecting_cps_indices(self) -> list:
        """Get CPs connecting two atoms."""
        # get the connected atoms for each CP
        # returns a list of tuples with the atom types and their indexes
        connection_cps = []
        for cp in self.cps:
            if cp.connected_atoms is not None:
                connection_cps.append(cp.cp_no)

        return connection_cps

    @property
    def connecting_cps(self):
        return self.get_connecting_cps_indices()

#################################################################
# --- Visualisation Code ---
#################################################################


    def get_distance(self, bcp: CP, atoms: Atoms) -> float:
       """Calculate the distance between a BCP and the atoms in the molecule.
       Args:
           bcp (CP): The BCP object.
           atoms (Atoms): The Atoms object representing the molecule.
       Returns:
           float: The distance between the BCP and the atoms in the molecule.
       """
       # Get the coordinates of the BCP
       bcp_coords = bcp.coordinates
       # get indices of the connected atoms
       index_1 = int(bcp.connected_atoms[0]) - 1
       index_2 = int(bcp.connected_atoms[2]) - 1
       atom_type_1 = bcp.connected_atoms[1]
       atom_type_2 = bcp.connected_atoms[3]
       assert atom_type_1 == atoms[index_1].symbol
       assert atom_type_2 == atoms[index_2].symbol
       # Get the coordinates of the atoms in the molecule
       coords_1 = atoms[index_1].position
       coords_2 = atoms[index_2].position
       # Calculate the distance between the BCP and the atoms
       bcp_coord1_dist = np.linalg.norm(bcp_coords - coords_1)
       bcp_coord2_dist = np.linalg.norm(bcp_coords - coords_2)
       # Calculate the sum of the distances
       sum_dist = bcp_coord1_dist + bcp_coord2_dist
       return bcp_coord1_dist, bcp_coord2_dist, sum_dist

    def visualise(self, xyz_file, show_cp=False, show_rho=False, show_lap=False,
                show_pos_lap=False, show_bond_lengths=False, aboveXangstrom: float = None,
                belowXangstrom: float = None, hide_ring_cage=False, show_only_same=False,
                show_only_different=False, show_atom_labels=False, connect_atoms_A_B=False,
                A: str = None, B: str = None, covalent=True, xyz_outline=False, print_parameters=False,
                legend=True, print_latex=False, fontsize: int=12, connection_radius: float=0.01,
                sphere_radius: float=0.1, show_nucleus_cp: bool=False, 
                ):
        """Visualise the bond critical points w/ Py3Dmol
            use xyz_file to get coordinates for connecting the CPs
            according to the connected_atoms list
        Args:
            xyz_file (str): The path to the xyz file.
            show_cp (bool): Whether to show the CPs or not.
            show_rho (bool): Whether to show the density or not.
            show_lap (bool): Whether to show the laplacian or not.
            show_pos_lap (bool): Whether to show the positive laplacian or not.
            show_bond_lengths (bool): Whether to show the bond lengths or not.
            aboveXangstrom (float): The minimum distance between BCPs and atoms.
            belowXangstrom (float): The maximum distance between BCPs and atoms.
            hide_ring_cage (bool): Whether to hide the ring and cage CPs or not.
            show_only_same (bool): Whether to show BCPs with the same atom type or not.
            show_only_different (bool): Whether to show only the BCPs with different atom types or not.
            show_atom_labels (bool): Whether to show the atom labels or not.
            connect_atoms_A_B (bool): Whether show BCPs between elements A and B or not.
            A (str): The first element to connect.
            B (str): The second element to connect.
            covalent (bool): Whether to show covalent BCPs or not.
            xyz_outline (bool): Whether to show the xyz outline or not.
            print_parameters (bool): Whether to print the parameters or not.
            legend (bool): Whether to show the legend or not.
            print_latex (bool): Whether to print the parameters in LaTeX format or not.
            fontsize (int): The font size of the labels.
            connection_radius (float): The radius of the connection lines.
        Returns:
            None
            """
        label_shift = fontsize * 0.025
        connection_radius = 0.01
        if print_latex:
            latex_headers_printed = False
        if print_parameters:
            headers_printed = False
            cp_of_interest_num = 0

        view = py3Dmol.view(width=1200, height=1000)
        visualised_cps = []
        atoms = read(xyz_file)
        legend_items_show = set()

        for cp in self.cps:
            label_off_set = 0.02
            # skip CPs that are not of interest
            if (cp.cp_type == 1 or cp.cp_type == 3) and hide_ring_cage:
                continue
            if show_only_same and not self.is_same_atom_type(cp.cp_no):
                continue
            if show_only_different and self.is_same_atom_type(cp.cp_no):
                continue
            if not show_nucleus_cp and cp.cp_type == -3:
                continue
            sphere_color = "black"
            font_color = "white"
            if not covalent and cp.cp_type == -1 and cp.laplacian < 0:
                continue
            elif covalent and cp.cp_type == -1 and cp.laplacian < 0:
                sphere_color = "blue"
                font_color = "white"
            elif cp.cp_type == -1 and cp.laplacian < 0:
                continue
            elif cp.cp_type == -3:
                sphere_color = "yellow"
                font_color = "black"
            elif cp.cp_type == -1 and cp.laplacian > 0:
                sphere_color = "lightblue"
                font_color = "black"
            elif cp.cp_type == +1:
                sphere_color = "red"
                font_color = "white"
            elif cp.cp_type == +3:
                sphere_color = "green"
                font_color = "white"

            # add sphere for the CP
            if cp.cp_type == -1:
                _,_, length = self.get_distance(cp, atoms)
                if aboveXangstrom and length > aboveXangstrom:
                    continue
                if belowXangstrom and length < belowXangstrom:
                    continue
                if connect_atoms_A_B and A and B:
                    assert A in atoms.get_chemical_symbols()
                    assert B in atoms.get_chemical_symbols()
                    atom_1 = atoms[int(cp.connected_atoms[0]) - 1].symbol
                    atom_2 = atoms[int(cp.connected_atoms[2]) - 1].symbol
                    if (atom_1 == A and atom_2 == B) or (atom_1 == B and atom_2 == A):
                        view.addCylinder({
                            'start': {'x': cp.coordinates[0], 'y': cp.coordinates[1], 'z': cp.coordinates[2]},
                            'end': {'x': atoms[int(cp.connected_atoms[0]) - 1].position[0],
                                    'y': atoms[int(cp.connected_atoms[0]) - 1].position[1],
                                    'z': atoms[int(cp.connected_atoms[0]) - 1].position[2]},
                            'color': sphere_color,
                            'radius': connection_radius
                        })
                        view.addCylinder({
                            'start': {'x': cp.coordinates[0], 'y': cp.coordinates[1], 'z': cp.coordinates[2]},
                            'end': {'x': atoms[int(cp.connected_atoms[2]) - 1].position[0],
                                    'y': atoms[int(cp.connected_atoms[2]) - 1].position[1],
                                    'z': atoms[int(cp.connected_atoms[2]) - 1].position[2]},
                            'color': sphere_color,
                            'radius': connection_radius
                        })
                    else:
                        continue

                if show_bond_lengths:
                    label = f"Length: {length:.2f} Å"
                    label_off_set += label_shift
                    view.addLabel(label, {
                        'position': {'x': cp.coordinates[0], 'y': cp.coordinates[1] + label_off_set, 'z': cp.coordinates[2]},
                        'backgroundColor': "black",
                        'backgroundOpacity': 0.3,
                        'fontSize': fontsize,
                        'fontColor': "white",
                        'fontWeight': 'bold'
                    })
                view.addCylinder({
                    'start': {'x': cp.coordinates[0], 'y': cp.coordinates[1], 'z': cp.coordinates[2]},
                    'end': {'x': atoms[int(cp.connected_atoms[0]) - 1].position[0],
                            'y': atoms[int(cp.connected_atoms[0]) - 1].position[1],
                            'z': atoms[int(cp.connected_atoms[0]) - 1].position[2]},
                    'color': sphere_color,
                    'radius': connection_radius
                })
                view.addCylinder({
                    'start': {'x': cp.coordinates[0], 'y': cp.coordinates[1], 'z': cp.coordinates[2]},
                    'end': {'x': atoms[int(cp.connected_atoms[2]) - 1].position[0],
                            'y': atoms[int(cp.connected_atoms[2]) - 1].position[1],
                            'z': atoms[int(cp.connected_atoms[2]) - 1].position[2]},
                    'color': sphere_color,
                    'radius': connection_radius
                })
                
            view.addSphere({'center': {'x': cp.coordinates[0], 'y': cp.coordinates[1], 'z': cp.coordinates[2]},
                            'radius': sphere_radius, 'color': sphere_color})
            visualised_cps.append(cp.cp_no)
            if cp.cp_type == -1 and cp.laplacian < 0:
                legend_items_show.add(0)
            elif cp.cp_type == -1 and cp.laplacian > 0:
                legend_items_show.add(1)
            elif cp.cp_type == +1:
                legend_items_show.add(2)
            elif cp.cp_type == +3:
                legend_items_show.add(3)
            elif cp.cp_type == -3:
                legend_items_show.add(4)
            if self.is_same_atom_type(cp.cp_no):
                legend_items_show.add(5)
            else:
                legend_items_show.add(6)

            if show_rho:
                label = f"ρ: {cp.rho:.2f}"
                label_off_set += label_shift
                view.addLabel(label, {
                    'position': {'x': cp.coordinates[0], 'y': cp.coordinates[1] + label_off_set, 'z': cp.coordinates[2]},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': fontsize,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })
            if show_lap:
                label = f"∇²ρ: {cp.laplacian:.2f}"
                label_off_set += label_shift
                view.addLabel(label, {
                    'position': {'x': cp.coordinates[0], 'y': cp.coordinates[1] + label_off_set, 'z': cp.coordinates[2]},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': fontsize,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })
            if show_pos_lap and cp.laplacian > 0:
                label = f"∇²ρ: {cp.laplacian:.2f}"
                label_off_set += label_shift
                view.addLabel(label, {
                    'position': {'x': cp.coordinates[0], 'y': cp.coordinates[1] + label_off_set, 'z': cp.coordinates[2]},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': fontsize,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })
            if show_pos_lap and cp.laplacian < 0:
                continue

        # legend for the different colors, positioned to avoid overlap
        if legend:
            legend_items = [
                {"label": "Bond Critical Point", "color": "blue", "mode": "sphere"},
                {"label": "Bond Critical Point (positive Laplacian)", "color": "lightblue", "mode": "sphere"},
                {"label": "Ring Critical Point", "color": "red", "mode": "sphere"},
                {"label": "Cage Critical Point", "color": "green", "mode": "sphere"},
                {"label": "Nuclear Critical Point", "color": "yellow", "mode": "sphere"},
                {"label": "Same Connecting Atoms", "color": "green", "mode": "line"},
                {"label": "Different Connecting Atoms", "color": "red", "mode": "line"}
            ]

            # set starting position and spacing offsets for the legend items
            start_x = min(atoms.get_positions()[:, 0]) - 1.5
            start_y = max(atoms.get_positions()[:, 1]) + 1.5
            spacing = 1.0
            for counter, i in enumerate(legend_items_show):
                y_pos = start_y - counter * spacing
                if legend_items[i]["mode"] == "sphere":
                    sphere_pos = {'x': start_x, 'y': y_pos, 'z': 0}
                    label_pos = {'x': start_x + 0.5, 'y': y_pos, 'z': 0}
                    view.addSphere({'center': sphere_pos, 'radius': 0.1, 'color': legend_items[i]["color"]})
                    view.addLabel(legend_items[i]["label"], {
                        'position': label_pos,
                        'backgroundColor': 'white',
                        'backgroundOpacity': 0.5,
                        'fontSize': fontsize,
                        'fontColor': 'black',
                        'fontWeight': 'bold'
                })
        if show_atom_labels:
            """show atom type labels for the atoms in the xyz file"""
            # extract atom types from the xyz file
            for atom in atoms:
                label = atom.symbol
                view.addLabel(label, {
                    'position': {'x': atom.position[0], 'y': atom.position[1], 'z': atom.position[2]},
                    'backgroundColor': 'white',
                    'backgroundOpacity': 0,
                    'fontSize': fontsize,
                    'fontColor': 'black',
                    'fontWeight': 'bold',
                    'attached': True  # Attach the label to the sphere
                })

        if xyz_outline:
            # add xyz molecule structure into the view
            view.addModel(open(xyz_file, 'r').read(), 'xyz')
            view.setStyle({'stick': {'radius': 0.03}})

        view.setBackgroundColor('white')
        print("QTAIM analysis complete.")
        print(f"Number of total CPs: {len(self.cps)}")
        print(f"Number of visualised CPs: {len(visualised_cps)}")
        view.show()
