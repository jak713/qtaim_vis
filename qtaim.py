import py3Dmol
import re
import math as m

class QTAIM:
    __version__ = "1.0.4"

    def __init__(self, file):
        self.file = file
        self.parameters = self.extract_qtaim(file)
	
    def extract_qtaim(self, file):
        parameters = {
            "CP_no": [],
            "type": [],
            "connected_atoms": [],
            "rho": [], 
            "laplacian": [],
            "energy_density": [],
            "potential_energy_density": [],
            "lagrangian_kinetic_energy": []
        }
        with open(file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if "CP" in line:
                parts = line.split()
                parameters["CP_no"].append(parts[2].replace(',', ''))
                parameters["type"].append(parts[4])
                parameters["connected_atoms"].append([])
            elif "Connected atoms" in line:
                connected_atoms = line.split()[2:]
                if parameters["connected_atoms"]:
                    parameters["connected_atoms"][-1] = connected_atoms
            elif "Density of all electrons" in line:
                parameters["rho"].append(float(line.split()[-1]))
            elif "Laplacian of electron density" in line:
                parameters["laplacian"].append(float(line.split()[-1]))
            elif "Energy density" in line:
                parameters["energy_density"].append(float(line.split()[-1]))
            elif "Potential energy density" in line:
                parameters["potential_energy_density"].append(float(line.split()[-1]))
            elif "Lagrangian kinetic energy" in line:
                parameters["lagrangian_kinetic_energy"].append(float(line.split()[-1]))

        for i, atoms in enumerate(parameters["connected_atoms"]):
            if len(atoms) >= 3:
                donor = atoms[0].replace('(', '').replace(')', '')
                acceptor = atoms[2].replace('(', '').replace(')', '')
                parameters["connected_atoms"][i] = donor + acceptor

        return parameters

    def is_same_atom_type(self, cp_no):
        #checks the two atom (ignoring numbering) types in connected_atoms for the given index (cp_no)
        #returns True if they are the same atom type, False otherwise
        atom_types = self.parameters["connected_atoms"][cp_no]
        if not atom_types:
            # print("No connected atoms found for CP", cp_no, "Type: " + self.parameters["type"][cp_no])
            return 
        elements = re.findall(r'[A-Za-z]+', atom_types)
        # print(elements)
        return elements[0] == elements[1] if len(elements) >= 2 else False
    
    def get_coordinates(self):
        """extract angstrom coordinates from self.file for each CP
        add x,y,z coordinates to the parameters dict"""
        parameters = self.parameters
        parameters["x_coord"] = []
        parameters["y_coord"] = []
        parameters["z_coord"] = []
        coordinates = []
        with open(self.file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if "Angstrom" in line:
                parts = line.split()
                if len(parts) >= 5:
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    coordinates.append((x, y, z))
        # add x_coord, y_coord, z_coord to parameters dict
        for i, cp_no in enumerate(parameters["CP_no"]):
            if i < len(coordinates):
                parameters["x_coord"].append(coordinates[i][0])
                parameters["y_coord"].append(coordinates[i][1])
                parameters["z_coord"].append(coordinates[i][2])
            else:
                parameters["x_coord"].append(None)
                parameters["y_coord"].append(None)
                parameters["z_coord"].append(None)
        return 


#################################################################
# --- Visualisation Code ---
#################################################################
    def visualise(self, xyz_file,view=None, display=True, show_cp = False, show_rho = False, show_lap = False, show_pos_lap = False, show_bond_lengths = False, aboveXangstrom = False, belowXangstrom =False, X = None, hide_ring_cage = False, show_only_same = False, show_only_different = False, show_atom_labels = False, connect_atoms_A_B = False, A=None, B=None, covalent = True, xyz_outline=False, print_parameters = False, legend = True, print_latex=False):
        """Visualise the bond critical points w/ Py3Dmol
            use xyz_file to get coordinates for connecting the CPs
            according to the connected_atoms list"""
        connection_indexes = []
        if print_latex:
            latex_headers_printed = False
        if print_parameters:
            headers_printed = False
            cp_of_interest_num = 0
        self.get_coordinates()
        if view is None:
            view = py3Dmol.view(width=1200, height=1000)


        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        coordinates = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                coordinates.append((x, y, z))
        
        # find indexes of connected atoms and use that to get the coordinates from the coordinates tuples list (offset index by 1 for python, connected_atoms is 1-indexed)
        for i, atoms in enumerate(self.parameters["connected_atoms"]):
            if connect_atoms_A_B and A is not None and B is not None:
                if A == B:
                    if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                        continue
                else:
                    if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                        continue  
            if show_only_different and self.is_same_atom_type(i):
                continue
            if show_only_same and not self.is_same_atom_type(i):
                continue
            if not atoms:
                continue
            atom_indexs = re.findall(r'\d+', atoms)
            if len(atom_indexs) >= 2:
                atom1 = int(atom_indexs[0]) - 1
                atom2 = int(atom_indexs[1]) - 1
                cp_no = int(self.parameters["CP_no"][self.parameters["connected_atoms"].index(atoms)]) - 1
                bcp_x = self.parameters["x_coord"][cp_no]
                bcp_y = self.parameters["y_coord"][cp_no]
                bcp_z = self.parameters["z_coord"][cp_no]
                coord1 = coordinates[atom1]
                coord2 = coordinates[atom2]
                dist1 = m.sqrt((coord1[0] - bcp_x)**2 + (coord1[1] - bcp_y)**2 + (coord1[2] - bcp_z)**2)
                dist2 = m.sqrt((coord2[0] - bcp_x)**2 + (coord2[1] - bcp_y)**2 + (coord2[2] - bcp_z)**2)
                length = dist1 + dist2
                if aboveXangstrom and X is not None:
                    if length < X:
                        continue
                if belowXangstrom and X is not None:
                    if length > X:
                        continue
                color = "green" if self.is_same_atom_type(cp_no) else "red"
                view.addLine({'start': {'x': coord1[0], 'y': coord1[1], 'z': coord1[2]},
                        'end':   {'x': bcp_x, 'y': bcp_y, 'z': bcp_z},
                        'color': color, 'radius': 5})
                view.addCylinder({'start': {'x': coord2[0], 'y': coord2[1], 'z': coord2[2]},
                        'end':   {'x': bcp_x, 'y': bcp_y, 'z': bcp_z},
                        'color': color, 'radius': 0.01})
                
                connection_indexes.append((atom1, atom2))

        # add bcps with different colors depending on the type of CP (blue for bond (3,-1), red for ring(3,+1), green for cage(3,+3), yellow for (3,-3))
        for i, cp_no in enumerate(self.parameters["CP_no"]):
            if self.parameters["x_coord"][i] is None:
                continue
            if connect_atoms_A_B and A is not None and B is not None:
                if A == B:
                    if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                        continue
                else:
                    if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                        continue  
            if hide_ring_cage and self.parameters["type"][i] in ["(3,+1)", "(3,+3)"]:
                continue
            if show_only_same and not self.is_same_atom_type(i):
                continue
            if show_only_different and self.is_same_atom_type(i):
                continue 
            x = self.parameters["x_coord"][i]
            y = self.parameters["y_coord"][i]
            z = self.parameters["z_coord"][i]
            connected = self.parameters["connected_atoms"][i]
            if not connected:
                continue
            atom_indexs = re.findall(r'\d+', connected)
            if len(atom_indexs) >= 2:
                atom1 = int(atom_indexs[0]) - 1
                atom2 = int(atom_indexs[1]) - 1
                cp_no = int(self.parameters["CP_no"][self.parameters["connected_atoms"].index(atoms)]) - 1
                # bcp_x = self.parameters["x_coord"][cp_no]
                # bcp_y = self.parameters["y_coord"][cp_no]
                # bcp_z = self.parameters["z_coord"][cp_no]
                coord1 = coordinates[atom1]
                coord2 = coordinates[atom2]
                dist1 = m.sqrt((coord1[0] - x)**2 + (coord1[1] - y)**2 + (coord1[2] - z)**2)
                dist2 = m.sqrt((coord2[0] - x)**2 + (coord2[1] - y)**2 + (coord2[2] - z)**2)
                length = dist1 + dist2
            if aboveXangstrom and X is not None:
                if length < X:
                    continue
            if belowXangstrom and X is not None:
                if length > X:
                    continue

            if not covalent:
                if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] < 0:
                    continue
                elif self.parameters["type"][i] == "(3,-3)":
                    continue
                color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0  else "red" if self.parameters["type"][i] == "(3,+1)" else "green" if  self.parameters["type"][i] == "(3,+3)" else "green" 
            else:    
                color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0 else "blue" if self.parameters["type"][i] == "(3,-1)" else "red" if self.parameters["type"][i] == "(3,+1)" else "green" if self.parameters["type"][i] == "(3,+3)" else "yellow"
            view.addSphere({'center': {'x': x, 'y': y, 'z': z}, 'radius': 0.1, 'color': color})
            
            if print_parameters:
                if not headers_printed:
                    print("{:<5} {:<10} {:<20} {:<10} {:<12} {:<15} {:<20} {:<20} {:<15}".format(
                        "CP", "Type", "ConnectedAtoms", "Rho", "Laplacian", "EnergyDensity", "PotentialEDensity", "LagrangianKcEnergy", "IntDist"
                    ))
                    headers_printed = True
                print("{:<5} {:<10} {:<20} {:<10.3f} {:<12.3f} {:<15.3f} {:<20.3f} {:<20.3f} {:<15}".format(
                    self.parameters["CP_no"][i],
                    self.parameters["type"][i],
                    self.parameters["connected_atoms"][i],
                    self.parameters["rho"][i],
                    self.parameters["laplacian"][i],
                    self.parameters["energy_density"][i],
                    self.parameters["potential_energy_density"][i],
                    self.parameters["lagrangian_kinetic_energy"][i],
                    f"{length:.2f}" if length else "N/A"
                ))
                cp_of_interest_num += 1

            if print_latex:
                # For LaTeX, print header if not printed already
                if not latex_headers_printed:
                    print("\\begin{table}[H]")
                    print("\\centering")
                    print("\\begin{tabular}{|c|c|c|c|c|c|}")
                    print("\\hline")
                    print("CP & Type & ConnectedAtoms & Rho & Laplacian & IntDist \\\\")
                    latex_headers_printed = True
                print("\\hline")
                print(f"{self.parameters['CP_no'][i]} & {self.parameters['type'][i]} & {self.parameters['connected_atoms'][i]} & {self.parameters['rho'][i]:.3f} & {self.parameters['laplacian'][i]:.3f} & {length:.2f} \\\\")
        
        if print_latex:
            print("\\hline")
            print("\\end{tabular}")
            print("\\caption{QTAIM parameters for CPs in the molecule.}")
            print("\\label{tab:qtaim}")
            print("\\end{table}")
        if print_parameters:
            print(f"Total number of interactions of interest: {cp_of_interest_num}")
        

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
            start_x = -5
            start_y = 5
            spacing = 1.0
            for i, item in enumerate(legend_items):
                y_pos = start_y - i * spacing
                if item["mode"] == "sphere":
                    sphere_pos = {'x': start_x, 'y': y_pos, 'z': 0}
                    label_pos = {'x': start_x + 0.5, 'y': y_pos, 'z': 0}
                    view.addSphere({'center': sphere_pos, 'radius': 0.1, 'color': item["color"]})
                    view.addLabel(item["label"], {
                        'position': label_pos,
                        'backgroundColor': 'white',
                        'backgroundOpacity': 0.5,
                        'fontSize': 10,
                        'fontColor': 'black',
                        'fontWeight': 'bold'
                })
                elif item["mode"] == "line":
                    line_start = {'x': start_x, 'y': y_pos, 'z': 0}
                    line_end = {'x': start_x + 0.5, 'y': y_pos, 'z': 0}
                    label_pos = {'x': start_x + 0.6, 'y': y_pos, 'z': 0}
                    view.addLine({'start': line_start, 'end': line_end, 'color': item["color"], 'radius': 5})
                    view.addLabel(item["label"], {
                        'position': label_pos,
                        'backgroundColor': 'white',
                        'backgroundOpacity': 0.5,
                        'fontSize': 10,
                        'fontColor': 'black',
                        'fontWeight': 'bold'
                    })
                
        if show_atom_labels:
            """show atom type labels for the atoms in the xyz file"""
            # extract atom types from the xyz file
            atom_types = []
            for line in lines[2:]:
                parts = line.split()
                if len(parts) >= 4:
                    atom_types.append(parts[0])
            # add labels for each atom in the xyz file
            for i, atom in enumerate(atom_types):
                if i < len(coordinates):
                    x = coordinates[i][0]
                    y = coordinates[i][1]
                    z = coordinates[i][2]
                    label = f"{atom}"
                    view.addLabel(label, {
                        'position': {'x': x, 'y': y, 'z': z},
                        'backgroundColor': 'white',
                        'backgroundOpacity': 0,
                        'fontSize': 10,
                        'fontColor': 'black',
                        'fontWeight': 'bold',
                        'attached': True  # Attach the label to the sphere
                    })

        if show_cp:
            for i, cp_no in enumerate(self.parameters["CP_no"]):
                if self.parameters["x_coord"][i] is None:
                    continue
                if hide_ring_cage and self.parameters["type"][i] in ["(3,+1)", "(3,+3)"]:
                    continue
                if show_only_same and not self.is_same_atom_type(i):
                    continue
                if show_only_different and self.is_same_atom_type(i):
                    continue
                if connect_atoms_A_B and A is not None and B is not None:
                    if A == B:
                        if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                            continue
                    else:
                        if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                            continue
                x = self.parameters["x_coord"][i]
                y = self.parameters["y_coord"][i]
                z = self.parameters["z_coord"][i]
                dist1 = m.sqrt((x - coordinates[0][0])**2 + (y - coordinates[0][1])**2 + (z - coordinates[0][2])**2)
                dist2 = m.sqrt((x - coordinates[1][0])**2 + (y - coordinates[1][1])**2 + (z - coordinates[1][2])**2)
                length = dist1 + dist2
                if aboveXangstrom and X is not None:
                    if length < X:
                        continue
                if belowXangstrom and X is not None:
                    if length > X:
                        continue
                label = f"{cp_no}"
                sphere_color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0 else \
                    "blue" if self.parameters["type"][i] == "(3,-1)" else \
                    "red" if self.parameters["type"][i] == "(3,+1)" else \
                    "green" if self.parameters["type"][i] == "(3,+3)" else "yellow"
                font_color = "black" if sphere_color == "yellow" or sphere_color == "lightblue"  else "white"
                view.addLabel(label, {
                    'position': {'x': x, 'y': y, 'z': z},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': 10,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })
                
        if show_rho:
            for i, cp_no in enumerate(self.parameters["CP_no"]):
                if self.parameters["x_coord"][i] is None:
                    continue
                if hide_ring_cage and self.parameters["type"][i] in ["(3,+1)", "(3,+3)"]:
                    continue
                if show_only_same and not self.is_same_atom_type(i):
                    continue
                if show_only_different and self.is_same_atom_type(i):
                    continue
                if connect_atoms_A_B and A is not None and B is not None:
                    if A == B:
                        if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                            continue
                    else:
                        if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                            continue  
                x = self.parameters["x_coord"][i]
                y = self.parameters["y_coord"][i]
                z = self.parameters["z_coord"][i]
                rho = self.parameters["rho"][i]
                label = f"ρ: {rho:.2f}"
                sphere_color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0 else \
                    "blue" if self.parameters["type"][i] == "(3,-1)" else \
                    "red" if self.parameters["type"][i] == "(3,+1)" else \
                    "green" if self.parameters["type"][i] == "(3,+3)" else "yellow"
                font_color = "black" if sphere_color == "yellow" or sphere_color == "lightblue" else "white"
                view.addLabel(label, {
                    'position': {'x': x, 'y': y, 'z': z},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': 10,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })
            
        if show_lap:
            for i, cp_no in enumerate(self.parameters["CP_no"]):
                if self.parameters["x_coord"][i] is None:
                    continue
                if hide_ring_cage and self.parameters["type"][i] in ["(3,+1)", "(3,+3)"]:
                    continue
                if show_only_same and not self.is_same_atom_type(i):
                    continue
                if show_only_different and self.is_same_atom_type(i):
                    continue
                if connect_atoms_A_B and A is not None and B is not None:
                    if A == B:
                        if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                            continue
                    else:
                        if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                            continue  
                x = self.parameters["x_coord"][i]
                y = self.parameters["y_coord"][i]
                z = self.parameters["z_coord"][i]
                laplacian = self.parameters["laplacian"][i]
                if laplacian < -1000:
                    continue
                label = f"∇²ρ: {laplacian:.2f}"
                sphere_color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0 else \
                    "blue" if self.parameters["type"][i] == "(3,-1)" else \
                    "red" if self.parameters["type"][i] == "(3,+1)" else \
                    "green" if self.parameters["type"][i] == "(3,+3)" else "yellow"
                font_color = "black" if sphere_color == "yellow" or sphere_color == "lightblue" else "white"
                view.addLabel(label, {
                    'position': {'x': x, 'y': y, 'z': z},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': 10,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })

        if show_pos_lap:
            for i, cp_no in enumerate(self.parameters["CP_no"]):
                if self.parameters["x_coord"][i] is None:
                    continue
                if hide_ring_cage and self.parameters["type"][i] in ["(3,+1)", "(3,+3)"]:
                    continue
                if show_only_same and not self.is_same_atom_type(i):
                    continue
                if show_only_different and self.is_same_atom_type(i):
                    continue
                if connect_atoms_A_B and A is not None and B is not None:
                    if A == B:
                        if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                            continue
                    elif A != B:
                        if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                            continue  
                x = self.parameters["x_coord"][i]
                y = self.parameters["y_coord"][i]
                z = self.parameters["z_coord"][i]
                laplacian = self.parameters["laplacian"][i]
                if laplacian < 0:
                    continue
                label = f"∇²ρ: {laplacian:.2f}"
                sphere_color = "lightblue" if self.parameters["type"][i] == "(3,-1)" and self.parameters["laplacian"][i] > 0 else \
                    "blue" if self.parameters["type"][i] == "(3,-1)" else \
                    "red" if self.parameters["type"][i] == "(3,+1)" else \
                    "green" if self.parameters["type"][i] == "(3,+3)" else "yellow"
                font_color = "black" if sphere_color == "yellow" or sphere_color == "lightblue" else "white"
                view.addLabel(label, {
                    'position': {'x': x, 'y': y, 'z': z},
                    'backgroundColor': sphere_color,
                    'backgroundOpacity': 0.3,
                    'fontSize': 10,
                    'fontColor': font_color,
                    'fontWeight': 'bold'
                })

        if show_bond_lengths:
            # calculate the sum of atom1-CP and CP-atom2 distances and print them at the CP
            for i, atoms in enumerate(self.parameters["connected_atoms"]):
                if show_only_different and self.is_same_atom_type(i):
                    continue
                if show_only_same and not self.is_same_atom_type(i):
                    continue
                if connect_atoms_A_B and A is not None and B is not None:
                    if A == B:
                        if set(re.findall(r'[A-Za-z]', str(self.parameters["connected_atoms"][i]))) != {A}:
                            continue
                    else:
                        if not (A in self.parameters["connected_atoms"][i] and B in self.parameters["connected_atoms"][i]):
                            continue  
                if not atoms:
                    continue
                atom_indexs = re.findall(r'\d+', atoms)
                if len(atom_indexs) >= 2:
                    atom1 = int(atom_indexs[0]) - 1
                    atom2 = int(atom_indexs[1]) - 1
                    cp_no = int(self.parameters["CP_no"][self.parameters["connected_atoms"].index(atoms)]) - 1
                    bcp_x = self.parameters["x_coord"][cp_no]
                    bcp_y = self.parameters["y_coord"][cp_no]
                    bcp_z = self.parameters["z_coord"][cp_no]
                    coord1 = coordinates[atom1]
                    coord2 = coordinates[atom2]
                    dist1 = m.sqrt((coord1[0] - bcp_x)**2 + (coord1[1] - bcp_y)**2 + (coord1[2] - bcp_z)**2)
                    dist2 = m.sqrt((coord2[0] - bcp_x)**2 + (coord2[1] - bcp_y)**2 + (coord2[2] - bcp_z)**2)
                    length = dist1 + dist2
                    if aboveXangstrom:
                        if length > X:
                            label = f"{length:.2f}"
                            view.addLabel(label, {
                                'position': {'x': bcp_x, 'y': bcp_y - 0.15, 'z': bcp_z},
                                'fontSize': 12,
                                'fontColor': 'black',
                                'fontWeight': 'bold',
                                'backgroundOpacity': 0.0
                            })
                    elif belowXangstrom:
                        if length < X:
                            label = f"{length:.2f}"
                            view.addLabel(label, {
                                'position': {'x': bcp_x, 'y': bcp_y - 0.15, 'z': bcp_z},
                                'fontSize': 12,
                                'fontColor': 'black',
                                'fontWeight': 'bold',
                                'backgroundOpacity': 0.0
                            })
                    else:
                        label = f"{length:.2f}"
                        view.addLabel(label, {
                            'position': {'x': bcp_x, 'y': bcp_y - 0.15, 'z': bcp_z},
                            'fontSize': 12,
                            'fontColor': 'black',
                            'fontWeight': 'bold',
                            'backgroundOpacity': 0.0
                        })
        if xyz_outline:
            # add xyz molecule structure into the view
            view.addModel(open(xyz_file, 'r').read(), 'xyz')
            view.setStyle({'stick': {'radius': 0.03}})
            

        view.setBackgroundColor('white')
        if display:
            view.zoomTo()
            view.show()
        connection_indexes = set(connection_indexes)
        return connection_indexes