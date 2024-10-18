#!/usr/bin/env python3
"""This module is used to extract neighbors of an atom up to the 4th neighbor.

Made by C. Cattin, 2024.
"""

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from sulley.symmetry import create_graph, split_disconnected_graphs
import networkx as nx

def load_molecule(smiles: str):
    """Load a molecule from a SMILES string.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol

def load_molecule_from_sdf(filename: str):
    """Load a molecule from a file.

    Parameters
    ----------
    filename : str
        Name of the .sdf file.

    Returns
    -------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    """
    mol = Chem.SDMolSupplier(filename, removeHs=False)[0]
    return mol

def parse_tinker_xyz(filename: str):
    """Parse a Tinker .xyz file.

    Parameters
    ----------
    filename : str
        Name of the .xyz file.

    Returns
    -------
    atom : list
        List of atoms.
    bonds : list
        List of bonds.

    Raises
    ------
    ValueError
        If the number of atoms in the .xyz file does not match the number of atoms in the header.
    """
    atom = []
    bonds = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                num_atoms = int(line.split()[0])
                continue
            if line == '\n':
                continue
            line = line.strip().split()
            atom_idx = int(line[0]) - 1
            atom_type = line[1]
            x, y, z = map(float, line[2:5])
            neighbors = [int(neigh) - 1 for neigh in line[6:]]
            atom.append((atom_type, (x, y, z), atom_idx))

            for neigh in neighbors:
                # Avoid duplicates
                if neigh > atom_idx:
                    bonds.append((atom_idx, neigh))
    
    if len(atom) != num_atoms:
        raise ValueError(f"Number of atoms in the .xyz file ({len(atom)}) does not match the number of atoms in the header ({num_atoms})")

    return atom, bonds

def load_molecule_from_tinker_xyz(filename: str):
    """Load a molecule from a Tinker .xyz file.

    Parameters
    ----------
    filename : str
        Name of the .xyz file.

    Returns
    -------
    mol : rdkit.Chem.Mol or list
        RDKit molecule object or list of RDKit molecule objects.
    """

    atom_type_to_atomic_num = {
        'H': 1, 'He': 2,
        'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
        'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
        'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
        'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
        'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
        'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
        'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
        'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
        'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
        'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
        'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
        'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
        'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
    }

    atom, bonds = parse_tinker_xyz(filename)
    mol = Chem.RWMol()

    # Add atoms to the molecule
    xyz_to_rdkit_idx = {}
    for atom_type, coords, xyz_idx in atom:
        atomic_num = atom_type_to_atomic_num[atom_type]
        idx = mol.AddAtom(Chem.Atom(atomic_num))
        xyz_to_rdkit_idx[xyz_idx] = idx

    # Add bonds to the molecule
    for idx1, idx2 in bonds:
        idx1 = xyz_to_rdkit_idx[idx1]
        idx2 = xyz_to_rdkit_idx[idx2]
        mol.AddBond(idx1, idx2, Chem.rdchem.BondType.SINGLE)

    # Set 3D coordinates
    conf = Chem.Conformer(len(atom))
    for i, (_, (x, y, z), index) in enumerate(atom):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.UpdatePropertyCache(strict=False)

    
    # If multiple molecules are present return list of molecules
    if not is_single_molecule(mol):
        graph = create_graph(mol)
        sub_graph = split_disconnected_graphs(graph)
        mols = []
        last_node_idx = 0
        for graph in sub_graph:
            mols.append(
                load_molecule_from_graph(graph, shift=last_node_idx)
            )
            last_node_idx += len(graph.nodes)

        return mols
        
    return mol.GetMol()

def is_single_molecule(mol):
    """Check if a molecule is a single molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.

    Returns
    -------
    bool
        True if the molecule is a single molecule.
    """


    graph = create_graph(mol)
    return nx.is_connected(graph)


def load_molecule_from_graph(graph, shift = 0):
    """Load a molecule from a graph.

    Parameters
    ----------
    graph : nx.Graph
        NetworkX graph object.

    Returns
    -------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    """
    mol = Chem.RWMol()
    for node in graph.nodes(data=True):
        atomic_num = node[1]['atomic_num']
        mol.AddAtom(Chem.Atom(atomic_num))
    for edge in graph.edges():
        mol.AddBond(
            edge[0] - shift,
            edge[1] - shift, 
            Chem.rdchem.BondType.SINGLE
        )
    return mol.GetMol()


def get_bonded_neighbors(mol, atom_idx):
    """Get the bonded neighbors of an atom.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    atom_idx : int
        Index of the atom.

    Returns
    -------
    neighbors : list
        List of bonded neighbors.
    """
    neighbors = []
    for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        neighbors.append(neighbor.GetIdx())
    return neighbors

def find_unique_non_repeating_neighbors(atom_neighbors, idx_to_sym_class):
    """Find all the unique non-repeating neighbors of an atom.
    
    Parameters
    ----------
    atom_neighbors : list
        List of bonded neighbors.
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.
    
    Returns
    -------
    sorted_unique_neighbors_no_repeat : list
        List of unique non-repeating neighbors.
    """
    
    dict_type_repeats = {}
    for neighbor in atom_neighbors:
        
        type_neighbor = idx_to_sym_class[neighbor.GetIdx()]
        
        # Count the number of times a neighbor type is repeated
        if type_neighbor not in dict_type_repeats:
            dict_type_repeats[type_neighbor] = 0
        dict_type_repeats[type_neighbor] += 1

    unique_neighbors = []
    # If the neighbor type is repeated only once, add it to the list
    for neighbor in atom_neighbors:
        type_neighbor = idx_to_sym_class[neighbor.GetIdx()]
        if dict_type_repeats[type_neighbor] == 1:
            unique_neighbors.append(neighbor.GetIdx())
    sorted_unique_neighbors_no_repeat = sorted(unique_neighbors, reverse=True)
    
    return sorted_unique_neighbors_no_repeat

def remove_from_list(atom_list, atom):
    """Remove an atom from a list of atoms.

    Parameters
    ----------
    atom_list : list
        List of atoms.
    atom : rdkit.Chem.Atom
        Atom to remove.
    
    Returns
    -------
    new_atom_list : list
        List of atoms without the atom to remove.
    """

    new_atom_list = []
    for new_atom in atom_list:
        if atom.GetIdx() != new_atom.GetIdx():
            new_atom_list.append(new_atom)
    return new_atom_list

def check_all_atom_same_class(class_list, idx_to_sym_class):
    """Check if all atoms in a list have the same symmetry class.

    Parameters
    ----------
    class_list : list
        List of atoms.
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.
    
    Returns
    -------
    all_symm : bool
        True if all atoms have the same symmetry class.
    """

    all_symm = True
    if len(class_list) >= 1:
        first_class = idx_to_sym_class[class_list[0].GetIdx()]
        for atom in class_list:
            if idx_to_sym_class[atom.GetIdx()] != first_class:
                all_symm = False
                break
    return all_symm

def check_neighbors_same_type(atom, neighbors, idx_to_sym_class):
    """Check if the neighbors of an atom have the same symmetry class.

    Parameters
    ----------
    atom : rdkit.Chem.Atom
        Atom.
    neighbors : list
        List of neighbors.
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.
    
    Returns
    -------
    same_type : bool
        True if the neighbors have the same symmetry class.
    """
    
    same_type = False
    ref_type = idx_to_sym_class[atom.GetIdx()]
    neighbors_type = [idx_to_sym_class[neighbor.GetIdx()] for neighbor in neighbors]
    if ref_type in neighbors_type:
        same_type = True
    return same_type

def at_least_heavy_neighbor(atom):
    """Check if an atom has at least one heavy neighbor.

    Parameters
    ----------
    atom : rdkit.Chem.Atom
        Atom.
    
    Returns
    -------
    found_heavy : bool
        True if the atom has at least one heavy neighbor.
    """

    found_heavy = False
    check_neighbors = atom.GetNeighbors()
    for neighbor in check_neighbors:
        if neighbor.GetAtomicNum() != 1:
            found_heavy = True
            break
    return found_heavy

def grab_index_from_unique_type_number(atom_list, type_num, idx_to_sym_class):
    """Grab the index of atoms that have the same symmetry class.

    Parameters
    ----------
    atom_list : list
        List of atoms.
    type_num : int
        Symmetry class.
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.

    Returns
    -------
    atom_index : list
        List of atom indices.
    """

    idx_list = []
    for atom in atom_list:
        idx = atom.GetIdx()
        if idx_to_sym_class[idx] == type_num:
            if idx not in idx_list:
                idx_list.append(idx)
    return idx_list
        




if __name__ == "__main__":
    
    # import symmetry

    # smiles = "C1=CC=CC=C1"
    # mol = load_molecule(smiles)
    # idx_to_sym_class, symmetry_class = symmetry.get_canonical_labels(mol)
    
    # atom_iter = mol.GetAtoms()
    
    # for atom in atom_iter:
            
    #     atom_index = atom.GetIdx()
    #     print('Atom index:',atom_index)
    #     atom_neighbors = atom.GetNeighbors()
    #     print('Atom neighbors:',[neighbor.GetIdx() for neighbor in atom_neighbors])
    
    #     sorted_unique_neighbors_no_repeat = find_unique_non_repeating_neighbors(atom_neighbors, idx_to_sym_class)

        
    #     print('Unique neighbors:',sorted_unique_neighbors_no_repeat)

    xyz = 'aspirin.xyz'
    mol = load_molecule_from_tinker_xyz(xyz)
    is_single_mol = is_single_molecule(mol)