#!/usr/bin/env python3
"""This module is used to extract neighbors of an atom up to the 4th neighbor.

Made by C. Cattin, 2024.
"""

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

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
    mol = Chem.Mol(mol)
    return mol

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
        if type_neighbor not in dict_type_repeats.keys():
            dict_type_repeats[type_neighbor] = 0
        dict_type_repeats[type_neighbor] += 1

        unique_neighbors = []
        # If the neighbor type is repeated only once, add it to the list
        for type_neighbor, repeats in dict_type_repeats.items():
            if repeats == 1:
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
        if atom.getIdx() != new_atom.getIdx():
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
        




if __name__ == "__main__":
    
    import symmetry

    smiles = "C1=CC=CC=C1"
    mol = load_molecule(smiles)
    idx_to_sym_class, symmetry_class = symmetry.get_canonical_labels(mol)
    
    atom_iter = mol.GetAtoms()
    
    for atom in atom_iter:
            
        atom_index = atom.GetIdx()
        print('Atom index:',atom_index)
        atom_neighbors = atom.GetNeighbors()
        print('Atom neighbors:',[neighbor.GetIdx() for neighbor in atom_neighbors])
    
        sorted_unique_neighbors_no_repeat = find_unique_non_repeating_neighbors(atom_neighbors, idx_to_sym_class)

        
        print('Unique neighbors:',sorted_unique_neighbors_no_repeat)