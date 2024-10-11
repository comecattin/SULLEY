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
    
    unique_neighbors_type = []
    unique_neighbors = []
    for neighbor in atom_neighbors:
        neighbor_index = neighbor.GetIdx()
        neighbor_sym_class = idx_to_sym_class[neighbor_index]
        if neighbor_sym_class not in unique_neighbors_type:
            unique_neighbors.append(neighbor_index)
            unique_neighbors_type.append(neighbor_sym_class)
    sorted_unique_neighbors_no_repeat = sorted(unique_neighbors, reverse=True)
    return sorted_unique_neighbors_no_repeat




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