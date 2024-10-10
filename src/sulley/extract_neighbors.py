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




if __name__ == "__main__":
    smiles = "C1=CC=CC=C1"
    mol = load_molecule(smiles)
    neighbors = get_bonded_neighbors(mol, 0)