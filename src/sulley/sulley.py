#!/usr/bin/env python3

""" This module intend to generate a local frame file."""

import symmetry


def generate_local_frame(mol):

    idx_to_bisec_then_z_bool={}
    idx_to_bisec_idx={}
    idx_to_trisec_bool={}
    idx_to_trisec_idx={}

    atom_iter = mol.GetAtoms()

    for atom in atom_iter:

        atom_index = atom.GetIdx()

        idx_to_bisec_then_z_bool[atom_index] = False
        idx_to_trisec_bool[atom_index] = False

        # Get the atom properties
        hybridization = atom.GetHybridization()
        atomic_num = atom.GetAtomicNum()
        atom_neighbors = atom.GetNeighbors()
        valence = len(atom_neighbors)

        # Get the symmetry class of the atom and its neighbors
        idx_to_sym_class, symmetry_class = symmetry.get_canonical_labels(mol)
        neighbors_type = list(
            [idx_to_sym_class[neighbor.GetIdx()] for neighbor in atom_neighbors]
        )
        unique_neighbors_type = list(set(neighbors_type))

        


if __name__ == "__main__":
    pass