#!/usr/bin/env python3

""" This module intend to generate a local frame file."""

import symmetry
import extract_neighbors

def generate_local_frame(mol):

    idx_to_bisec_then_z_bool={}
    idx_to_bisec_idx={}
    idx_to_trisec_bool={}
    idx_to_trisec_idx={}

    local_frame1 = []
    local_frame2 = []

    atom_iter = mol.GetAtoms()

    # Get the symmetry class of the atoms
    idx_to_sym_class, symmetry_class = symmetry.get_canonical_labels(mol)

    for atom in atom_iter:

        is_found_case = False

        atom_index = atom.GetIdx()

        idx_to_bisec_then_z_bool[atom_index] = False
        idx_to_trisec_bool[atom_index] = False

        # Get the atom properties
        hybridization = atom.GetHybridization()
        atomic_num = atom.GetAtomicNum()
        atom_neighbors = atom.GetNeighbors()
        valence = len(atom_neighbors)
        num_hydrogens = atom.GetTotalNumHs()
        is_in_a_ring = atom.IsInRing()

        # Get the neighbors type
        neighbors_type = list(
            [idx_to_sym_class[neighbor.GetIdx()] for neighbor in atom_neighbors]
        )
        unique_neighbors_type = list(set(neighbors_type))
        sorted_unique_neighbors_no_repeat = extract_neighbors.find_unique_non_repeating_neighbors(
            atom_neighbors, idx_to_sym_class
        )
        neighbors_idx = [neighbor.GetIdx() for neighbor in atom_neighbors]

        # Their is at least one unique neighbor
        if len(sorted_unique_neighbors_no_repeat) != 0:
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~#
            #    Atomic environment    #
            #~~~~~~~~~~~~~~~~~~~~~~~~~~#

            # Take the highest symmetry neighbor
            highest_sym_neighbor_idx = sorted_unique_neighbors_no_repeat[0]
            highest_sym_neighbor = mol.GetAtomWithIdx(highest_sym_neighbor_idx)
            
            # Get the properties of the highest symmetry neighbor
            highest_sym_neighbor_atomic_num = highest_sym_neighbor.GetAtomicNum()
            highest_sym_neighbor_hybridization = highest_sym_neighbor.GetHybridization()
            highest_sym_neighbor_num_hydrogens = highest_sym_neighbor.GetTotalNumHs()
            
            # Neighbors of the highest symmetry neighbor
            highest_sym_neighbor_neighbors = highest_sym_neighbor.GetNeighbors()
            highest_sym_neighbor_valence = len(highest_sym_neighbor_neighbors)
            highest_sym_neighbor_neighbors_type = list(
                [idx_to_sym_class[neighbor.GetIdx()]
                for neighbor in highest_sym_neighbor_neighbors]
            )
            highest_sym_neighbor_unique_neighbors_type = list(set(highest_sym_neighbor_neighbors_type))
            highest_sym_neighbor_neighbors_without_atom = extract_neighbors.remove_from_list(highest_sym_neighbor_neighbors, atom)
            highest_sym_neighbor_highest_sym_neighbor = list(
                set(
                    [idx_to_sym_class[neighbor.GetIdx()]
                    for neighbor in highest_sym_neighbor_neighbors_without_atom]
                )
            )

            # Neighbors without the atom
            atom_neighbors_without_atom = extract_neighbors.remove_from_list(atom_neighbors, atom)


            #~~~~~~~~~~~~~~~~~~~~~~~#
            #    Check atom Type    #
            #~~~~~~~~~~~~~~~~~~~~~~~#

            #---------------#
            # Special cases #
            #---------------#

            # Special case for N_2 => Z-only
            if atomic_num == 7 and valence == 1:
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0]
                local_frame2[atom_index] = 0
                is_found_case = True
            
            # Like H on amonia => Z-then-bisec
            elif (
                highest_sym_neighbor_valence == 3
                and extract_neighbors.check_all_atom_same_class(
                    highest_sym_neighbor_neighbors, idx_to_sym_class
                )
                and highest_sym_neighbor_num_hydrogens == 3
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0]
                idx_to_bisec_then_z_bool[atom_index] = True
                bisec_idx = [neighbors.GetIdx()
                            for neighbors in highest_sym_neighbor_neighbors_without_atom]
                idx_to_bisec_idx[atom_index] = bisec_idx
                found_case = True




        

        


if __name__ == "__main__":
    smiles = "C1=CC=CC=C1"
    mol = extract_neighbors.load_molecule(smiles)
    generate_local_frame(mol)