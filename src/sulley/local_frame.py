#!/usr/bin/env python3

""" This module intend to generate a local frame file.

Made by C. Cattin, 2024.
"""

from sulley import symmetry
from sulley import extract_neighbors
from sulley import write_file
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import AllChem
from collections import Counter

def generate_local_frame(mol, filename="local_frame.txt", use_ecfp=False, radius=3):

    # Get the conformer
    mol.UpdatePropertyCache()
    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer(0)

    idx_to_bisec_then_z_bool={}
    idx_to_bisec_idx={}
    idx_to_trisec_bool={}
    idx_to_trisec_idx={}

    local_frame1 = [None] * mol.GetNumAtoms()
    local_frame2 = [None] * mol.GetNumAtoms()

    atom_iter = mol.GetAtoms()

    # Get the symmetry class of the atoms
    idx_to_sym_class, _symmetry_class = symmetry.get_canonical_labels(mol, use_ecfp=use_ecfp, radius=radius)

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
        num_hydrogens = atom.GetTotalNumHs(includeNeighbors=True)
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
            highest_sym_neighbor_num_hydrogens = highest_sym_neighbor.GetTotalNumHs(includeNeighbors=True)
            
            # Neighbors of the highest symmetry neighbor
            highest_sym_neighbor_neighbors = highest_sym_neighbor.GetNeighbors()
            highest_sym_neighbor_valence = len(highest_sym_neighbor_neighbors)
            highest_sym_neighbor_neighbors_type = list(
                [idx_to_sym_class[neighbor.GetIdx()]
                for neighbor in highest_sym_neighbor_neighbors]
            )
            highest_sym_neighbor_unique_neighbors_type = list(set(highest_sym_neighbor_neighbors_type))
            highest_sym_neighbor_neighbors_without_atom = extract_neighbors.remove_from_list(highest_sym_neighbor_neighbors, atom)
            highest_sym_neighbor_unique_neighbors_type_without_atom = list(
                set(
                    [idx_to_sym_class[neighbor.GetIdx()]
                    for neighbor in highest_sym_neighbor_neighbors_without_atom]
                )
            )

            #~~~~~~~~~~~~~~~~~~~~~~~#
            #    Check atom Type    #
            #~~~~~~~~~~~~~~~~~~~~~~~#

            #---------------#
            # Special cases #
            #---------------#

            # Special case for N_2 => Z-only
            if atomic_num == 7 and valence == 1:
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                is_found_case = True
            
            # Like H on amonia => Z-then-bisec
            elif (
                highest_sym_neighbor_valence == 3
                and extract_neighbors.check_all_atom_same_class(
                    highest_sym_neighbor_neighbors, idx_to_sym_class
                )
                and highest_sym_neighbor_num_hydrogens == 3
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                idx_to_bisec_then_z_bool[atom_index] = True
                bisec_idx = [neighbors.GetIdx()
                            for neighbors in highest_sym_neighbor_neighbors_without_atom]
                idx_to_bisec_idx[atom_index] = bisec_idx
                is_found_case = True

            # Like CH3PO3 => Z-only
            elif (
                len(unique_neighbors_type) == 2
                and valence == 4
                and valence == highest_sym_neighbor_valence
                and len(highest_sym_neighbor_unique_neighbors_type) == 2
                and extract_neighbors.check_neighbors_same_type(
                    atom,
                    highest_sym_neighbor_neighbors_without_atom,
                    idx_to_sym_class
                ) is False
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                is_found_case = True

            # Like methyl-amine => Z-then-bisec
            elif (
                (
                    (
                        valence == 4
                        and len(unique_neighbors_type) == 2
                        and highest_sym_neighbor_valence == 3
                        and highest_sym_neighbor_hybridization != 2
                    )
                    or (
                        valence == 3
                        and hybridization != 2
                        and len(unique_neighbors_type) == 2
                        and (
                            highest_sym_neighbor_valence == 3
                            or highest_sym_neighbor_valence == 4
                        )
                    )
                )
                and len(highest_sym_neighbor_unique_neighbors_type) == 2
                and len(highest_sym_neighbor_unique_neighbors_type_without_atom) == 1
                and (
                    num_hydrogens == 2
                    or highest_sym_neighbor_num_hydrogens == 2
                )

            ):
                idx_to_bisec_then_z_bool[atom_index] = True
                if atomic_num == 6:
                    bisec_idx = [
                        neighbors.GetIdx()
                        for neighbors in highest_sym_neighbor_neighbors_without_atom
                    ]
                else:
                    neighbors_without_atom = extract_neighbors.remove_from_list(
                        atom_neighbors, highest_sym_neighbor
                    )
                    bisec_idx = [
                        neighbor.GetIdx()
                        for neighbor in neighbors_without_atom
                    ]
                
                idx_to_bisec_idx[atom_index] = bisec_idx
                local_frame1[atom_index] = highest_sym_neighbor_idx + 1
                is_found_case = True

            # Like ethylamine => Z-then-bisec
            elif (
                valence == 3
                and len(unique_neighbors_type) == 2
                and highest_sym_neighbor_valence == 4
                and hybridization != 2
                and len(highest_sym_neighbor_unique_neighbors_type) == 3
                and len(highest_sym_neighbor_unique_neighbors_type_without_atom) == 2
            ):
                neighbors_without_atom = extract_neighbors.remove_from_list(
                    atom_neighbors, highest_sym_neighbor
                )
                bisec_idx = [
                    neighbor.GetIdx()
                    for neighbor in neighbors_without_atom
                ]
                idx_to_bisec_idx[atom_index] = bisec_idx
                local_frame1[atom_index] = highest_sym_neighbor_idx + 1
                is_found_case = True
                idx_to_bisec_then_z_bool[atom_index] = True

            # Like dimethylamine => Z-then-bisec
            elif(
                valence == 1
                and highest_sym_neighbor_valence == 3
                and highest_sym_neighbor_num_hydrogens == 1
                and len(unique_neighbors_type) == 1
                and len(highest_sym_neighbor_unique_neighbors_type) == 2
                and highest_sym_neighbor_atomic_num == 7
            ):
                idx_to_bisec_then_z_bool[atom_index] = True
                bisec_idx = [
                    neighbors.GetIdx()
                    for neighbors in highest_sym_neighbor_neighbors_without_atom
                ]
                idx_to_bisec_idx[atom_index] = bisec_idx
                local_frame1[atom_index] = highest_sym_neighbor_idx + 1
                is_found_case = True

            # Like dimethylamine => Z-then-bisec
            elif (
                valence == 3
                and highest_sym_neighbor_valence == 1
                and highest_sym_neighbor_num_hydrogens == 0
                and len(unique_neighbors_type) == 2
                and len(highest_sym_neighbor_unique_neighbors_type) == 1
                and atomic_num == 7
            ):
                idx_to_bisec_then_z_bool[atom_index] = True
                neighbors_without_atom = extract_neighbors.remove_from_list(
                    atom_neighbors, highest_sym_neighbor
                )
                bisec_idx = [
                    neighbors.GetIdx()
                    for neighbors in neighbors_without_atom
                ]
                idx_to_bisec_idx[atom_index] = bisec_idx
                local_frame1[atom_index] = highest_sym_neighbor_idx + 1
                is_found_case = True
            
            # Like H on benzene => Z-only
            elif (
                valence == 1
                and highest_sym_neighbor_valence == 3
                and len(unique_neighbors_type) == 1
                and len(highest_sym_neighbor_unique_neighbors_type) == 2
                and len(highest_sym_neighbor_unique_neighbors_type_without_atom) ==1
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                is_found_case = True
            
            # Like H on iodine => Z-only
            elif (
                valence == 1
                and highest_sym_neighbor_valence == 1
                and len(unique_neighbors_type) == 1
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                is_found_case = True

            # Like N(CH3)(CH3)(CH3)H => Z-only
            elif (
                valence == 4
                and len(unique_neighbors_type) == 2
                and (
                    len(highest_sym_neighbor_unique_neighbors_type_without_atom) == 0
                    or len(highest_sym_neighbor_unique_neighbors_type_without_atom) == 1
                )
            ):
                local_frame1[atom_index] = highest_sym_neighbor_idx + 1
                is_found_case = True
        

        #---------------#
        # General cases #
        #---------------#
        if is_found_case is False:
            first_neighbors = atom_neighbors[0].GetNeighbors()

            # Like H in methane => Z-only
            if (
                valence == 1
                and extract_neighbors.check_all_atom_same_class(
                    first_neighbors, idx_to_sym_class
                )
                and len(first_neighbors) == 4
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                is_found_case = True
            
            # Like C in methane => Z-only
            elif (
                valence == 4
                and extract_neighbors.check_all_atom_same_class(
                    atom_neighbors, idx_to_sym_class
                )
                and not extract_neighbors.at_least_heavy_neighbor(atom)
            ):
                idx_list = extract_neighbors.grab_index_from_unique_type_number(atom_neighbors, unique_neighbors_type[0], idx_to_sym_class)
                local_frame1[atom_index] = idx_list[0] + 1
            
            # Like N in ammonia => trisect
            elif (
                valence == 3
                and num_hydrogens == 3
                and extract_neighbors.check_all_atom_same_class(
                    atom_neighbors, idx_to_sym_class
                )
            ):
                idx_to_trisec_bool[atom_index] = True
                trisec_idx = [neighbors.GetIdx() for neighbors in atom_neighbors]
                idx_to_trisec_idx[atom_index] = trisec_idx

            # Like middle C in propane or O in H2O => Bisector
            elif (
                (
                    valence == 2
                    and extract_neighbors.check_all_atom_same_class(
                        atom_neighbors, idx_to_sym_class
                    )
                )
                or (
                    valence == 4
                    and len(unique_neighbors_type) == 2
                )
                and len(sorted_unique_neighbors_no_repeat) == 0
            ):
                idx_list = extract_neighbors.grab_index_from_unique_type_number(atom_neighbors, unique_neighbors_type[0], idx_to_sym_class)
                local_frame1[atom_index] = -1 * (idx_list[0] + 1)
                local_frame2[atom_index] = -1 * (idx_list[1] + 1)
            
            # Like ethane, ethene, ... => Z-only
            elif (
                len(unique_neighbors_type) == 2
                and extract_neighbors.check_neighbors_same_type(
                    atom, atom_neighbors, idx_to_sym_class
                )
                and not is_in_a_ring
            ):
                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
            
            # Like benzene, C on aniline => Bisector
            elif (
                len(unique_neighbors_type) == 2
                and valence == 3
            ):
                type_num_to_count = {}
                for type_num in neighbors_type:
                    count = neighbors_type.count(type_num)
                    type_num_to_count[type_num] = count
                max_count = max(type_num_to_count.values())
                for type_num in unique_neighbors_type:
                    count = type_num_to_count[type_num]
                    if count == max_count:
                        break
                idx_list = extract_neighbors.grab_index_from_unique_type_number(atom_neighbors, type_num, idx_to_sym_class)
                local_frame1[atom_index] = -1 * (idx_list[0] + 1)
                local_frame2[atom_index] = -1 * (idx_list[1] + 1)

            # No special case found, more general processing
            else:
                # Only one unique neighbor: take neighbors of the highest symmetry neighbor
                if len(sorted_unique_neighbors_no_repeat) == 1:
                    neighbors_of_first_neighbor = highest_sym_neighbor.GetNeighbors()
                    neighbors_of_first_neighbor_without_atom = extract_neighbors.remove_from_list(
                        neighbors_of_first_neighbor, atom
                    )
                    sorted_unique_neighbors_no_repeat_new = extract_neighbors.find_unique_non_repeating_neighbors(
                        neighbors_of_first_neighbor_without_atom, idx_to_sym_class
                    )
                    sorted_unique_neighbors_no_repeat += sorted_unique_neighbors_no_repeat_new


                # Check if the molecule is linear
                if len(sorted_unique_neighbors_no_repeat) >= 2:

                    a_atom_idx = atom_index
                    b_atom_idx = sorted_unique_neighbors_no_repeat[0]
                    c_atom_idx = sorted_unique_neighbors_no_repeat[1]
                    indices = [a_atom_idx, b_atom_idx, c_atom_idx]

                    angle = rdMolTransforms.GetAngleDeg(conf,*indices)
                    angle = angle%180
                    angle_tolerance = 3.5
                    is_linear = abs(angle) <= angle_tolerance

                    # The molecule is linear => Z-only
                    if is_linear:
                        local_frame1[atom_index] = sorted_unique_neighbors_no_repeat[0] + 1
                        is_found_case = True
                
                # No unique neighbors found => Z-only
                else:
                    local_frame1[atom_index] = neighbors_idx[0] + 1
                    is_found_case = True
                    continue

                
                # Most general case => Z-then-X
                # Use the first most unique neighbor
                
                sorted_unique_neighbors_no_repeat_types = [
                    idx_to_sym_class[i] for i in sorted_unique_neighbors_no_repeat
                ]
                # Atoms that have the same type as the current atom
                atom_type = idx_to_sym_class[atom_index]
                idx_of_the_same_type_to_move = [
                    index
                    for i, index in enumerate(
                        sorted_unique_neighbors_no_repeat
                    )
                    if sorted_unique_neighbors_no_repeat_types[i] == atom_type
                ]
                # Reorder the list to have the same type at the end
                sorted_unique_neighbors_no_repeat_new = [
                    index for index in sorted_unique_neighbors_no_repeat
                    if index not in idx_of_the_same_type_to_move
                ]
                sorted_unique_neighbors_no_repeat_new.extend(idx_of_the_same_type_to_move)
                # Count the number of each type
                type_num_to_count = Counter(sorted_unique_neighbors_no_repeat_types)
                # Get the type that is the less repeated
                min_count = min(type_num_to_count.values())
                if min_count == 1:
                    sym_type = next(
                        sym_type for sym_type, count in type_num_to_count.items()
                        if count == min_count
                    )
                    special_index = next(
                        index for i, index in enumerate(sorted_unique_neighbors_no_repeat)
                        if sorted_unique_neighbors_no_repeat_types[i] == sym_type
                    )

                    # Reorder the list to have the most rare unique type at the beginning
                    sorted_unique_neighbors_no_repeat_new = [special_index] + [
                        index for index in sorted_unique_neighbors_no_repeat
                        if index != special_index
                    ]
                    sorted_unique_neighbors_no_repeat_new = sorted_unique_neighbors_no_repeat_new[:]

                local_frame1[atom_index] = sorted_unique_neighbors_no_repeat_new[0] + 1
                local_frame2[atom_index] = sorted_unique_neighbors_no_repeat_new[1] + 1
    
    # Write the local frame file
    local_frame = write_file.write_peditin_file(
        mol,
        idx_to_bisec_then_z_bool, idx_to_trisec_bool,
        idx_to_bisec_idx, idx_to_trisec_idx,
        local_frame1, local_frame2,
        filename = filename
    )

    return local_frame




        

        


if __name__ == "__main__":
    molecules = {
        "Benzene": "c1ccccc1",
        "Lactic Acid": "CC(C(=O)O)O",
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ethanol": "CCO",
        "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "Cholesterol": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCCC4)C)O",
        "Styrene": "C=CC1=CC=CC=C1",
        "Acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
        "Serotonin": "C1=CC2=C(C=C1CCN)NC=C2O",
        "Ammonia": "N",
        "Methane": "C",
        "Phosphine": "P",
        "Arsine": "[AsH3]",
        "Trimethylamine": "N(C)(C)(C)"
    }
    mol = extract_neighbors.load_molecule(molecules["Trimethylamine"])
    generate_local_frame(mol)