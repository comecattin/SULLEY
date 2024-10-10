#!/usr/bin/env python3
""" This module intend to extract the symmetry of a molecule."""

from rdkit import Chem
from extract_neighbors import load_molecule
import networkx as nx
import networkx.algorithms.isomorphism as iso

def get_canonical_labels(mol, start_idx:int = 0):
    """Get the canonical labels of the atoms in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    start_idx : int, optional
        Starting index for the symmetry classes, by default 0
    
    Returns
    -------
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.
    symmetry_class : list
        List of symmetry classes.
    """

    index_to_matching_indices = compute_symmetry_type(mol)

    groups = []
    group_to_heavy_atom = {}
    for index, matching_indices in index_to_matching_indices.items():
        atom = mol.GetAtomWithIdx(index)
        atomic_num = atom.GetAtomicNum()
        heavy = False

        # Group atoms if they're heavy or not
        if atomic_num > 1:
            heavy = True
        if set(matching_indices) not in groups:
            groups.append(set(matching_indices))
            group_to_heavy_atom[tuple(set(matching_indices))] = heavy
    
    # Sort the groups
    sorted_groups = []
    for group, heavy in group_to_heavy_atom.items():
        if heavy:
            if group not in sorted_groups:
                sorted_groups.append(group)
    for group, heavy in group_to_heavy_atom.items():
        if not heavy:
            if group not in sorted_groups:
                sorted_groups.append(group)

    sym_class_to_group = {}
    idx_to_sym_class = {}
    symclass = start_idx

    # Assign a symmetry class to each atom
    for group in sorted_groups:
        sym_class_to_group[symclass] = group
        symclass += 1
    for sym_class, group in sym_class_to_group.items():
        for idx in group:
            idx_to_sym_class[idx+1] = sym_class
    
    symmetry_class = idx_to_sym_class.values()

    return idx_to_sym_class, symmetry_class

def compute_symmetry_type(mol):
    """Define the symmetry type of a molecule by a graph invariant vector.
    
    If two atoms have the same graph invariant vector then they have the same type.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    
    Returns
    -------

    """

    n_atoms = mol.GetNumAtoms()
    index_to_matching_indices = {i: [i] for i in range(n_atoms)}

    # Create the graph
    graph = nx.Graph()

    # Add nodes
    for atom in mol.GetAtoms():
        graph.add_nodes_from(
            [(atom.GetIdx(), {"atomic_num": atom.GetAtomicNum()})]
        )

    # Add edges
    for bond in mol.GetBonds():
        graph.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
        )
    
    # Get nodes attributes
    atomic_nums = nx.get_node_attributes(graph, "atomic_num")
    atom_match = iso.numerical_node_match("atomic_num", 0)
    center = nx.center(graph)

    # Compute the node fingerprint
    # The node fingerprint is a list of 3 elements:
    #     - atomic number
    #     - number of neighbors
    #     - distance to the center
    node_fingerprint = {}
    for node in graph.nodes():
        dist_to_center = min(
            [nx.shortest_path_length(
                graph,
                source=node,
                target=center_node
            ) for center_node in center])
        node_fingerprint[node] = [
            atomic_nums[node],
            graph.degree[node],
            dist_to_center
        ]
    
    node_list = sorted(list(graph.nodes()))

    assert node_list == list(range(0, len(node_list)))

    # Create graph list without each node
    graph_list = []
    for node in node_list:
        graph_copy = graph.copy()
        graph_copy.remove_node(node)
        graph_list.append(graph_copy)

    # Compare each subgraph
    for i1 in range(len(node_list)):
        node1 = node_list[i1]
        for i2 in range(i1+1, len(node_list)):
            node2 = node_list[i2]

            # The node fingerprint must be the same
            if node_fingerprint[node1] != node_fingerprint[node2]:
                continue
            # The graph must be isomorphic
            if not nx.faster_could_be_isomorphic(
                graph_list[i1],
                graph_list[i2]
            ):
                continue
            
            # Check if the two graphs are isomorphic
            is_isomorphic = nx.is_isomorphic(
                graph_list[i1],
                graph_list[i2],
                node_match=atom_match
            )
            if is_isomorphic:
                if node2 not in index_to_matching_indices[node1]:
                    index_to_matching_indices[node1].append(node2)
                if node1 not in index_to_matching_indices[node2]:
                    index_to_matching_indices[node2].append(node1)
    
    return index_to_matching_indices



if __name__ == "__main__":
    smiles = "C1=CC=CC=C1"
    mol = load_molecule(smiles)
    index_to_matching_indices = compute_symmetry_type(mol)