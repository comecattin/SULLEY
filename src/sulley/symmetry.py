#!/usr/bin/env python3
""" This module intend to extract the symmetry of a molecule."""

import networkx as nx
import networkx.algorithms.isomorphism as iso
from collections import defaultdict

def get_canonical_labels(
        mol,
        use_ecfp:bool = False,
        radius:int = None,
        start_idx:int = 0
    ):
    """Get the canonical labels of the atoms in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    use_ecfp : bool, optional
        Use the ECFP method to compute the symmetry classes, by default False
    radius : int, optional
        Radius of the ECFP, by default 3
    start_idx : int, optional
        Starting index for the symmetry classes, by default 0
    
    Returns
    -------
    idx_to_sym_class : dict
        Dictionary of the atom index to the symmetry class.
    symmetry_class : list
        List of symmetry classes.
    """

    if use_ecfp:
        index_to_matching_indices = compute_symmetry_type_ecfp(mol, radius=radius)
    else:
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
            idx_to_sym_class[idx] = sym_class
    
    symmetry_class = idx_to_sym_class.values()

    return idx_to_sym_class, symmetry_class

def create_graph(mol):
    """Create a graph from a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    
    Returns
    -------
    nx.Graph
        NetworkX graph object.
    """

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
    
    return graph

def split_disconnected_graphs(graph):
    """Split a graph into disconnected subgraphs.

    Parameters
    ----------
    graph : nx.Graph
        NetworkX graph object.
    
    Returns
    -------
    list
        List of disconnected subgraphs.
    """
    original_node_order = list(graph.nodes)
    connected = nx.connected_components(graph)
    subgraphs = []

    for nodes in connected:
        sorted_nodes = sorted(nodes, key=original_node_order.index)
        subgraph = nx.Graph()
        for node in sorted_nodes:
            subgraph.add_node(
                node,
                **graph.nodes[node]
            )
        subgraph.add_edges_from(
            (u, v) for u, v in graph.edges(sorted_nodes)
        )
        subgraphs.append(subgraph)

    return subgraphs
    


def compute_symmetry_type(mol):
    """Define the symmetry type of a molecule by a graph invariant vector.
    
    If two atoms have the same graph invariant vector then they have the same type.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    
    Returns
    -------

    index_to_matching_indices : dict
        Dictionary of the atom index to the list of matching indices.

    """

    n_atoms = mol.GetNumAtoms()
    index_to_matching_indices = {i: [i] for i in range(n_atoms)}

    graph = create_graph(mol)
    
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
    
    index_to_matching_indices = {k: sorted(v) for k, v in index_to_matching_indices.items()}

    return index_to_matching_indices

def compute_symmetry_type_ecfp(mol, radius=3):
    """Define the symmetry type of a molecule by a ECFP.
    
    If two atoms have the same ECFP then they have the same type.
    This method is much faster than the graph invariant vector method.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.

    radius : int, optional
        Radius of the ECFP, by default 3.
    
    Returns
    -------

    index_to_matching_indices : dict
        Dictionary of the atom index to the list of matching indices.

    """

    # Get ECFP atom properties
    atom_bonded = []
    species = []
    attached_hydrogens = []
    is_in_ring = []

    for atom in mol.GetAtoms():
        species.append(atom.GetAtomicNum())
        attached_hydrogens.append(atom.GetTotalNumHs())
        is_in_ring.append(atom.IsInRing())
        neighbors_index = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
        atom_bonded.append(neighbors_index)

    # Compute the ECFP
    ecfp = compute_ecfp(species, attached_hydrogens, is_in_ring, atom_bonded, radius=radius)
    ecfp_set = set(ecfp)
    d_ecfp = {k: v for v, k in enumerate(ecfp_set)}
    sym = [d_ecfp[k] for k in ecfp]
    n_atoms = len(species)
    index_to_matching_indices = {}
    for i in range(n_atoms):
        index_to_matching_indices[i] = sorted([j for j in range(n_atoms) if sym[j] == sym[i]])
    
    return index_to_matching_indices
    

def compute_ecfp(species, attached_hydrogens, is_in_ring, neighbors, radius=3):
    """Compute the ECFP invariants of a molecule.

    Parameters
    ----------
    species : list
        List of atomic numbers.
    attached_hydrogens : list
        List of the number of attached hydrogens.
    is_in_ring : list
        List of booleans indicating if the atom is in a ring.
    neighbors : list
        List of list of neighbors.
    radius : int, optional
        Radius of the ECFP, by default 3.
    """
    identifiers = []
    # Iteration 0
    for i in range(len(species)):
        identifier_init = get_hash(
            [
                len(neighbors[i]),
                species[i],
                attached_hydrogens[i],
                is_in_ring[i],
            ]
        )
        identifiers.append(identifier_init)

    # Iteration 1 to radius
    for layer in range(radius):
        new_identifiers = []
        for i in range(len(species)):
            atom_new_identifier = [layer,identifiers[i]]+[identifiers[j] for j in sorted(neighbors[i])]
            new_identifiers.append(get_hash(atom_new_identifier))
        identifiers = new_identifiers

    return identifiers

def get_hash(args):
    """Get the hash of a tuple."""
    return hash(tuple(args))

def symetrize_charges(mol, sym):
    """Symetrize the charges of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object.
    sym : list
        List of symmetry classes.

    Returns
    -------
    qsym : list
        List of symetrized charges.
    formal_charges : list
        List of formal charges.
    """
    
    qsum = defaultdict(float)
    qcount = defaultdict(int)
    formal_charges = [0] * len(sym)
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        qsum[sym[atom.GetIdx()]] += charge
        qcount[sym[atom.GetIdx()]] += 1
        formal_charges[atom.GetIdx()] = charge
    
    qsym = [0] * len(sym)
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        qsym[i] = qsum[sym[i]] / qcount[sym[i]]

    return qsym, formal_charges



if __name__ == "__main__":
    
    from extract_neighbors import load_molecule

    smiles = "C1=CC=CC=C1"
    mol = load_molecule(smiles)
    index_to_matching_indices = compute_symmetry_type(mol)
    idx_to_sym_class, symmetry_class = get_canonical_labels(mol)