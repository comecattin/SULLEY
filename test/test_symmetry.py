from sulley import symmetry
from sulley import extract_neighbors

def test_symmetry_type_benzene():
    smiles = "C1=CC=CC=C1"
    mol = extract_neighbors.load_molecule(smiles)
    index_to_matching_indices = symmetry.compute_symmetry_type(mol)
    index_to_matching_indices_ref = {
        0: [0, 1, 2, 3, 4, 5],
        1: [1, 0, 2, 3, 4, 5],
        2: [2, 0, 1, 3, 4, 5],
        3: [3, 0, 1, 2, 4, 5],
        4: [4, 0, 1, 2, 3, 5],
        5: [5, 0, 1, 2, 3, 4],
        6: [6, 7, 8, 9, 10, 11],
        7: [7, 6, 8, 9, 10, 11],
        8: [8, 6, 7, 9, 10, 11],
        9: [9, 6, 7, 8, 10, 11],
        10: [10, 6, 7, 8, 9, 11],
        11: [11, 6, 7, 8, 9, 10]
    }
    index_to_matching_indices_ref = {k: sorted(v) for k, v in index_to_matching_indices_ref.items()}
    assert index_to_matching_indices == index_to_matching_indices_ref

def test_symmetry_type_benzene_ecfp():
    smiles = "C1=CC=CC=C1"
    mol = extract_neighbors.load_molecule(smiles)
    index_to_matching_indices = symmetry.compute_symmetry_type_ecfp(mol, radius=3)
    index_to_matching_indices_ref = {
        0: [0, 1, 2, 3, 4, 5],
        1: [1, 0, 2, 3, 4, 5],
        2: [2, 0, 1, 3, 4, 5],
        3: [3, 0, 1, 2, 4, 5],
        4: [4, 0, 1, 2, 3, 5],
        5: [5, 0, 1, 2, 3, 4],
        6: [6, 7, 8, 9, 10, 11],
        7: [7, 6, 8, 9, 10, 11],
        8: [8, 6, 7, 9, 10, 11],
        9: [9, 6, 7, 8, 10, 11],
        10: [10, 6, 7, 8, 9, 11],
        11: [11, 6, 7, 8, 9, 10]
    }
    index_to_matching_indices_ref = {k: sorted(v) for k, v in index_to_matching_indices_ref.items()}
    assert index_to_matching_indices == index_to_matching_indices_ref

def test_get_canonical_labels_benzene():
    smiles = "C1=CC=CC=C1"
    mol = extract_neighbors.load_molecule(smiles)
    (
        idx_to_sym_class,
        symmetry_class
    ) = symmetry.get_canonical_labels(mol)

    assert idx_to_sym_class == {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1}
    assert list(symmetry_class) == [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]