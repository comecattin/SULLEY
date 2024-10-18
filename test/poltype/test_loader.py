from poltype import Poltype
from multipole import gen_peditinfile
from sulley.local_frame import generate_local_frame
from sulley.extract_neighbors import load_molecule_from_tinker_xyz
from test_with_poltype import open_sdf_convert_to_mol, is_same_local_frame
import os


def test_aspirin_xyz():

    test_path = os.getcwd() + "/test/poltype/"
    
    # Poltype
    mol = open_sdf_convert_to_mol(
        test_path + "structures/aspirin.sdf",
        test_path + "structures/aspirin.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/aspirin_lf_poltype.txt",
        molstructfname=test_path + "structures/aspirin.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)
    
    # Sulley
    mol = load_molecule_from_tinker_xyz(test_path + "structures/aspirin.xyz")
    generate_local_frame(
        mol=mol,
        filename=test_path + "output/aspirin_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/aspirin_lf_sulley.txt",
        test_path + "output/aspirin_lf_poltype.txt"
    )

def test_aspirin_ethanol_multiple_xyz():

    test_path = os.getcwd() + "/test/poltype/"
    
    # Poltype
    # Aspirin
    mol = open_sdf_convert_to_mol(
        test_path + "structures/aspirin.sdf",
        test_path + "structures/aspirin.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/aspirin_lf_poltype.txt",
        molstructfname=test_path + "structures/aspirin.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)

    # Ethanol
    mol = open_sdf_convert_to_mol(
        test_path + "structures/ethanol.sdf",
        test_path + "structures/ethanol.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/ethanol_lf_poltype.txt",
        molstructfname=test_path + "structures/ethanol.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)
    
    # Sulley
    mols = load_molecule_from_tinker_xyz(test_path + "structures/aspirin_ethanol.xyz")
    # Aspirin
    generate_local_frame(
        mol=mols[0],
        filename=test_path + "output/aspirin_lf_sulley.txt"
    )
    # Ethanol
    generate_local_frame(
        mol=mols[1],
        filename=test_path + "output/ethanol_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/aspirin_lf_sulley.txt",
        test_path + "output/aspirin_lf_poltype.txt"
    ) and is_same_local_frame(
        test_path + "output/ethanol_lf_sulley.txt",
        test_path + "output/ethanol_lf_poltype.txt"
    )
