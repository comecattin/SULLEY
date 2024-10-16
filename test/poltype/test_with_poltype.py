#!/usr/bin/env python3
from poltype import Poltype
from openbabel import openbabel
from multipole import gen_peditinfile
from sulley.local_frame import generate_local_frame
from sulley.extract_neighbors import load_molecule, load_molecule_from_sdf
import os

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




def test_trimethylamine():

    test_path = os.getcwd() + "/test/poltype/"
    
    # Poltype
    mol = open_sdf_convert_to_mol(
        test_path + "structures/trimethylamine.sdf",
        test_path + "structures/trimethylamine.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/trimethylamine_lf_poltype.txt",
        molstructfname=test_path + "structures/trimethylamine.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)
    
    # Sulley
    mol = load_molecule(molecules["Trimethylamine"])
    generate_local_frame(
        mol=mol,
        filename=test_path + "output/trimethylamine_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/trimethylamine_lf_sulley.txt",
        test_path + "output/trimethylamine_lf_poltype.txt"
    )

def test_benzene():

    test_path = os.getcwd() + "/test/poltype/"
    
    # Poltype
    mol = open_sdf_convert_to_mol(
        test_path + "structures/benzene.sdf",
        test_path + "structures/benzene.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/benzene_lf_poltype.txt",
        molstructfname=test_path + "structures/benzene.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)
    
    # Sulley
    mol = load_molecule_from_sdf(test_path + "structures/benzene.sdf")
    generate_local_frame(
        mol=mol,
        filename=test_path + "output/benzene_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/benzene_lf_sulley.txt",
        test_path + "output/benzene_lf_poltype.txt"
    )

def test_lactic():

    test_path = os.getcwd() + "/test/poltype/"
    
    # Poltype
    mol = open_sdf_convert_to_mol(
        test_path + "structures/lactic.sdf",
        test_path + "structures/lactic.mol"
    )
    poltype = Poltype(
        mol,
        peditinfile=test_path + "output/lactic_lf_poltype.txt",
        molstructfname=test_path + "structures/lactic.mol",
        paramfile=test_path + "polarize.prm"
    )

    gen_peditinfile(poltype, mol)
    
    # Sulley
    mol = load_molecule_from_sdf(test_path + "structures/lactic.sdf")
    generate_local_frame(
        mol=mol,
        filename=test_path + "output/lactic_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/lactic_lf_sulley.txt",
        test_path + "output/lactic_lf_poltype.txt"
    )

def test_aspirin():

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
    mol = load_molecule_from_sdf(test_path + "structures/aspirin.sdf")
    generate_local_frame(
        mol=mol,
        filename=test_path + "output/aspirin_lf_sulley.txt"
    )

    assert is_same_local_frame(
        test_path + "output/aspirin_lf_sulley.txt",
        test_path + "output/aspirin_lf_poltype.txt"
    )


def open_sdf_convert_to_mol(sdf_file, mol_file):
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, sdf_file)
    mol.AddHydrogens()

    poltype = Poltype(mol)
    obConversion.SetOutFormat("mol")
    obConversion.WriteFile(mol, mol_file)

    return mol

def is_same_local_frame(sully_file, poltype_file):
    
    with open(sully_file, "r") as sully:
        sully_lines = sully.readlines()
    with open(poltype_file, "r") as poltype:
        poltype_lines = poltype.readlines()

    for sully_line, poltype_line in zip(sully_lines, poltype_lines):
        if sully_line != poltype_line:
            return False

    return True