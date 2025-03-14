#!/usr/bin/env python3

"""This is the CLI for the Sulley package.

Made by C. Cattin 2024
"""

import argparse
from sulley.extract_neighbors import load_molecule, load_molecule_from_sdf, load_molecule_from_tinker_xyz
from sulley.local_frame import generate_local_frame

def main(args=None):

    # Check arguments
    if not args.smiles and not args.sdf and not args.xyz:
        raise ValueError("You must provide a SMILES or SDF or XYZ file.")
    if args.use_ecfp and not args.radius:
        raise ValueError("You must provide a radius for the ECFP method.")


    # Sulley
    if args.smiles:
        mol = load_molecule(args.smiles)
    if args.sdf:
        mol = load_molecule_from_sdf(args.sdf)
    if args.xyz:
        mol = load_molecule_from_tinker_xyz(
            filename=args.xyz,
            return_single_molecule=args.return_single_molecule
        )
        if type(mol) is list:
            for i, m in enumerate(mol):
                print(f"Generating local frame for molecule {i+1}...")
                local_frame = generate_local_frame(
                    mol=m,
                    filename=args.output.replace(".txt", f"_{i+1}.txt"),
                    use_ecfp=args.use_ecfp,
                    radius=args.radius,
                )
                print(f"Local frame written to {args.output.replace('.txt', f'_{i+1}.txt')}.")
            return local_frame

    local_frame = generate_local_frame(
        mol=mol, filename=args.output,
        use_ecfp=args.use_ecfp,
        radius=args.radius
    )

    print(f"Local frame written to {args.output}.")

    if args.verbose:
        with open(args.output, 'r') as f:
            print(f.read())
    
    if args.debug:
        import os
        import sys
        sys.path.append(os.getcwd() + "/test/poltype/")
        from poltype import Poltype
        from multipole import gen_peditinfile
        from test_with_poltype import open_sdf_convert_to_mol

        mol = open_sdf_convert_to_mol(args.sdf, "temp.mol")
        poltype = Poltype(
            mol,
            peditinfile="temp_poltype.txt",
            molstructfname="temp.mol",
            paramfile=os.getcwd()+"/test/poltype/polarize.prm"
        )

        gen_peditinfile(poltype, mol)

        if args.verbose:
            print("-------------")
            with open("temp_poltype.txt", 'r') as f:
                print(f.read())
        

    return local_frame

def cli():
    parser = argparse.ArgumentParser(
        description='Generate local frames for a molecule.'
    )
    parser.add_argument(
        '--smiles',
        type=str,
        help='The molecule SMILES to generate the local frame for.'
    )
    parser.add_argument(
        '--sdf',
        type=str,
        help='The molecule SDF file to generate the local frame for.'
    )
    parser.add_argument(
        '--xyz',
        type=str,
        help='The molecule XYZ Tinker file to generate the local frame for.',
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='The file to write the local frame to. Default is local_frame.txt.',
        default='local_frame.txt',
    )
    parser.add_argument(
        '--use_ecfp',
        action='store_true',
        help='Use ECFP for faster symmetry class computation.'
    )
    parser.add_argument(
        '--radius',
        type=int,
        help='The ECFP radius. Default is 3.',
        default=3,
    )
    parser.add_argument(
        '--return_single_molecule',
        action='store_true',
        help='Return the local frame of a single molecule in the case of a multi-molecule XYZ file. Use with --xyz and --use_ecfp.'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Print the local frame to the console.'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Print debug information. Compare to the Poltype local frame.',
    )

    args = parser.parse_args()

    main(args=args)

    

if __name__ == '__main__':
    pass