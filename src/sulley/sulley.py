#!/usr/bin/env python3

"""This is the CLI for the Sulley package.

Made by C. Cattin 2024
"""

import argparse
from sulley.extract_neighbors import load_molecule, load_molecule_from_sdf
from sulley.local_frame import generate_local_frame

def main():
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
        '-o',
        '--output',
        type=str,
        help='The file to write the local frame to. Default is local_frame.txt.',
        default='local_frame.txt',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Print the local frame to the console.'
    )

    args = parser.parse_args()

    # Check arguments
    if not args.smiles and not args.sdf:
        raise ValueError("You must provide a SMILES or SDF file.")
    
    if args.smiles and args.sdf:
        raise ValueError("You must provide either a SMILES or SDF file, not both.")


    # Sulley
    if args.smiles:
        mol = load_molecule(args.smiles)
    if args.sdf:
        mol = load_molecule_from_sdf(args.sdf)

    generate_local_frame(mol=mol, filename=args.output)

    print(f"Local frame written to {args.output}.")

    if args.verbose:
        with open(args.output, 'r') as f:
            print(f.read())
    

if __name__ == '__main__':
    main()