#!/usr/bin/env python3

"""This is the CLI for the Sulley package.

Made by C. Cattin 2024
"""

import argparse
from sulley.extract_neighbors import load_molecule
from sulley.local_frame import generate_local_frame

def main():
    parser = argparse.ArgumentParser(
        description='Generate local frames for a molecule.'
    )
    parser.add_argument(
        'smiles',
        type=str,
        help='The molecule SMILES to generate the local frame for.'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='The file to write the local frame to. Default is local_frame.txt.',
        default='local_frame.txt',
    )

    args = parser.parse_args()

    # Sulley
    mol = load_molecule(args.smiles)
    generate_local_frame(mol=mol, filename=args.output)

    print(f"Local frame written to {args.output}.")

if __name__ == '__main__':
    main()