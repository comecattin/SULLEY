# SULLEY 🚪 ⚡
<p align='center'>
  <img src="https://github.com/user-attachments/assets/d2298efa-b9c7-43e1-b3c3-857f15c4cd98" alt="drawing SULLEY" width="400"/>
</p>



**SULLEY** stands for **Symmetry Understanding of Local Frames for Learning Equivariant geometrY**. This project focuses on developing a framework to find local frames for each atom in a molecule, allowing for the accurate description of their electrostatic multipoles.

## Objective

The primary goal of SULLEY is to identify the local frames of each atom, facilitating the description of their multipole characteristics. By leveraging these local frames, the project aims to enhance the learning of multipoles in an equivariant manner.

## Methodology

SULLEY draws inspiration from existing methodologies, particularly the Poltype framework, utilizing a possibility tree approach to model the local frames effectively. This will enable a structured exploration of atomic environments and their interactions.

## Key Features

- **Equivariant Learning**: The local frames will be used to train neural networks in a way that respects the symmetries of the molecular system.
- **Efficiency**: By optimizing the representation of local frames, SULLEY aims to create a much faster neural network for learning multipoles, significantly improving computational performance.

## Applications

SULLEY's approach is intended for applications in fields such as computational chemistry and molecular dynamics, where understanding electrostatic interactions at the atomic level is crucial.

## Installation

To install SULLEY:

1. Clone current repository

   ```bash
   git clone git@github.com:comecattin/SULLEY.git .
   ```

2. Navigate to the cloned repository

   ```bash
   cd SULLEY
   ```

3. Install application

   ```bash
   pip install -e .
   ```

## Using Sulley

To generate local frames for a molecule, use the following arguments:

### Arguments

- `--smiles` : The molecule SMILES to generate the local frame for.
- `--sdf` : The molecule SDF file to generate the local frame for.
- `--xyz` : The molecule Tinker XYZ file to generate the local frame for. Multiple molecules can be given in the same XYZ file.
- `-o`, `--output` : The file to write the local frame to. Default is `local_frame.txt`.
- `-v`, `--verbose` : Print the local frame to the console.
- `--debug` : Print debug information. Compare to the Poltype local frame.
- `-h`, `--help`: Display the CLI help.

### Examples

1. Generate a local frame from a SMILES and write it to the default file:

    ```bash
    sulley --smiles "CCO"
    ```

2. Generate a local frame from an SDF file and specify an output file:

    ```bash
    sulley --sdf molecule.sdf -o output_frame.txt
    ```

## Project Status

The project is currently under development and ongoing research. Contributions and suggestions are welcome!
