#!/usr/bin/env python3
"""This module is used to write the peditin file

Made by C. Cattin, 2024.
"""

def write_peditin_file(
        mol,
        idx_to_bisec_then_z_bool, idx_to_trisec_bool,
        idx_to_bisec_idx, idx_to_trisec_idx,
        local_frame1, local_frame2,
        filename:str
    ):

    f = open(filename, 'w')

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()

        # Two atoms define the local frame
        if (
            not idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
            f.write(
                str(atom_idx) + " " +
                str(local_frame1[atom_idx]) + " " +
                str(local_frame2[atom_idx]) + "\n"
            )
        
        # Bisection then z-axis
        elif (
            idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
            bisec_idx = idx_to_bisec_idx[atom_idx]
            f.write(
                str(atom_idx) + " " +
                str(local_frame1[atom_idx]) + " -" +
                str(bisec_idx[0]) + " -" +
                str(bisec_idx[1]) + "\n"
            )
        
        # Trisection
        else:
            trisec_idx = idx_to_trisec_idx[atom_idx]
            f.write(
                str(atom_idx) + " -" +
                str(trisec_idx[0]) + " -" +
                str(trisec_idx[1]) + " -" +
                str(trisec_idx[2]) + "\n"
            )




if __name__ == "__main__":
    pass