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
    """Write the local frame file.

    Parameters
    ----------
    mol : mol : rdkit.Chem.Mol
        RDKit molecule object.
    idx_to_bisec_then_z_bool : dict
        Index associated to a bisection then z-axis.
    idx_to_trisec_bool : dict
        Index associated to a trisection.
    idx_to_bisec_idx : dict
        Index associated to the bisection index.
    idx_to_trisec_idx : dict
        Index associated to the trisection index.
    local_frame1 : list
        Local frame first vector.
    local_frame2 : list
        Local frame second vector.
    filename : str
        Name of the file.
    """

    local_frame1 = sanitize_local_frame(local_frame1)
    local_frame2 = sanitize_local_frame(local_frame2)

    f = open(filename, 'w')

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()

        # Two atoms define the local frame
        if (
            not idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
            f.write(
                str(atom_idx + 1) + " " +
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
                str(atom_idx + 1) + " " +
                str(local_frame1[atom_idx]) + " -" +
                str(bisec_idx[0] + 1) + " -" +
                str(bisec_idx[1] + 1) + "\n"
            )
        
        # Trisection
        else:
            trisec_idx = idx_to_trisec_idx[atom_idx]
            print(trisec_idx)
            f.write(
                str(atom_idx + 1) + " -" +
                str(trisec_idx[0] + 1) + " -" +
                str(trisec_idx[1] + 1) + " -" +
                str(trisec_idx[2] + 1) + "\n"
            )

def sanitize_local_frame(local_frame):
    """Shift the local frame indices to start at 1.
    
    Parameters
    ----------
    local_frame : list
        Local frame.
    
    Returns
    -------
    local_frame : list
        Local frame shifted.
    """
    for i, atom in enumerate(local_frame):
        if atom is None:
            local_frame[i] = 0
        # elif atom < 0:
        #     local_frame[i] -= 1
        # else:
        #     local_frame[i] += 1
    return local_frame 
            



if __name__ == "__main__":
    pass