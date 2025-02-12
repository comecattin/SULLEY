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
        Name of the file. If None, no file is written.
    
    Returns
    -------
    output_local_frame : list
        Local frame for every atom.
    """

    local_frame1 = sanitize_local_frame(local_frame1)
    local_frame2 = sanitize_local_frame(local_frame2)

    output_local_frame = []

    if filename is not None:
        f = open(filename, 'w')

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()

        # Two atoms define the local frame
        if (
            not idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
            if filename is not None:
                f.write(
                    str(atom_idx + 1) + " " +
                    str(local_frame1[atom_idx]) + " " +
                    str(local_frame2[atom_idx]) + "\n"
                )
            output_local_frame.append(
                [
                    atom_idx + 1,
                    local_frame1[atom_idx],
                    local_frame2[atom_idx]
                ]
            )
        
        # Bisection then z-axis
        elif (
            idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
            bisec_idx = idx_to_bisec_idx[atom_idx]
            if filename is not None:
                f.write(
                    str(atom_idx + 1) + " " +
                    str(local_frame1[atom_idx]) + " -" +
                    str(bisec_idx[0] + 1) + " -" +
                    str(bisec_idx[1] + 1) + "\n"
                )
            output_local_frame.append(
                [
                    atom_idx + 1,
                    local_frame1[atom_idx],
                    -(bisec_idx[0] + 1),
                    -(bisec_idx[1] + 1)
                ]
            )
        
        # Trisection
        else:
            trisec_idx = idx_to_trisec_idx[atom_idx]
            if filename is not None:
                f.write(
                    str(atom_idx + 1) + " -" +
                    str(trisec_idx[0] + 1) + " -" +
                    str(trisec_idx[1] + 1) + " -" +
                    str(trisec_idx[2] + 1) + "\n"
                )
            output_local_frame.append(
                [
                    atom_idx + 1,
                    -(trisec_idx[0] + 1),
                    -(trisec_idx[1] + 1),
                    -(trisec_idx[2] + 1)
                ]
            )
    
    if filename is not None:
        f.close()
    return output_local_frame

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

def shift_multiple_local_frame(local_frame):
    """
    Shift the local frame if the output local frame is from multiple molecules.
    
    Parameters
    ----------
    local_frame : list
        Local frame.
    
    Returns
    -------
    local_frame_shifted : list
        Local frame shifted
    """
    shift = 0
    local_frame_shifted = []
    for i, lf in enumerate(local_frame):
        if lf[0] == 1 and i != 0:
            shift = i
        shifted = [
            idx + shift if idx > 0 else
            idx - shift if idx < 0 else
            0
            for idx in lf
        ]
        local_frame_shifted.append(shifted)
    return local_frame_shifted
            



if __name__ == "__main__":
    local_frame = [
        [1, 2, 3, 0],
        [2, 3, 1, 0],
        [3, 2, 0, 0],
        [1, 2, 5, 0],
        [2, 5, 0, 0],
        [3, 2, 5, 0],
        [4, 2, 5, 0],
        [5, -6, -7, 0],
        [6, 5, 7, 0],
        [7, 5, 6, 0]
    ]
    local_frame_shifted = shift_multiple_local_frame(local_frame)
    print(local_frame_shifted)