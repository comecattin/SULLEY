#!/usr/bin/env python3
"""This module is used to write the peditin file

Made by C. Cattin, 2024.
"""

def local_frame_to_output(
        mol,
        idx_to_bisec_then_z_bool, idx_to_trisec_bool,
        idx_to_bisec_idx, idx_to_trisec_idx,
        local_frame1, local_frame2,
    ):

    local_frame1 = sanitize_local_frame(local_frame1)
    local_frame2 = sanitize_local_frame(local_frame2)

    output_local_frame = []

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()

        # Two atoms define the local frame
        if (
            not idx_to_bisec_then_z_bool[atom_idx]
            and not idx_to_trisec_bool[atom_idx]
        ):
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
            output_local_frame.append(
                [
                    atom_idx + 1,
                    -(trisec_idx[0] + 1),
                    -(trisec_idx[1] + 1),
                    -(trisec_idx[2] + 1)
                ]
            )
    return output_local_frame


def write_peditin_file(
        output_local_frame,
        filename
    ):
    """Write the local frame file.

    Parameters
    ----------
    output_local_frame : list
        Local frame for every atom.
    filename : str
        Name of the file.
    """
    with open(filename, 'w') as f:
        for atom in output_local_frame:
            f.write(' '.join([str(i) for i in atom]) + '\n')
    

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
    pass