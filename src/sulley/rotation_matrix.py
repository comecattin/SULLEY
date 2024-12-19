#!/usr/bin/python

"""Rotation matrix from cartesian frame to local frame.

Made by C. Cattin, 2024.
"""

import numpy as np

def compute_rotation_matrix(local_frame, coords):
    """Compute the rotation matrix from the cartesian frame to the local frame."""

    rot_mat = np.zeros((len(local_frame), 3, 3))

    # Making sure the local frame has the correct shape
    for lf in local_frame:
        if len(lf) == 3:
            lf.append(0)
    local_frame = np.array(local_frame)

    #--------------------#
    # Z-only local frame #
    #--------------------#
    

    return rot_mat





if __name__ == "__main__":
    from sulley.extract_neighbors import load_molecule_from_tinker_xyz
    from sulley.local_frame import generate_local_frame 
    xyz = '/home/ccattin/dev/SULLEY/test/poltype/structures/aspirin.xyz'
    mol = load_molecule_from_tinker_xyz(xyz)
    local_frame = generate_local_frame(mol=mol, filename='local_frame.txt')
    rot_mat = compute_rotation_matrix(local_frame, None)