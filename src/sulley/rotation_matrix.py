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
    z_only_mask = (
        (local_frame[:, -2] == 0) &
        (local_frame[:, -1] == 0)
    )
    src = np.where(z_only_mask)[0]
    dst = local_frame[src, 1] - 1 # as local frame is 1-indexed
    rot_mat[src] = z_only_rotation(src, dst, coords)


    #----------------------#
    # Z-then-X local frame #
    #----------------------#
    z_then_x_mask = (
        (local_frame[:, -2] > 0) &
        (local_frame[:, -1] == 0)
    )
    src = np.where(z_then_x_mask)[0]
    dst = np.abs(local_frame[src, 1:3]) - 1 # as local frame is 1-indexed
    rot_mat[src] = z_then_x_rotation(src, dst, coords)
    

    return rot_mat

def z_only_rotation(src, dst, coords):
    """Compute the rotation matrix for the z-only local frame.
    
    Parameters
    -----------
    src: np.array
        Source atoms.
    dst: np.array
        Destination atoms.
    coords: np.array
        Coordinates of the atoms.

    Returns
    -------
    rot_mat: np.array
        Rotation matrix.
    """

    vec_z = coords[dst] - coords[src]
    vec_z = vec_z / np.linalg.norm(vec_z)

    # ux = random - (random . uz) uz
    vec_x = np.random.rand(vec_z.shape)
    vec_x = vec_x - np.sum(vec_x * vec_z) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

def z_then_x_rotation(src, dst, coords):
    """Compute the rotation matrix for the z-then-x local frame.
    
    Parameters
    -----------
    src: np.array
        Source atoms.
    dst: np.array
        Destination atoms.
    coords: np.array
        Coordinates of the atoms.

    Returns
    -------
    rot_mat: np.array
        Rotation matrix.
    """
    
    z_dst = dst[:, 0]
    vec_z = coords[z_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_z = np.where(
        np.linalg.norm(vec_z) == 0,
        np.array([0, 0, 1]), 
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z)
    
    # ux = atom1 - (atom1 . uz) uz
    x_dst = dst[:, 1]
    vec_x = coords[x_dst] - coords[src]
    # Remplace the zero vectors by [1, 0, 0]
    vec_x = np.where(
        np.linalg.norm(vec_x) == 0,
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x - np.sum(vec_x * vec_z) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat





if __name__ == "__main__":
    from sulley.extract_neighbors import load_molecule_from_tinker_xyz
    from sulley.local_frame import generate_local_frame 
    xyz = '/home/ccattin/dev/SULLEY/test/poltype/structures/aspirin.xyz'
    mol = load_molecule_from_tinker_xyz(xyz)
    local_frame = generate_local_frame(mol=mol, filename='local_frame.txt')
    rot_mat = compute_rotation_matrix(local_frame, None)