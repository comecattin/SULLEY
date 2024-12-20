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

    #----------------------#
    # Bisector local frame #
    #----------------------#
    bisector_mask = (
        (local_frame[:, -2] < 0) &
        (local_frame[:, -1] == 0)
    )
    src = np.where(bisector_mask)[0]
    dst = np.abs(local_frame[src, 1:3]) - 1 # as local frame is 1-indexed
    rot_mat[src] = bisector_rotation(src, dst, coords)

    #------------------------#
    # Z-Bisector local frame #
    #------------------------#
    z_bisector_mask = (
        (local_frame[:, -3] > 0) &
        (local_frame[:, -2] < 0) &
        (local_frame[:, -1] < 0)
    )
    src = np.where(z_bisector_mask)[0]
    dst = np.abs(local_frame[src, 1:4]) - 1 # as local_frame is 1-indexed
    rot_mat[src] = z_then_bisector_rotation(src, dst, coords)

    #-----------------------#
    # Trisector local frame #
    #-----------------------#
    trisector_mask = (
        (local_frame[:, -3] < 0) &
        (local_frame[:, -2] < 0) &
        (local_frame[:, -1] < 0)
    )
    src = np.where(trisector_mask)[0]
    dst = np.abs(local_frame[src, 1:4]) - 1 # as local_frame is 1-indexed
    rot_mat[src] = trisector_rotation(src, dst, coords)

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

    if dst.shape[0] == 0:
        return np.eye(3)

    vec_z = coords[dst] - coords[src]
    vec_z = vec_z / np.linalg.norm(vec_z, axis=-1, keepdims=True)

    # ux = random - (random . uz) uz
    vec_x = np.random.rand(vec_z.shape[0], vec_z.shape[1])
    vec_x = vec_x - np.sum(vec_x * vec_z, axis=-1, keepdims=True) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)

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

    if dst.shape[0] == 0:
        return np.eye(3)
    
    z_dst = dst[:, 0]
    vec_z = coords[z_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_z = np.where(
        np.linalg.norm(vec_z, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]), 
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z, axis=-1, keepdims=True)
    
    # ux = atom1 - (atom1 . uz) uz
    x_dst = dst[:, 1]
    vec_x = coords[x_dst] - coords[src]
    # Remplace the zero vectors by [1, 0, 0]
    vec_x = np.where(
        np.linalg.norm(vec_x, axis=-1, keepdims=True) == 0,
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x - np.sum(vec_x * vec_z, axis=-1, keepdims=True) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

def bisector_rotation(src, dst, coords):
    """Compute the rotation matrix for the bisector local frame.
        
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

    if dst.shape[0] == 0:
        return np.eye(3)

    # uz1 = atom0
    z_dst = dst[:, 0]
    vect_1 = coords[z_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vect_1 = np.where(
        np.linalg.norm(vect_1, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]),
        vect_1
    )
    vect_1 = vect_1 / np.linalg.norm(vect_1, axis=-1, keepdims=True)
    
    # ux1 = atom1
    x_dst = dst[:, 1]
    vect_2 = coords[x_dst] - coords[src]
    # Replace the zero vectors by [1, 0, 0]
    vect_2 = np.where(
        np.linalg.norm(vect_2, axis=-1, keepdims=True) == 0,
        np.array([1, 0, 0]),
        vect_2
    )
    vect_2 = vect_2 / np.linalg.norm(vect_2, axis=-1, keepdims=True)

    # uz = uz1 + ux1
    vec_z = vect_1 + vect_2
    # Replace non-physical vectors by [0, 0, 1]
    vec_z = np.where(
        (vec_z == np.array([1, 0, 1])).all(axis=1, keepdims=True),
        np.array([0, 0, 1]),
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z, axis=-1, keepdims=True)

    # ux = ux1 - (ux1 . uz) uz
    vec_x = vect_2 - np.sum(vect_2 * vec_z, axis=-1, keepdims=True) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

def z_then_bisector_rotation(src, dst, coords):
    """Compute the rotation matrix for the Z-then-bisector local frame.
        
    Parameters
    -----------
    src: jnp.array
        Source atoms.
    dst: jnp.array
        Destination atoms.
    coords: jnp.array
        Coordinates of the atoms.

    Returns
    -------
    rot_mat: jnp.array
        Rotation matrix.
    """

    if dst.shape[0] == 0:
        return np.eye(3)

    # uz = atom0
    z_dst = dst[:, 0]
    vec_z = coords[z_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_z = np.where(
        np.linalg.norm(vec_z) == 0,
        np.array([0, 0, 1]),
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z)

    # First bisector vector
    bisec1_dst = dst[:, 1]
    vec_bisec1 = coords[bisec1_dst] - coords[src]
    # Replace the zero vectors by [1, 0, 0]
    vec_bisec1 = np.where(
        np.linalg.norm(vec_bisec1) == 0,
        np.array([1, 0, 0]),
        vec_bisec1
    )
    vec_bisec1 = vec_bisec1 / np.linalg.norm(vec_bisec1)

    # Second bisector vector
    bisec2_dst = dst[:, 2]
    vec_bisec2 = coords[bisec2_dst] - coords[src]
    # Replace the zero vectors by [1, 0, 0]
    vec_bisec2 = np.where(
        np.linalg.norm(vec_bisec2) == 0,
        np.array([1, 0, 0]),
        vec_bisec2
    )
    vec_bisec2 = vec_bisec2 / np.linalg.norm(vec_bisec2)

    # ux = bisec1 + bisec2 and then orthogonalized
    vec_x = vec_bisec1 + vec_bisec2
    # Replace non-physical vectors by [1, 0, 0]
    vec_x = np.where(
        (vec_x == np.array([2, 0, 0])).all(),
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x / np.linalg.norm(vec_x)
    vec_x = vec_x - np.sum(vec_x * vec_z) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

def trisector_rotation(src, dst, coords):
    """Compute the rotation matrix for the trisector local frame.
    
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

    if dst.shape[0] == 0:
        return np.eye(3)

    # Trisector vector 1
    trisec1_dst = dst[:, 0]
    vec_trisec1 = coords[trisec1_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_trisec1 = np.where(
        np.linalg.norm(vec_trisec1) == 0,
        np.array([0, 0, 1]),
        vec_trisec1
    )
    vec_trisec1 = vec_trisec1 / np.linalg.norm(vec_trisec1)

    # Trisector vector 2
    trisec2_dst = dst[:, 1]
    vec_trisec2 = coords[trisec2_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_trisec2 = np.where(
        np.linalg.norm(vec_trisec2) == 0,
        np.array([0, 0, 1]),
        vec_trisec2
    )
    vec_trisec2 = vec_trisec2 / np.linalg.norm(vec_trisec2)

    # Trisector vector 3
    trisec3_dst = dst[:, 2]
    vec_trisec3 = coords[trisec3_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_trisec3 = np.where(
        np.linalg.norm(vec_trisec3) == 0,
        np.array([0, 0, 1]),
        vec_trisec3
    )
    vec_trisec3 = vec_trisec3 / np.linalg.norm(vec_trisec3)

    # uz = trisec1 + trisec2 + trisec3
    vec_z = vec_trisec1 + vec_trisec2 + vec_trisec3
    # Replace non-physical vectors by [0, 0, 1]
    vec_z = np.where(
        (vec_z == np.array([0, 0, 3])).all(),
        np.array([0, 0, 1]),
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z)

    # ux = trisec2 - (trisec2 . uz) uz
    vec_x = vec_trisec2 - np.sum(vec_trisec2 * vec_z) * vec_z
    # Replace non-physical vectors by [1, 0, 0]
    vec_x = np.where(
        np.linalg.norm(vec_x) == 0,
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x / np.linalg.norm(vec_x)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

if __name__ == "__main__":
    # from sulley.extract_neighbors import load_molecule_from_tinker_xyz
    # from sulley.local_frame import generate_local_frame 
    # xyz = '/home/ccattin/dev/SULLEY/test/poltype/structures/aspirin.xyz'
    # mol = load_molecule_from_tinker_xyz(xyz)
    # local_frame = generate_local_frame(mol=mol, filename='local_frame.txt')
    # rot_mat = compute_rotation_matrix(local_frame, None)

    # Water test for bisector and z-then-x
    local_frame = [
        [1, -2, -3],
        [2, 1, 3],
        [3, 1, 2]
    ]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
    water = compute_rotation_matrix(local_frame, coords)
    water_ref = np.array(
        [
            [
                [-1, 0, 1],
                [1, 0, 1],
                [0, 1, 0]
            ],
            [
                [0, 0, -1],
                [1, 0, 0],
                [0, -1, 0]
            ],
            [
                [1, 0, 0],
                [0, 0, -1],
                [0, 1, 0]
            ]
        ],
        dtype=np.float32
    )
    water_ref[0,0:2,:] *= 1/np.sqrt(2)

    # HCl test for z-only
    local_frame = [
            [1, 2, 0],
            [2, 1, 0],
        ]
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    hcl = compute_rotation_matrix(local_frame, coordinates)
    # No test available for Z-only as the rotation matrix is random
    # Can just check if the last column is the vector between the two atoms

    
