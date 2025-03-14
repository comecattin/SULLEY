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

    # uz = atom0
    z_dst = dst[:, 0]
    vec_z = coords[z_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_z = np.where(
        np.linalg.norm(vec_z, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]),
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z, axis=-1, keepdims=True)

    # First bisector vector
    bisec1_dst = dst[:, 1]
    vec_bisec1 = coords[bisec1_dst] - coords[src]
    # Replace the zero vectors by [1, 0, 0]
    vec_bisec1 = np.where(
        np.linalg.norm(vec_bisec1, axis=-1, keepdims=True) == 0,
        np.array([1, 0, 0]),
        vec_bisec1
    )
    vec_bisec1 = vec_bisec1 / np.linalg.norm(vec_bisec1, axis=-1, keepdims=True)

    # Second bisector vector
    bisec2_dst = dst[:, 2]
    vec_bisec2 = coords[bisec2_dst] - coords[src]
    # Replace the zero vectors by [1, 0, 0]
    vec_bisec2 = np.where(
        np.linalg.norm(vec_bisec2, axis=-1, keepdims=True) == 0,
        np.array([1, 0, 0]),
        vec_bisec2
    )
    vec_bisec2 = vec_bisec2 / np.linalg.norm(vec_bisec2, axis=-1, keepdims=True)

    # ux = bisec1 + bisec2 and then orthogonalized
    vec_x = vec_bisec1 + vec_bisec2
    # Replace non-physical vectors by [1, 0, 0]
    vec_x = np.where(
        (vec_x == np.array([2, 0, 0])).all(axis=1, keepdims=True),
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)
    vec_x = vec_x - np.sum(vec_x * vec_z, axis=-1, keepdims=True) * vec_z
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)

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
        np.linalg.norm(vec_trisec1, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]),
        vec_trisec1
    )
    vec_trisec1 = vec_trisec1 / np.linalg.norm(vec_trisec1, axis=-1, keepdims=True)

    # Trisector vector 2
    trisec2_dst = dst[:, 1]
    vec_trisec2 = coords[trisec2_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_trisec2 = np.where(
        np.linalg.norm(vec_trisec2, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]),
        vec_trisec2
    )
    vec_trisec2 = vec_trisec2 / np.linalg.norm(vec_trisec2, axis=-1, keepdims=True)

    # Trisector vector 3
    trisec3_dst = dst[:, 2]
    vec_trisec3 = coords[trisec3_dst] - coords[src]
    # Replace the zero vectors by [0, 0, 1]
    vec_trisec3 = np.where(
        np.linalg.norm(vec_trisec3, axis=-1, keepdims=True) == 0,
        np.array([0, 0, 1]),
        vec_trisec3
    )
    vec_trisec3 = vec_trisec3 / np.linalg.norm(vec_trisec3, axis=-1, keepdims=True)

    # uz = trisec1 + trisec2 + trisec3
    vec_z = vec_trisec1 + vec_trisec2 + vec_trisec3
    # Replace non-physical vectors by [0, 0, 1]
    vec_z = np.where(
        (vec_z == np.array([0, 0, 3])).all(axis=1, keepdims=True),
        np.array([0, 0, 1]),
        vec_z
    )
    vec_z = vec_z / np.linalg.norm(vec_z, axis=-1, keepdims=True)

    # ux = trisec2 - (trisec2 . uz) uz
    vec_x = vec_trisec2 - np.sum(vec_trisec2 * vec_z, axis=-1, keepdims=True) * vec_z
    # Replace non-physical vectors by [1, 0, 0]
    vec_x = np.where(
        np.linalg.norm(vec_x, axis=-1, keepdims=True) == 0,
        np.array([1, 0, 0]),
        vec_x
    )
    vec_x = vec_x / np.linalg.norm(vec_x, axis=-1, keepdims=True)

    # uy = uz x ux
    vec_y = np.cross(vec_z, vec_x)

    rot_mat = np.stack([vec_x, vec_y, vec_z], axis=-1)

    return rot_mat

if __name__ == "__main__":
    pass