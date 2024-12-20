#!/usr/bin/python
from sulley.rotation_matrix import compute_rotation_matrix
import numpy as np

def test_H2O():
    """Water test for bisector and z-then-x"""
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
    assert np.allclose(water, water_ref)

def test_hcl():
    """HCl test for z-only
    
    No test available for Z-only as the rotation matrix is random
    Can just check if the last column is the vector between the two atoms
    """
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
    assert np.allclose(
        hcl[:,:,2], np.array([[0, 0, 1], [0, 0, -1]])
    )