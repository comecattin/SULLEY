from sulley import write_file

def test_shift_multiple_local_frame():
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
    local_frame_shifted = write_file.shift_multiple_local_frame(local_frame)
    local_frame_shifted_ref = [
        [1, 2, 3, 0],
        [2, 3, 1, 0],
        [3, 2, 0, 0],
        [4, 5, 8, 0],
        [5, 8, 0, 0],
        [6, 5, 8, 0],
        [7, 5, 8, 0],
        [8, -9, -10, 0],
        [9, 8, 10, 0],
        [10, 8, 9, 0]
    ]
    assert local_frame_shifted == local_frame_shifted_ref