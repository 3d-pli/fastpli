import numpy as np


def _replace_mat_entry(B, r, d):
    return np.linalg.det(np.delete(np.vstack([r, B]), (d + 1), axis=0))


def nearest_neighbors(i, j, image, M):
    x_max, y_max = image.shape[0] - 1, image.shape[1] - 1
    x, y, _ = M @ np.array([i, j, 1.0])
    x = max(0, min(x_max, int(round(x))))
    y = max(0, min(y_max, int(round(y))))
    return image[x, y]


def calc_back_affine_transformation(p_in, p_out):
    p_in = np.array(p_in)
    p_out = np.array(p_out)

    if not np.all(np.equal(np.array(p_in.shape), np.array(p_out.shape))):
        raise TypeError("in and out not the same shape")

    if not np.all(np.equal(np.array(p_in.shape), np.array([3, 2]))):
        print(p_in.shape)
        raise TypeError("shape error: [3x2], [3x2]")

    l = p_in.shape[0]
    B = np.vstack([np.transpose(p_in), np.ones(l)])
    D = 1.0 / np.linalg.det(B)
    M = np.array([[(-1)**i * D * _replace_mat_entry(B, R, i)
                   for i in range(l)]
                  for R in np.transpose(p_out)])

    return np.vstack([M, [0, 0, 1]])
    # return M[:, :l - 1], M[:, -1]


# # input data
# ins = [[1, 1], [2, 3], [3, 2]]  # <- points
# out = [[0, 2], [1, 2], [-2, -1]]  # <- mapped to
# # calculations

# for i in np.array(ins):
#     print(i)
# print(ins)

# # A, b = calc_back_affine_transformation(ins, out)
# # for i in ins:
# #     print(np.dot(A, i) + b)
