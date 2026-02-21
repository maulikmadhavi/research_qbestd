"""
This custom code for Vector quantiztion followed by 4 cluster k-means clustering.
"""

import numpy as np
from utils.gen_data import generated_simulated_data
import einops

x, y = generated_simulated_data()  # 2D array of shape (2, 4000)
N_CLUSTERS = 1
N_SAMPLES = 4000
FEAT_DIM = 2
EPS = 0.000001
x = x.T  # Input features 4000 x 2 (Number of samples=4000), feature dim=2
N_CLUSTERS_REQ = 4


MAX_ITER = 100  # Limit for iterations to prevent infinite loops
CONVERGENCE_THRESHOLD = 1e-4  # Threshold to check if centers have moved significantly


c = np.zeros([1, 2])  # Initial code vector
assignments = np.zeros([N_SAMPLES])


def closeness(x, c):
    """
    Compute the closeness between x and o.
    x: [N_SAMPLES, FEAT_DIM]
    c: [N_CLUSTERS, FEAT_DIM]
    """
    NS = x.shape[0]
    NC = c.shape[0]
    # [c1,c1,c1... N_sample, c2,c2, ... N_sample, c3,c3, ... N_sample, c4,c4, ... N_sample]
    c1 = einops.repeat(c, "k d -> (k n) d", n=NS)
    # [x1,x2,x3,...,xN (1st time), x1, x2, x3 .... xN (2nd time), ..... x1, x2, x3 .... xN (4th time)]
    x1 = einops.repeat(x, "k d -> (n k) d", n=NC)  # x1[:3],x1[4000:4003],x1[8000:8003],x1[12000:12003]

    xc2 = np.sqrt((x1**2).sum(1) - 2 * (x1 * c1).sum(1) + (c1**2).sum(1))
    xc2 = xc2.reshape(NC, NS)

    # [
    # dist(x1,c1), dist(x2,c1), dist(x3,c1), dist(x4,c1), .... dist(xN, c1),
    # dist(x1,c2), dist(x2,c2), dist(x3,c2), dist(x4,c2), .... dist(xN, c2),
    # dist(x1,c3), dist(x2,c3), dist(x3,c3), dist(x4,c3), .... dist(xN, c3),
    # ...
    # ...
    # dist(x1,cM), dist(x2,cM), dist(x3,cM), dist(x4,cM), .... dist(xN, cM),
    #    ]
    return xc2.argmin(axis=0), xc2  # N_CLUSTERS x N_SAMPLES


def new_center_compute(x, assignments, NC):
    """
    Compute new centers.
    x: [N_SAMPLES, FEAT_DIM]
    assignments: [N_SAMPLES]
    """
    for i in range(NC):
        c[i] = x[assignments == i].mean(0)

    # Sort the vector

    return c


def not_change_detect(assignment_prev, assignment_new):
    """
    Detect if the assignment has no changed.
    """
    print(np.mean(assignment_prev == assignment_new))
    return np.mean(assignment_prev == assignment_new) > 0.95


# sourcery skip: avoid-builtin-shadow
iter = 0
while iter < MAX_ITER:
    assignment_prev = assignments.copy()
    # Compute the closeness between x and o
    assignments, _ = closeness(x, c)
    # Compute new centers
    c = new_center_compute(x, assignments, N_CLUSTERS)
    print(f"New centers from assignment: {c}, iter: {iter}")
    iter += 1
    # Detect if the assignment has changed
    if not_change_detect(assignment_prev, assignments):
        if N_CLUSTERS == N_CLUSTERS_REQ:
            print("Converged")
            break
        N_CLUSTERS = N_CLUSTERS * 2
        print(f"Number of clusters increased to {N_CLUSTERS}")
        c = np.concatenate((c - EPS, c + EPS), axis=0)
        # Reset the iter counter
        iter = 0
        print(f"New centers from split: {c}, iter: {iter}")
    if iter == MAX_ITER:
        print(f" iter={iter}: Max iterations reached to prevent infinite loops")
        break

# Compute distances by broadcasting (no heavy tiling)
dists = np.sum((x[:, None, :] - c[None, :, :]) ** 2, axis=2)  # shape (NS, NC) squared distances
assignments = dists.argmin(axis=1)
expected_assignments = y
print(f"Final centers: {c},\n Final assignments: {assignments}")

# Check if how many final assignments match the expected assignments
n_matches = np.sum(assignments == expected_assignments)
print(f"Number of matches: {n_matches}/{N_SAMPLES}")


# Robust recompute of centers handling empty clusters
def recompute_centers(x, assignments, NC):
    centers = np.zeros((NC, x.shape[1]))
    for i in range(NC):
        pts = x[assignments == i]
        if pts.size == 0:
            # reinitialize empty cluster (example: small random offset from global mean)
            centers[i] = x.mean(0) + np.random.randn(x.shape[1]) * 1e-2
        else:
            centers[i] = pts.mean(0)
    return centers
