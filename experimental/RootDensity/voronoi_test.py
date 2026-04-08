import numpy as np
from scipy.spatial import ConvexHull, Voronoi

points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]])

vor = Voronoi(points)

print(vor)

print("points", points.shape)
print("pioint_region", vor.point_region.shape, np.min(vor.point_region), np.max(vor.point_region))
print(vor.point_region)
