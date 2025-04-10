import numpy as np

def prepare(polygons):
    '''
    This is where we will store logic for cleaning and preparing the dxf files. This is going to include things like adding vertices to edges, snapping close lines, etc.
    '''
    add_vertices_to_edges(polygons)

def add_vertices_to_edges(polygons, tolerance=1e-1):
    """
    Ensures that every vertex that lies on another polygon’s edge is added as a vertex to that edge. This is necessary to draw boundary conditions properly later.
    """
    new_polygons = [polygon for polygon in polygons]  # Deep copy to avoid modifying input

    for poly_idx, polygon in enumerate(polygons):
        for vertex in polygon[0]:
            for other_idx, other_polygon in enumerate(polygons):
                if poly_idx == other_idx:
                    continue  # Don't check against itself

                for edge_idx in range(len(other_polygon[0])):  # Only iterate over the coordinates
                    edge_start = other_polygon[0][edge_idx]
                    edge_end = other_polygon[0][(edge_idx + 1) % len(other_polygon[0])]

                    if point_on_segment(vertex, edge_start, edge_end, tolerance):
                        # Insert the vertex if it isn't already present
                        if not np.allclose(vertex, edge_start, atol=tolerance) and \
                           not np.allclose(vertex, edge_end, atol=tolerance):
                            new_polygons[other_idx][0].insert(edge_idx + 1, vertex)

def point_on_segment(point, seg_start, seg_end, tolerance=1e-1):
    """
    Checks if a point lies on a line segment within a given tolerance.
    """
    p = np.array(point, dtype=float)
    a = np.array(seg_start, dtype=float)
    b = np.array(seg_end, dtype=float)
    
    # Vector cross-product should be near zero (collinear)
    cross_product = np.cross(b - a, p - a)
    if abs(cross_product) > tolerance:
        return False

    # Check if the point is within segment bounds
    dot_product = np.dot(p - a, b - a)
    if dot_product < 0 or dot_product > np.dot(b - a, b - a):
        return False
    
    return True
