import math
import numpy as np
import matplotlib.pyplot as plt
from difflib import get_close_matches

def is_counterclockwise(polygon_points):
    area = 0
    for i in range(len(polygon_points)):
        x1, y1 = polygon_points[i]
        x2, y2 = polygon_points[(i + 1) % len(polygon_points)]
        area += (x2 - x1) * (y2 + y1)
    return area < 0  # Negative means CCW

def is_touching(edge1, polygon2, tolerance=1e-1):
    """
    Determines if an edge (line segment) is touching any edge in another polygon.
    This ensures that the two segments actually overlap, not just share a single point.
    """
    def line_segments_match(seg1, seg2, tol):
        """Checks if two line segments overlap within a tolerance."""
        (x1, y1), (x2, y2) = seg1
        (x3, y3), (x4, y4) = seg2

        # Sort endpoints to make comparisons consistent
        seg1 = sorted(seg1)
        seg2 = sorted(seg2)

        # Check if both endpoints are approximately equal
        return (
            (math.isclose(seg1[0][0], seg2[0][0], abs_tol=tol) and
             math.isclose(seg1[0][1], seg2[0][1], abs_tol=tol) and
             math.isclose(seg1[1][0], seg2[1][0], abs_tol=tol) and
             math.isclose(seg1[1][1], seg2[1][1], abs_tol=tol))
        )

    # Get edges from polygon2
    polygon_edges = [(polygon2[i], polygon2[(i+1) % len(polygon2)]) for i in range(len(polygon2))]

    # Check if any edge in polygon2 matches edge1
    for other_edge in polygon_edges:
        if line_segments_match(edge1, other_edge, tolerance):
            return True

    return False

def plot_polygon_debug(polygon_points, touching_edges):
    """
    Plots the polygon being checked. 
    - Touching edges are drawn in green.
    - Non-touching edges are drawn in red.
    """
    fig, ax = plt.subplots()

    # Loop through each edge of the polygon
    for side in range(len(polygon_points)):
        p1 = polygon_points[side]
        p2 = polygon_points[(side + 1) % len(polygon_points)]  # Wrap around

        # Determine edge color
        color = "green" if (p1, p2) in touching_edges or (p2, p1) in touching_edges else "red"

        # Plot the edge
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=2)

    # Scatter plot the vertices
    x_vals, y_vals = zip(*polygon_points)
    ax.scatter(x_vals, y_vals, color="blue", label="Vertices", zorder=3)

    # Formatting
    ax.set_aspect("equal")
    ax.set_title("Polygon Debug Visualization")
    ax.legend()
    plt.show()

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

def add_vertices_to_edges(polygons, tolerance=1e-1):
    """
    Ensures that every vertex that lies on another polygon’s edge is added as a vertex to that edge.
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

    return new_polygons

def find_touching_edges(polygons_element, polygon_id, polygon_points):
    touching_edges = {}

    # Check for touching polygons
    for other_polygon_element in polygons_element:  # Loop over polygons_element instead of polygons
        other_polygon_id = int(other_polygon_element.attrib.get('ID'))
        
        if polygon_id != other_polygon_id:
            other_polygon_points = []
            for point_element in other_polygon_element.findall("Point"):
                x = float(point_element.attrib["x"])
                y = float(point_element.attrib["y"])
                other_polygon_points.append((x, y))

            other_material_name = other_polygon_element.attrib.get('Material')

            for side in range(len(polygon_points)):  # Loop through each side
                p1 = polygon_points[side]
                p2 = polygon_points[(side + 1) % len(polygon_points)]

                if is_touching([p1, p2], other_polygon_points):
                    touching_edges[(p1, p2)] = other_material_name  # Store touching material
    return touching_edges

def find_closest_material(layer_name, materials):
    # Find the closest matching material name in the CSV file
    material_names = [material['Name'] for material in materials]

    closest_matches = get_close_matches(layer_name.replace("_", " ").lower().strip(), material_names, n=1, cutoff=0.6)
    if closest_matches:
        closest_name = closest_matches[0]
        for material in materials:
            if material['Name'] == closest_name:
                return material
    return None