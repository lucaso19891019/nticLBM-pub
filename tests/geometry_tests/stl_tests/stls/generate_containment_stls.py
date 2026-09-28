#!/usr/bin/env python3

import math
import os
import struct


OUTPUT_DIR = os.path.join(
    "tests",
    "geometry_tests",
    "stl_tests",
    "stls",
)

SUBDIVISIONS = 2


#=============================================================================
# Vector operations
#=============================================================================

def add(a, b):
    return (
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
    )


def subtract(a, b):
    return (
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2],
    )


def multiply(a, scalar):
    return (
        a[0] * scalar,
        a[1] * scalar,
        a[2] * scalar,
    )


def dot(a, b):
    return (
        a[0] * b[0] +
        a[1] * b[1] +
        a[2] * b[2]
    )


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a):
    return math.sqrt(
        dot(a, a)
    )


def normalize(a):
    length = norm(a)

    if length == 0.0:
        raise RuntimeError(
            "Cannot normalize zero-length vector."
        )

    return (
        a[0] / length,
        a[1] / length,
        a[2] / length,
    )


#=============================================================================
# Icosahedron
#=============================================================================

def create_icosahedron():
    phi = (
        1.0 +
        math.sqrt(5.0)
    ) / 2.0

    vertices = [
        (-1.0,  phi,  0.0),
        ( 1.0,  phi,  0.0),
        (-1.0, -phi,  0.0),
        ( 1.0, -phi,  0.0),

        ( 0.0, -1.0,  phi),
        ( 0.0,  1.0,  phi),
        ( 0.0, -1.0, -phi),
        ( 0.0,  1.0, -phi),

        ( phi,  0.0, -1.0),
        ( phi,  0.0,  1.0),
        (-phi,  0.0, -1.0),
        (-phi,  0.0,  1.0),
    ]

    vertices = [
        normalize(vertex)
        for vertex in vertices
    ]

    faces = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),

        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),

        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),

        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ]

    return vertices, faces


#=============================================================================
# Icosphere subdivision
#=============================================================================

def subdivide(vertices, faces):
    new_vertices = list(vertices)

    midpoint_cache = {}


    def midpoint_index(index0, index1):
        edge = (
            min(index0, index1),
            max(index0, index1),
        )

        if edge in midpoint_cache:
            return midpoint_cache[edge]

        vertex0 = new_vertices[index0]
        vertex1 = new_vertices[index1]

        midpoint = normalize(
            (
                0.5 * (vertex0[0] + vertex1[0]),
                0.5 * (vertex0[1] + vertex1[1]),
                0.5 * (vertex0[2] + vertex1[2]),
            )
        )

        index = len(new_vertices)

        new_vertices.append(
            midpoint
        )

        midpoint_cache[edge] = index

        return index


    new_faces = []


    for index0, index1, index2 in faces:
        midpoint01 = midpoint_index(
            index0,
            index1,
        )

        midpoint12 = midpoint_index(
            index1,
            index2,
        )

        midpoint20 = midpoint_index(
            index2,
            index0,
        )


        new_faces.extend(
            [
                (
                    index0,
                    midpoint01,
                    midpoint20,
                ),
                (
                    index1,
                    midpoint12,
                    midpoint01,
                ),
                (
                    index2,
                    midpoint20,
                    midpoint12,
                ),
                (
                    midpoint01,
                    midpoint12,
                    midpoint20,
                ),
            ]
        )


    return new_vertices, new_faces


#=============================================================================
# Sphere
#=============================================================================

def create_sphere(
    center,
    radius,
    subdivisions=SUBDIVISIONS,
):
    vertices, faces = create_icosahedron()


    for _ in range(subdivisions):
        vertices, faces = subdivide(
            vertices,
            faces,
        )


    vertices = [
        add(
            center,
            multiply(
                vertex,
                radius,
            ),
        )
        for vertex in vertices
    ]


    return vertices, faces


#=============================================================================
# Combine components
#=============================================================================

def combine_components(components):
    vertices = []

    faces = []


    for component_vertices, component_faces in components:
        offset = len(vertices)

        vertices.extend(
            component_vertices
        )


        for face in component_faces:
            faces.append(
                (
                    face[0] + offset,
                    face[1] + offset,
                    face[2] + offset,
                )
            )


    return vertices, faces


#=============================================================================
# STL writer
#=============================================================================

def write_binary_stl(
    filename,
    vertices,
    faces,
):
    header_text = (
        "nticLBM STL containment test"
    )

    header = header_text.encode(
        "ascii"
    )

    header = header.ljust(
        80,
        b"\0",
    )[:80]


    with open(filename, "wb") as file:
        file.write(
            header
        )

        file.write(
            struct.pack(
                "<I",
                len(faces),
            )
        )


        for face in faces:
            vertex0 = vertices[
                face[0]
            ]

            vertex1 = vertices[
                face[1]
            ]

            vertex2 = vertices[
                face[2]
            ]


            edge0 = subtract(
                vertex1,
                vertex0,
            )

            edge1 = subtract(
                vertex2,
                vertex0,
            )


            normal = normalize(
                cross(
                    edge0,
                    edge1,
                )
            )


            file.write(
                struct.pack(
                    "<12fH",
                    normal[0],
                    normal[1],
                    normal[2],

                    vertex0[0],
                    vertex0[1],
                    vertex0[2],

                    vertex1[0],
                    vertex1[1],
                    vertex1[2],

                    vertex2[0],
                    vertex2[1],
                    vertex2[2],

                    0,
                )
            )


#=============================================================================
# Write test case
#=============================================================================

def write_case(
    filename,
    spheres,
):
    components = []


    for center, radius in spheres:
        components.append(
            create_sphere(
                center,
                radius,
            )
        )


    vertices, faces = combine_components(
        components
    )


    path = os.path.join(
        OUTPUT_DIR,
        filename,
    )


    write_binary_stl(
        path,
        vertices,
        faces,
    )


    print(
        f"[GENERATED] {path}"
    )

    print(
        f"            components = {len(spheres)}, "
        f"facets = {len(faces)}"
    )


#=============================================================================
# Main
#=============================================================================

def main():
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )


    #-------------------------------------------------------------------------
    # 1. Single root
    #
    # One isolated component.
    #-------------------------------------------------------------------------

    write_case(
        "single_root.stl",
        [
            (
                (0.0, 0.0, 0.0),
                5.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 2. Two separated roots
    #
    # The two spheres and their bounding boxes are separated.
    #-------------------------------------------------------------------------

    write_case(
        "two_roots.stl",
        [
            (
                (-10.0, 0.0, 0.0),
                3.0,
            ),
            (
                (10.0, 0.0, 0.0),
                3.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 3. AABB overlap without containment
    #
    # Radius = 2 for both spheres.
    #
    # Center distance:
    #
    #     sqrt(3^2 + 3^2) = 4.2426 > 4
    #
    # Therefore the spheres are separated.
    #
    # Their x and y bounding intervals overlap, but neither bounding box
    # contains the other.
    #-------------------------------------------------------------------------

    write_case(
        "aabb_overlap.stl",
        [
            (
                (0.0, 0.0, 0.0),
                2.0,
            ),
            (
                (3.0, 3.0, 0.0),
                2.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 4. True containment
    #
    # Concentric spheres.
    #-------------------------------------------------------------------------

    write_case(
        "nested.stl",
        [
            (
                (0.0, 0.0, 0.0),
                10.0,
            ),
            (
                (0.0, 0.0, 0.0),
                3.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 5. Pseudo containment
    #
    # Outer sphere:
    #
    #     center = (0, 0, 0)
    #     radius = 10
    #
    # Small sphere:
    #
    #     center = (8, 8, 0)
    #     radius = 1
    #
    # Small-sphere AABB:
    #
    #     x = [7, 9]
    #     y = [7, 9]
    #     z = [-1, 1]
    #
    # is completely inside the large-sphere AABB:
    #
    #     [-10, 10]^3
    #
    # However:
    #
    #     center distance = sqrt(8^2 + 8^2)
    #                     = 11.3137
    #
    #     radius sum = 11
    #
    # so the two sphere surfaces are separated and the small sphere is not
    # geometrically contained by the large sphere.
    #-------------------------------------------------------------------------

    write_case(
        "pseudo_containment.stl",
        [
            (
                (0.0, 0.0, 0.0),
                10.0,
            ),
            (
                (8.0, 8.0, 0.0),
                1.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 6. One root with two children
    #
    # Root radius = 10.
    #
    # Two radius-2 children are fully contained by the root and are mutually
    # separated.
    #-------------------------------------------------------------------------

    write_case(
        "two_children.stl",
        [
            (
                (0.0, 0.0, 0.0),
                10.0,
            ),
            (
                (-4.0, 0.0, 0.0),
                2.0,
            ),
            (
                (4.0, 0.0, 0.0),
                2.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 7. Two independent nested roots
    #
    # Tree 1:
    #
    #     sphere at (-15, 0, 0), R = 6
    #         └── sphere at (-15, 0, 0), R = 2
    #
    # Tree 2:
    #
    #     sphere at (15, 0, 0), R = 6
    #         └── sphere at (15, 0, 0), R = 2
    #
    # The two trees are completely separated.
    #-------------------------------------------------------------------------

    write_case(
        "two_nested_roots.stl",
        [
            (
                (-15.0, 0.0, 0.0),
                6.0,
            ),
            (
                (-15.0, 0.0, 0.0),
                2.0,
            ),
            (
                (15.0, 0.0, 0.0),
                6.0,
            ),
            (
                (15.0, 0.0, 0.0),
                2.0,
            ),
        ],
    )


    #-------------------------------------------------------------------------
    # 8. Three containment levels
    #
    # This geometry is valid as an STL, but unsupported by the current
    # Geometry containment model.
    #
    # Expected tree:
    #
    #     R = 10
    #         └── R = 6
    #             └── R = 2
    #
    # analyzeSTLContainment() must reject level 2.
    #-------------------------------------------------------------------------

    write_case(
        "three_levels.stl",
        [
            (
                (0.0, 0.0, 0.0),
                10.0,
            ),
            (
                (0.0, 0.0, 0.0),
                6.0,
            ),
            (
                (0.0, 0.0, 0.0),
                2.0,
            ),
        ],
    )


    print()
    print(
        "All STL containment test geometries generated."
    )


if __name__ == "__main__":
    main()
