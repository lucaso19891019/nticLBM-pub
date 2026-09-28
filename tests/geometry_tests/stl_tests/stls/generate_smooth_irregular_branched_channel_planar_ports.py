import os
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.ndimage import gaussian_filter
from skimage import measure
import trimesh
import matplotlib.pyplot as plt

# ============================================================
# Fixed smooth irregular branched channel generator
#
# Output:
#   smooth_irregular_branched_channel.stl
#   smooth_irregular_branched_channel_preview.png
#
# Units: mm
#
# Dependencies:
#   pip install numpy scipy scikit-image trimesh matplotlib
#
# Run:
#   python generate_smooth_irregular_branched_channel.py
# ============================================================


OUTPUT_STL = "smooth_irregular_branched_channel.stl"
OUTPUT_PREVIEW = "smooth_irregular_branched_channel_preview.png"

# Marching-cubes grid spacing.
# Smaller -> smoother/finer STL, but more memory/time.
GRID_SPACING = 1.15

# Light smoothing of the implicit field.
SMOOTH_SIGMA = 0.65

# Number of spline samples along each centerline.
CENTERLINE_SAMPLES = 220

# Length of the axis-aligned straight section added before each exposed port.
# It only affects the terminal portion of the tube and keeps the bifurcations smooth.
PORT_STRAIGHT_LENGTH = 12.0

# Small extra extension beyond the final cutting plane.  Marching cubes first
# creates a closed rounded end; we then slice that rounded end away and cap it
# exactly with a planar face.
PORT_CUT_EXTENSION = 4.0


# ------------------------------------------------------------
# Geometry definition
# ------------------------------------------------------------
#
# Each branch is a smooth variable-radius tube.
# The full fluid domain is the union of all tubes.
#
# "pts" are spline control points.
# r0/r1 are branch start/end radii.
#
# Branches overlap near bifurcations so the implicit union creates
# smooth junctions without explicit CAD boolean operations.
# ------------------------------------------------------------

BRANCHES = [
    # Main inlet trunk
    {
        "pts": np.array([
            [-75, -10,  -2],
            [-58,  -7,   6],
            [-40,   2,  12],
            [-22,  15,   6],
            [ -5,  20,  -4],
            [ 12,  12, -10],
            [ 28,   0,  -6],
        ], dtype=float),
        "r0": 9.5,
        "r1": 8.0,
    },

    # Major upper daughter
    {
        "pts": np.array([
            [ 24,   2,  -6],
            [ 34,  10,   0],
            [ 44,  24,   8],
            [ 52,  39,  14],
            [ 66,  48,  10],
            [ 79,  53,   2],
        ], dtype=float),
        "r0": 7.6,
        "r1": 5.3,
    },

    # Major lower daughter
    {
        "pts": np.array([
            [ 24,  -1,  -6],
            [ 34, -12,  -9],
            [ 42, -27, -14],
            [ 55, -39, -10],
            [ 69, -46,   0],
            [ 82, -50,   9],
        ], dtype=float),
        "r0": 7.4,
        "r1": 5.0,
    },

    # Upper daughter, out-of-plane branch
    {
        "pts": np.array([
            [ 47,  30,  10],
            [ 50,  34,  19],
            [ 54,  37,  30],
            [ 61,  34,  39],
            [ 70,  29,  44],
        ], dtype=float),
        "r0": 5.2,
        "r1": 3.6,
    },

    # Upper daughter, lateral secondary branch
    {
        "pts": np.array([
            [ 57,  42,  13],
            [ 63,  47,  18],
            [ 69,  49,  26],
            [ 76,  46,  34],
            [ 84,  41,  39],
        ], dtype=float),
        "r0": 4.8,
        "r1": 3.2,
    },

    # Lower daughter, downward out-of-plane branch
    {
        "pts": np.array([
            [ 49, -33, -12],
            [ 54, -39, -21],
            [ 58, -40, -31],
            [ 65, -35, -40],
            [ 75, -28, -45],
        ], dtype=float),
        "r0": 5.0,
        "r1": 3.4,
    },

    # Lower daughter, lateral secondary branch
    {
        "pts": np.array([
            [ 60, -43,  -5],
            [ 66, -47,  -9],
            [ 72, -47, -17],
            [ 79, -43, -25],
            [ 87, -37, -30],
        ], dtype=float),
        "r0": 4.5,
        "r1": 3.0,
    },
]


def build_smooth_centerline(points, n_samples):
    """
    Fit a cubic parametric spline through control points and sample it densely.

    Chord-length parameterization gives better behavior than uniformly assigning
    parameters to unevenly spaced control points.
    """
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    u = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    u /= u[-1]

    k = min(3, len(points) - 1)

    # Small positive smoothing avoids tiny local wiggles while remaining close
    # to the specified geometry.
    tck, _ = splprep(points.T, u=u, s=0.25, k=k)

    uu = np.linspace(0.0, 1.0, n_samples)
    xyz = np.vstack(splev(uu, tck)).T

    return xyz, uu


def build_radius_profile(u, r0, r1):
    """
    Smoothly vary radius from r0 to r1.

    A cosine easing avoids slope jumps at branch endpoints.
    A small smooth modulation makes the channel slightly irregular without
    introducing abrupt diameter changes.
    """
    w = 0.5 - 0.5 * np.cos(np.pi * u)
    radius = r0 * (1.0 - w) + r1 * w

    # Mild smooth irregularity. Zero at both endpoints.
    modulation = (
        1.0
        + 0.045
        * np.sin(2.0 * np.pi * u + 0.7)
        * np.sin(np.pi * u)
    )

    return radius * modulation


def _dominant_axis_direction(tangent):
    """
    Snap a tangent to the closest Cartesian axis (+/-X, +/-Y, +/-Z).

    The sign is preserved so the terminal tube continues outward rather than
    turning back into the geometry.
    """
    tangent = np.asarray(tangent, dtype=float)
    axis = int(np.argmax(np.abs(tangent)))

    direction = np.zeros(3, dtype=float)
    direction[axis] = 1.0 if tangent[axis] >= 0.0 else -1.0

    return direction


def _add_axis_aligned_terminal(centerline, radius, at_start=False):
    """
    Replace only the exposed terminal portion by a short Cartesian-axis-aligned
    straight section.

    Returns
    -------
    new_centerline, new_radius, port
        port is a dict containing the exact plane point and outward normal used
        later to cut/cap the marching-cubes mesh.
    """
    centerline = np.asarray(centerline, dtype=float)
    radius = np.asarray(radius, dtype=float)

    if at_start:
        endpoint = centerline[0].copy()

        # Outward tangent points from interior toward the exposed start.
        tangent = centerline[0] - centerline[min(8, len(centerline) - 1)]
        direction = _dominant_axis_direction(tangent)

        # Keep a point inside the original tube, then transition to an
        # axis-aligned terminal section.
        anchor_idx = min(
            len(centerline) - 2,
            max(2, int(0.10 * len(centerline)))
        )
        anchor = centerline[anchor_idx].copy()

        # The final port plane lies PORT_STRAIGHT_LENGTH away from the anchor
        # along a Cartesian axis. Add a little material beyond it so slicing
        # produces a clean cap.
        port_point = anchor + direction * PORT_STRAIGHT_LENGTH
        extended_end = port_point + direction * PORT_CUT_EXTENSION

        n_terminal = max(
            10,
            int((PORT_STRAIGHT_LENGTH + PORT_CUT_EXTENSION) / GRID_SPACING * 3)
        )
        terminal = np.linspace(extended_end, anchor, n_terminal, endpoint=False)

        r_port = radius[0]
        r_anchor = radius[anchor_idx]
        rr_terminal = np.linspace(r_port, r_anchor, n_terminal, endpoint=False)

        new_centerline = np.vstack([terminal, centerline[anchor_idx:]])
        new_radius = np.concatenate([rr_terminal, radius[anchor_idx:]])

    else:
        endpoint = centerline[-1].copy()

        # Outward tangent points from interior toward the exposed end.
        tangent = centerline[-1] - centerline[max(0, len(centerline) - 9)]
        direction = _dominant_axis_direction(tangent)

        anchor_idx = max(
            1,
            min(len(centerline) - 3, int(0.90 * len(centerline)))
        )
        anchor = centerline[anchor_idx].copy()

        port_point = anchor + direction * PORT_STRAIGHT_LENGTH
        extended_end = port_point + direction * PORT_CUT_EXTENSION

        n_terminal = max(
            10,
            int((PORT_STRAIGHT_LENGTH + PORT_CUT_EXTENSION) / GRID_SPACING * 3)
        )
        terminal = np.linspace(anchor, extended_end, n_terminal + 1)[1:]

        r_anchor = radius[anchor_idx]
        r_port = radius[-1]
        rr_terminal = np.linspace(r_anchor, r_port, n_terminal + 1)[1:]

        new_centerline = np.vstack([centerline[:anchor_idx + 1], terminal])
        new_radius = np.concatenate([radius[:anchor_idx + 1], rr_terminal])

    axis = int(np.argmax(np.abs(direction)))
    plane_name = ("X", "Y", "Z")[axis]

    port = {
        "point": port_point,
        "normal": direction,
        "plane": plane_name,
        "coordinate": float(port_point[axis]),
    }

    return new_centerline, new_radius, port


def prepare_branches():
    sampled = []
    ports = []

    for ib, branch in enumerate(BRANCHES):
        centerline, u = build_smooth_centerline(
            branch["pts"],
            CENTERLINE_SAMPLES
        )

        radius = build_radius_profile(
            u,
            branch["r0"],
            branch["r1"]
        )

        # Only exposed terminals are planar:
        #   branch 0 start = main inlet
        #   branches 1..6 end = outlets
        # Branch starts inside bifurcations remain untouched.
        if ib == 0:
            centerline, radius, port = _add_axis_aligned_terminal(
                centerline,
                radius,
                at_start=True
            )
            ports.append(port)
        else:
            centerline, radius, port = _add_axis_aligned_terminal(
                centerline,
                radius,
                at_start=False
            )
            ports.append(port)

        sampled.append((centerline, radius))

    return sampled, ports


def compute_bounds(sampled_branches):
    all_points = np.vstack([centerline for centerline, _ in sampled_branches])
    all_radii = np.concatenate([radius for _, radius in sampled_branches])

    margin = float(all_radii.max()) + 7.0

    mins = all_points.min(axis=0) - margin
    maxs = all_points.max(axis=0) + margin

    return mins, maxs


def build_implicit_field(sampled_branches, ports, mins, maxs):
    """
    Construct an implicit union of variable-radius tubes.

    For each branch sample:
        phi = radius - distance_to_centerline_sample

    Positive phi -> inside fluid domain
    Negative phi -> outside

    The union is obtained with max(phi).
    """
    h = GRID_SPACING

    xs = np.arange(mins[0], maxs[0] + h, h)
    ys = np.arange(mins[1], maxs[1] + h, h)
    zs = np.arange(mins[2], maxs[2] + h, h)

    nx, ny, nz = len(xs), len(ys), len(zs)

    print("Grid size:")
    print(f"  nx = {nx}")
    print(f"  ny = {ny}")
    print(f"  nz = {nz}")
    print(f"  voxels = {nx * ny * nz:,}")

    field = np.full(
        (nx, ny, nz),
        -1.0e9,
        dtype=np.float32
    )

    X, Y = np.meshgrid(xs, ys, indexing="ij")

    # Use every second centerline sample for field evaluation.
    # This is still substantially finer than the Cartesian grid spacing
    # while keeping runtime and memory moderate.
    evaluation_sets = []

    for ib, (centerline, radius) in enumerate(sampled_branches):
        idx = np.arange(0, len(centerline), 2)
        evaluation_sets.append((centerline[idx], radius[idx], ports[ib]))

    print("Building implicit field...")

    for k, z in enumerate(zs):
        slab = np.full(
            (nx, ny),
            -1.0e9,
            dtype=np.float32
        )

        for centerline, radius, port in evaluation_sets:

            # First build the implicit field of this branch alone.
            branch_phi = np.full(
                (nx, ny),
                -1.0e9,
                dtype=np.float32
            )

            # Chunk the centerline to limit temporary-array memory.
            for j0 in range(0, len(centerline), 32):
                c = centerline[j0:j0 + 32]
                r = radius[j0:j0 + 32]

                dx = X[..., None] - c[:, 0]
                dy = Y[..., None] - c[:, 1]
                dz = z - c[:, 2]

                distance = np.sqrt(
                    dx * dx
                    + dy * dy
                    + dz * dz
                )

                local_phi = r[None, None, :] - distance

                branch_phi = np.maximum(
                    branch_phi,
                    local_phi.max(axis=2).astype(np.float32)
                )

            # Intersect this terminal branch with the inward half-space.
            # cap_phi > 0 on the kept/interior side and cap_phi = 0 exactly
            # on an X/Y/Z-aligned port plane.
            origin = np.asarray(port["point"], dtype=float)
            outward = np.asarray(port["normal"], dtype=float)

            cap_phi = -(
                (X - origin[0]) * outward[0]
                + (Y - origin[1]) * outward[1]
                + (z - origin[2]) * outward[2]
            )

            branch_phi = np.minimum(
                branch_phi,
                cap_phi.astype(np.float32)
            )

            # Union of all already-capped branches.
            slab = np.maximum(slab, branch_phi)

        field[:, :, k] = slab

        if (k + 1) % 10 == 0 or k == nz - 1:
            print(f"  z-slice {k + 1}/{nz}")

    print("Smoothing implicit field...")

    field = gaussian_filter(
        field,
        sigma=SMOOTH_SIGMA
    )

    return field


def extract_surface(field, mins):
    """
    Extract phi=0 isosurface with marching cubes.
    """
    print("Running marching cubes...")

    vertices, faces, normals, values = measure.marching_cubes(
        field,
        level=0.0,
        spacing=(
            GRID_SPACING,
            GRID_SPACING,
            GRID_SPACING
        )
    )

    vertices += mins

    mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=faces,
        process=True
    )

    # If numerical discretization somehow creates multiple disconnected pieces,
    # retain the largest connected component.
    components = mesh.split(only_watertight=False)

    if len(components) > 1:
        mesh = max(
            components,
            key=lambda m: m.volume if m.is_volume else m.area
        )

    mesh.remove_unreferenced_vertices()

    trimesh.repair.fix_normals(mesh)

    try:
        trimesh.repair.fill_holes(mesh)
    except Exception:
        pass

    return mesh


def _cap_cut_boundary(mesh, plane_origin, outward_normal, tol=1.0e-6):
    """
    Cap boundary loops which lie on one known cutting plane.

    This avoids trimesh's optional polygon-triangulation dependencies.  Port
    sections are close to circular/star-shaped, so a centroid fan is robust
    here and keeps every cap vertex exactly on the requested Cartesian plane.
    """
    vertices = mesh.vertices
    faces = mesh.faces

    # Boundary edges are mesh edges used by exactly one triangle.
    edges = np.sort(
        np.vstack([
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ]),
        axis=1
    )

    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = unique_edges[counts == 1]

    n = np.asarray(outward_normal, dtype=float)
    n /= np.linalg.norm(n)
    o = np.asarray(plane_origin, dtype=float)

    d0 = np.abs((vertices[boundary_edges[:, 0]] - o) @ n)
    d1 = np.abs((vertices[boundary_edges[:, 1]] - o) @ n)
    plane_edges = boundary_edges[(d0 < tol) & (d1 < tol)]

    if len(plane_edges) == 0:
        raise RuntimeError("No boundary loop found on requested port plane.")

    # Build adjacency of boundary edges, then walk each closed loop.
    adjacency = {}
    for a, b in plane_edges:
        a = int(a)
        b = int(b)
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    unused = {tuple(sorted((int(a), int(b)))) for a, b in plane_edges}
    loops = []

    while unused:
        first_edge = next(iter(unused))
        start_v, next_v = first_edge
        loop = [start_v]
        prev = None
        cur = start_v

        while True:
            candidates = adjacency[cur]
            chosen = None

            for nb in candidates:
                e = tuple(sorted((cur, nb)))
                if e in unused and nb != prev:
                    chosen = nb
                    break

            if chosen is None:
                # At the final vertex, the only unused edge may close to start.
                for nb in candidates:
                    e = tuple(sorted((cur, nb)))
                    if e in unused:
                        chosen = nb
                        break

            if chosen is None:
                raise RuntimeError("Open/non-manifold boundary encountered while capping port.")

            e = tuple(sorted((cur, chosen)))
            unused.remove(e)

            prev, cur = cur, chosen

            if cur == start_v:
                break

            loop.append(cur)

        if len(loop) >= 3:
            loops.append(loop)

    new_vertices = vertices.tolist()
    new_faces = faces.tolist()

    for loop in loops:
        pts = vertices[np.asarray(loop, dtype=int)]
        center = pts.mean(axis=0)

        # Project the centroid back to the plane exactly (roundoff protection).
        center = center - n * np.dot(center - o, n)

        center_idx = len(new_vertices)
        new_vertices.append(center.tolist())

        # Determine loop orientation from first non-degenerate fan triangle.
        sign = 0.0
        for j in range(len(loop)):
            a = vertices[loop[j]] - center
            b = vertices[loop[(j + 1) % len(loop)]] - center
            sign = float(np.dot(np.cross(a, b), n))
            if abs(sign) > 1.0e-14:
                break

        forward = sign > 0.0

        for j in range(len(loop)):
            v0 = loop[j]
            v1 = loop[(j + 1) % len(loop)]

            if forward:
                new_faces.append([center_idx, v0, v1])
            else:
                new_faces.append([center_idx, v1, v0])

    capped = trimesh.Trimesh(
        vertices=np.asarray(new_vertices, dtype=float),
        faces=np.asarray(new_faces, dtype=np.int64),
        process=True
    )
    capped.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(capped)

    return capped


def cut_and_cap_planar_ports(mesh, ports):
    """
    Remove each rounded marching-cubes terminal and replace it with an exactly
    planar cap.  No optional triangulation package is required.
    """
    print("Cutting and capping planar ports...")

    result = mesh

    for i, port in enumerate(ports, start=1):
        outward = np.asarray(port["normal"], dtype=float)
        origin = np.asarray(port["point"], dtype=float)

        # slice_mesh_plane keeps the positive side of plane_normal.  The tube
        # interior is opposite to the outward port direction.
        result = trimesh.intersections.slice_mesh_plane(
            result,
            plane_normal=-outward,
            plane_origin=origin,
            cap=False
        )

        result.remove_unreferenced_vertices()
        result = _cap_cut_boundary(
            result,
            plane_origin=origin,
            outward_normal=outward,
            tol=max(1.0e-7, GRID_SPACING * 1.0e-5)
        )

        print(
            f"  port {i}: {port['plane']} = {port['coordinate']:.3f} mm"
        )

    return result


def save_preview(mesh, sampled_branches):
    """
    Save a lightweight 3D preview showing centerlines and sampled STL surface.
    """
    print(f"Writing preview: {OUTPUT_PREVIEW}")

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection="3d")

    for centerline, radius in sampled_branches:
        ax.plot(
            centerline[:, 0],
            centerline[:, 1],
            centerline[:, 2],
            linewidth=2
        )

    stride = max(
        1,
        len(mesh.vertices) // 5000
    )

    surface_points = mesh.vertices[::stride]

    ax.scatter(
        surface_points[:, 0],
        surface_points[:, 1],
        surface_points[:, 2],
        s=0.4,
        alpha=0.18
    )

    bounds = mesh.bounds
    size = bounds[1] - bounds[0]

    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")
    ax.set_title("Smooth irregular branched channel")
    ax.set_box_aspect(size)

    plt.tight_layout()
    plt.savefig(
        OUTPUT_PREVIEW,
        dpi=180
    )
    plt.close(fig)


def main():
    print("===============================================")
    print(" Smooth irregular branched channel generator")
    print("===============================================")

    sampled_branches, ports = prepare_branches()

    mins, maxs = compute_bounds(sampled_branches)

    print("Domain bounds used for implicit grid:")
    print(f"  min = {mins}")
    print(f"  max = {maxs}")

    field = build_implicit_field(
        sampled_branches,
        ports,
        mins,
        maxs
    )

    mesh = extract_surface(
        field,
        mins
    )

    print(f"Writing STL: {OUTPUT_STL}")
    mesh.export(OUTPUT_STL)

    save_preview(
        mesh,
        sampled_branches
    )

    print("")
    print("Finished.")
    print(f"Vertices   : {len(mesh.vertices):,}")
    print(f"Triangles  : {len(mesh.faces):,}")
    print(f"Watertight : {mesh.is_watertight}")
    print(f"Volume     : {mesh.volume:.3f} mm^3")
    print("")
    print("Mesh bounds [mm]:")
    print(mesh.bounds)
    print("")
    print("Generated files:")
    print(f"  {os.path.abspath(OUTPUT_STL)}")
    print(f"  {os.path.abspath(OUTPUT_PREVIEW)}")


if __name__ == "__main__":
    main()
