import argparse
import numpy as np
import plyfile
import polyscope as ps

def load_ply(path):
    data = plyfile.PlyData.read(path)
    verts = data["vertex"]
    vertices = np.column_stack([verts["x"], verts["y"], verts["z"]])
    vectors = np.column_stack([verts["vf_0"], verts["vf_1"], verts["vf_2"]])
    scalars = np.array(verts["f"])
    face_data = data["face"]["vertex_indices"]
    faces = np.vstack(face_data).astype(np.int32)
    return vertices, faces, vectors, scalars

def compute_geodesic(vertices, source_id):
    diff = vertices - vertices[source_id]
    return np.linalg.norm(diff, axis=1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", nargs="?", default="mesh.ply")
    parser.add_argument("source_vertex", type=int,
                        help="Vertex id to compute Euclidean geodesic distances from")
    parser.add_argument("--true-scale", action="store_true",
                        help="Display vectors at their original world-space magnitude instead of polyscope's auto-scaled length")
    parser.add_argument("--direction", action="store_true",
                        help="Normalize all vectors to unit length, showing only direction")
    parser.add_argument("--divide-by-scalar", action="store_true",
                        help="Divide each vector by the scalar value at its vertex")
    args = parser.parse_args()

    vertices, faces, vectors, scalars = load_ply(args.path)

    if args.source_vertex < 0 or args.source_vertex >= len(vertices):
        raise ValueError(f"source_vertex {args.source_vertex} is out of range [0, {len(vertices) - 1}]")

    geodesic = compute_geodesic(vertices, args.source_vertex)

    if args.divide_by_scalar:
        safe_scalars = np.where(scalars != 0, scalars, 1.0)
        vectors = vectors / safe_scalars[:, np.newaxis]

    if args.direction:
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        vectors = np.where(norms > 0, vectors / norms, vectors)

    vector_kwargs = {}
    if args.true_scale:
        bbox_diagonal = np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0))
        if bbox_diagonal > 0:
            vector_kwargs["length"] = 1.0 / bbox_diagonal

    mag = np.linalg.norm(vectors, axis=1)

    ps.init()
    mesh = ps.register_surface_mesh("mesh", vertices, faces, smooth_shade=True)
    mesh.add_scalar_quantity("f", scalars, enabled=True)
    mesh.add_scalar_quantity("mag", mag, enabled=False)
    mesh.add_scalar_quantity("geodesic", geodesic, enabled=False)
    mesh.add_vector_quantity("vf", vectors, enabled=True, **vector_kwargs)
    ps.show()

if __name__ == "__main__":
    main()
