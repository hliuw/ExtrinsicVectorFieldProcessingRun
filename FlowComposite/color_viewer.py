import sys
import numpy as np
import polyscope as ps
from plyfile import PlyData


def load_ply(path):
    ply = PlyData.read(path)
    v = ply["vertex"]
    vertices = np.stack([v["x"], v["y"], v["z"]], axis=-1).astype(np.float64)
    colors = np.stack([v["red"], v["green"], v["blue"]], axis=-1).astype(np.float64) / 255.0
    faces = np.vstack(list(ply["face"]["vertex_indices"]))
    return vertices, faces, colors


def main():
    paths = sys.argv[1:]
    if not paths:
        print("Usage: python color_viewer.py file1.ply [file2.ply ...]")
        sys.exit(1)

    ps.init()

    for path in paths:
        name = path.rsplit("/", 1)[-1].rsplit("\\", 1)[-1].replace(".ply", "")
        vertices, faces, colors = load_ply(path)
        mesh = ps.register_surface_mesh(name, vertices, faces)
        mesh.add_color_quantity("color", colors, defined_on="vertices", enabled=True)

    ps.show()


if __name__ == "__main__":
    main()
