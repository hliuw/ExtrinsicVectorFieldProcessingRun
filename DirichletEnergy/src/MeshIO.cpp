#include "MeshIO.h"

namespace dirichlet::io {

bool ReadPlyMeshWithFaceVectors(const std::string& path, TriangleMesh& mesh, std::string& error_message) {
    (void)path;
    mesh = {};
    error_message =
        "PLY reading is not implemented yet. Add a parser that fills vertices, faces, and one tangent vector per face.";
    return false;
}

}  // namespace dirichlet::io
