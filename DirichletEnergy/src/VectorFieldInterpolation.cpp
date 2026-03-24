#include "VectorFieldInterpolation.h"

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace dirichlet::interpolation {
namespace {

std::string LowercaseCopy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool InterpolateByVertexAverage(const TriangleMesh& mesh,
                                const std::vector<Vec3>& vertex_vectors,
                                std::vector<Vec3>& face_vectors,
                                std::string& error_message) {
    if (vertex_vectors.size() != mesh.vertices.size()) {
        error_message = "Vertex vector count does not match the number of mesh vertices.";
        return false;
    }

    face_vectors.clear();
    face_vectors.reserve(mesh.faces.size());
    for (const Triangle& face : mesh.faces) {
        for (std::size_t vertex_index : face.vertex_indices) {
            if (vertex_index >= vertex_vectors.size()) {
                error_message = "Triangle references an out-of-range vertex vector index during interpolation.";
                face_vectors.clear();
                return false;
            }
        }

        const Vec3 averaged_vector =
            (vertex_vectors[face.vertex_indices[0]] + vertex_vectors[face.vertex_indices[1]] +
             vertex_vectors[face.vertex_indices[2]]) /
            3.0;
        face_vectors.push_back(averaged_vector);
    }

    return true;
}

}  // namespace

bool TryParseVertexToFaceInterpolationMethod(const std::string& method_name,
                                             VertexToFaceInterpolationMethod& method) {
    const std::string normalized = LowercaseCopy(method_name);
    if (normalized == "none") {
        method = VertexToFaceInterpolationMethod::kNone;
        return true;
    }
    if (normalized == "averaging" || normalized == "average" || normalized == "vertex-average" ||
        normalized == "vertex-averaging") {
        method = VertexToFaceInterpolationMethod::kVertexAverage;
        return true;
    }
    return false;
}

bool InterpolateVertexVectorsToFaceVectors(const TriangleMesh& mesh,
                                           const std::vector<Vec3>& vertex_vectors,
                                           VertexToFaceInterpolationMethod method,
                                           std::vector<Vec3>& face_vectors,
                                           std::string& error_message) {
    switch (method) {
        case VertexToFaceInterpolationMethod::kVertexAverage:
            return InterpolateByVertexAverage(mesh, vertex_vectors, face_vectors, error_message);
        case VertexToFaceInterpolationMethod::kNone:
            face_vectors.clear();
            error_message = "No interpolation method was selected.";
            return false;
        default:
            face_vectors.clear();
            error_message = "Unsupported interpolation method.";
            return false;
    }
}

}  // namespace dirichlet::interpolation
