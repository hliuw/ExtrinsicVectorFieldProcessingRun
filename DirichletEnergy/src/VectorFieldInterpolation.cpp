#include "VectorFieldInterpolation.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

#include "MathUtils.h"

namespace dirichlet::interpolation {
namespace {

std::string LowercaseCopy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool ValidateVertexVectors(const TriangleMesh& mesh,
                           const std::vector<Vec3>& vertex_vectors,
                           std::string& error_message) {
    if (vertex_vectors.size() != mesh.vertices.size()) {
        error_message = "Vertex vector count does not match the number of mesh vertices.";
        return false;
    }
    return true;
}

bool ComputeCornerAngle(const TriangleMesh& mesh,
                        const Triangle& face,
                        std::size_t local_vertex_index,
                        double& corner_angle,
                        std::string& error_message) {
    const std::size_t i0 = face.vertex_indices[local_vertex_index];
    const std::size_t i1 = face.vertex_indices[(local_vertex_index + 1) % 3];
    const std::size_t i2 = face.vertex_indices[(local_vertex_index + 2) % 3];

    if (i0 >= mesh.vertices.size() || i1 >= mesh.vertices.size() || i2 >= mesh.vertices.size()) {
        error_message = "Triangle references an out-of-range vertex while computing corner angles.";
        return false;
    }

    const Vec3 edge1 = mesh.vertices[i1] - mesh.vertices[i0];
    const Vec3 edge2 = mesh.vertices[i2] - mesh.vertices[i0];
    const double norm1 = Norm(edge1);
    const double norm2 = Norm(edge2);
    if (norm1 < kGeometricEpsilon || norm2 < kGeometricEpsilon) {
        error_message = "Degenerate triangle encountered while computing corner angles.";
        return false;
    }

    const double cosine_angle = Clamp(Dot(edge1, edge2) / (norm1 * norm2), -1.0, 1.0);
    corner_angle = std::acos(cosine_angle);
    return true;
}

bool ComputeAngleWeightedVertexNormals(const TriangleMesh& mesh,
                                       std::vector<Vec3>& vertex_normals,
                                       std::string& error_message) {
    vertex_normals.assign(mesh.vertices.size(), {});

    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const Triangle& face = mesh.faces[face_index];
        const Vec3 face_normal = TriangleNormal(mesh, face_index);
        if (Norm(face_normal) < kGeometricEpsilon) {
            error_message = "Degenerate triangle encountered while computing face normals for rotation interpolation.";
            vertex_normals.clear();
            return false;
        }

        for (std::size_t local_vertex_index = 0; local_vertex_index < 3; ++local_vertex_index) {
            double corner_angle = 0.0;
            if (!ComputeCornerAngle(mesh, face, local_vertex_index, corner_angle, error_message)) {
                vertex_normals.clear();
                return false;
            }

            const std::size_t vertex_index = face.vertex_indices[local_vertex_index];
            vertex_normals[vertex_index] = vertex_normals[vertex_index] + face_normal * corner_angle;
        }
    }

    for (std::size_t vertex_index = 0; vertex_index < vertex_normals.size(); ++vertex_index) {
        const double normal_length = Norm(vertex_normals[vertex_index]);
        if (normal_length < kGeometricEpsilon) {
            error_message = "Could not compute a valid angle-weighted normal for vertex " + std::to_string(vertex_index) + ".";
            vertex_normals.clear();
            return false;
        }
        vertex_normals[vertex_index] = vertex_normals[vertex_index] / normal_length;
    }

    return true;
}

bool InterpolateByVertexAverage(const TriangleMesh& mesh,
                                const std::vector<Vec3>& vertex_vectors,
                                std::vector<Vec3>& face_vectors,
                                std::string& error_message) {
    if (!ValidateVertexVectors(mesh, vertex_vectors, error_message)) {
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

bool InterpolateByRotation(const TriangleMesh& mesh,
                           const std::vector<Vec3>& vertex_vectors,
                           std::vector<Vec3>& face_vectors,
                           std::string& error_message) {
    if (!ValidateVertexVectors(mesh, vertex_vectors, error_message)) {
        return false;
    }

    std::vector<Vec3> vertex_normals;
    if (!ComputeAngleWeightedVertexNormals(mesh, vertex_normals, error_message)) {
        return false;
    }

    face_vectors.clear();
    face_vectors.reserve(mesh.faces.size());
    for (std::size_t face_index = 0; face_index < mesh.faces.size(); ++face_index) {
        const Triangle& face = mesh.faces[face_index];
        const std::array<std::size_t, 3>& indices = face.vertex_indices;
        for (std::size_t vertex_index : indices) {
            if (vertex_index >= vertex_vectors.size()) {
                error_message = "Triangle references an out-of-range vertex vector index during rotation interpolation.";
                face_vectors.clear();
                return false;
            }
        }

        Vec3 phong_normal = (vertex_normals[indices[0]] + vertex_normals[indices[1]] + vertex_normals[indices[2]]) / 3.0;
        if (Norm(phong_normal) < kGeometricEpsilon) {
            phong_normal = TriangleNormal(mesh, face_index);
        } else {
            phong_normal = Normalize(phong_normal);
        }

        Vec3 accumulated_face_vector{};
        for (std::size_t vertex_index : indices) {
            const Vec3 rotated_vector = RotateVectorBetweenNormals(
                vertex_vectors[vertex_index], vertex_normals[vertex_index], phong_normal);
            accumulated_face_vector = accumulated_face_vector + rotated_vector;
        }

        face_vectors.push_back(accumulated_face_vector / 3.0);
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
    if (normalized == "rotation" || normalized == "rotational") {
        method = VertexToFaceInterpolationMethod::kRotation;
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
        case VertexToFaceInterpolationMethod::kRotation:
            return InterpolateByRotation(mesh, vertex_vectors, face_vectors, error_message);
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
